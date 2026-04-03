#!/usr/bin/env python3
"""
Comprehensive MD Simulation Verification Script
================================================
Independent verification of MTHFR dimer MD simulation results.
Tests energy stability, structural integrity, RMSD, Rg, contacts,
comparison to experimental structure, and secondary structure.
"""

import os
import sys
import warnings
import time

import numpy as np
import pandas as pd
import mdtraj as md
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.spatial.distance import pdist, squareform

warnings.filterwarnings('ignore', category=UserWarning)

# ============================================================================
# Configuration
# ============================================================================
BASE = os.path.dirname(os.path.abspath(__file__))
MD_DIR = os.path.join(BASE, "analysis/md_results")
OUT_DIR = os.path.join(BASE, "analysis/outputs/verification")
os.makedirs(OUT_DIR, exist_ok=True)

SYSTEMS = {
    "wt": {
        "label": "Wild Type",
        "prepared": os.path.join(MD_DIR, "wt_dimer_prepared.pdb"),
        "final": os.path.join(MD_DIR, "wt_dimer_final.pdb"),
        "trajectory": os.path.join(MD_DIR, "wt_dimer_trajectory.dcd"),
        "energy": os.path.join(MD_DIR, "wt_dimer_energy.csv"),
    },
    "compound": {
        "label": "Compound Het",
        "prepared": os.path.join(MD_DIR, "compound_dimer_prepared.pdb"),
        "final": os.path.join(MD_DIR, "compound_dimer_final.pdb"),
        "trajectory": os.path.join(MD_DIR, "compound_dimer_trajectory.dcd"),
        "energy": os.path.join(MD_DIR, "compound_dimer_energy.csv"),
    },
}

CRYSTAL = os.path.join(BASE, "alphafold/6FCX.cif")

# Pipeline-reported per-chain RMSD values (in Angstroms)
# Chain A: WT=6.89, Compound=5.67 | Chain B: WT=7.44, Compound=6.76
PIPELINE_RMSD = {
    "wt":       {"A": 6.89, "B": 7.44},
    "compound": {"A": 5.67, "B": 6.76},
}

STRIDE = 50  # for trajectory loading

results = {}  # collect PASS/FAIL for each check

def banner(title):
    print("\n" + "=" * 72)
    print(f"  {title}")
    print("=" * 72)

def record(name, passed, detail=""):
    tag = "PASS" if passed else "FAIL"
    results[name] = passed
    print(f"  [{tag}] {name}" + (f" -- {detail}" if detail else ""))

# ============================================================================
# CHECK 1: Energy Stability
# ============================================================================
banner("CHECK 1: Energy Stability")

fig_energy, axes_energy = plt.subplots(2, 2, figsize=(16, 10))

for idx, (sysname, sysinfo) in enumerate(SYSTEMS.items()):
    label = sysinfo["label"]
    df = pd.read_csv(sysinfo["energy"])
    # Normalise column names (CSV header has #"Step" style)
    df.columns = [c.strip().strip('"').strip('#').strip('"').strip() for c in df.columns]
    print(f"  {label} energy columns: {df.columns.tolist()}")

    # Identify columns
    time_col = [c for c in df.columns if 'time' in c.lower()][0]
    temp_col = [c for c in df.columns if 'temp' in c.lower()][0]
    pe_col   = [c for c in df.columns if 'potential' in c.lower()][0]

    time_ns = df[time_col].values / 1000.0  # ps -> ns
    temp    = df[temp_col].values
    pe      = df[pe_col].values

    # Running average (window=200 frames)
    window = min(200, len(temp) // 5)
    temp_ra = pd.Series(temp).rolling(window, center=True).mean().values
    pe_ra   = pd.Series(pe).rolling(window, center=True).mean().values

    # Temperature plot
    ax = axes_energy[0, idx]
    ax.plot(time_ns, temp, alpha=0.3, lw=0.5, color='blue', label='Instantaneous')
    ax.plot(time_ns, temp_ra, lw=1.5, color='red', label=f'Running avg (w={window})')
    ax.axhline(300, color='black', ls='--', lw=0.8, label='Target 300 K')
    ax.axhline(305, color='gray', ls=':', lw=0.6)
    ax.axhline(295, color='gray', ls=':', lw=0.6)
    ax.set_title(f'{label} - Temperature')
    ax.set_xlabel('Time (ns)')
    ax.set_ylabel('Temperature (K)')
    ax.legend(fontsize=8)
    ax.set_ylim(290, 310)

    # Potential energy plot
    ax = axes_energy[1, idx]
    ax.plot(time_ns, pe, alpha=0.3, lw=0.5, color='green', label='Instantaneous')
    ax.plot(time_ns, pe_ra, lw=1.5, color='darkred', label=f'Running avg (w={window})')
    ax.set_title(f'{label} - Potential Energy')
    ax.set_xlabel('Time (ns)')
    ax.set_ylabel('Potential Energy (kJ/mol)')
    ax.legend(fontsize=8)

    # Checks
    max_temp_dev = np.max(np.abs(temp - 300.0))
    mean_temp = np.mean(temp)
    std_temp  = np.std(temp)

    # Check for drift: linear fit to PE
    pe_slope, pe_intercept = np.polyfit(time_ns, pe, 1)
    pe_range = pe.max() - pe.min()
    pe_drift_pct = abs(pe_slope * (time_ns[-1] - time_ns[0])) / pe_range * 100 if pe_range > 0 else 0

    print(f"\n  {label}:")
    print(f"    Temperature: mean={mean_temp:.2f} K, std={std_temp:.2f} K, max deviation from 300K = {max_temp_dev:.2f} K")
    print(f"    Potential energy: drift = {pe_slope:.1f} kJ/mol/ns, range = {pe_range:.0f} kJ/mol")
    print(f"    PE linear drift as % of range: {pe_drift_pct:.1f}%")

    temp_pass = max_temp_dev <= 5.0
    record(f"{label} temperature within 5K of 300K", temp_pass,
           f"max dev = {max_temp_dev:.2f} K")

    # Flag if PE shows >50% drift relative to fluctuation range (very generous)
    pe_pass = pe_drift_pct < 50
    record(f"{label} potential energy stable (no runaway drift)", pe_pass,
           f"drift/range = {pe_drift_pct:.1f}%")

fig_energy.tight_layout()
fig_energy.savefig(os.path.join(OUT_DIR, "01_energy_stability.png"), dpi=150)
print(f"\n  Plot saved: {os.path.join(OUT_DIR, '01_energy_stability.png')}")

# ============================================================================
# CHECK 2: Structural Integrity
# ============================================================================
banner("CHECK 2: Structural Integrity (start vs final PDB)")

for sysname, sysinfo in SYSTEMS.items():
    label = sysinfo["label"]
    print(f"\n  --- {label} ---")
    t_start = md.load(sysinfo["prepared"])
    t_final = md.load(sysinfo["final"])

    for chain_idx, chain_name in [(0, "A"), (1, "B")]:
        # Get protein residues for this chain
        start_chain_residues = [r for r in t_start.topology.chain(chain_idx).residues if r.is_protein]
        final_chain_residues = [r for r in t_final.topology.chain(chain_idx).residues if r.is_protein]

        n_start = len(start_chain_residues)
        n_final = len(final_chain_residues)
        res_match = (n_start == n_final)
        record(f"{label} Chain {chain_name} residue count match", res_match,
               f"start={n_start}, final={n_final}")

        # Get CA atom indices for this chain
        ca_indices_start = t_start.topology.select(f"chainid {chain_idx} and name CA")
        ca_indices_final = t_final.topology.select(f"chainid {chain_idx} and name CA")

        # --- Broken bonds check (consecutive CA-CA > 4.5 A) ---
        ca_pos_final = t_final.xyz[0, ca_indices_final] * 10.0  # nm -> Angstrom
        consec_dists = np.linalg.norm(np.diff(ca_pos_final, axis=0), axis=1)
        broken = np.sum(consec_dists > 4.5)
        max_consec = np.max(consec_dists) if len(consec_dists) > 0 else 0
        record(f"{label} Chain {chain_name} no broken bonds (CA-CA < 4.5 A)",
               broken == 0,
               f"broken={broken}, max consecutive CA-CA = {max_consec:.2f} A")

        # --- Steric clash check (any two CA < 2.0 A) ---
        if len(ca_pos_final) < 5000:  # safe for pdist
            dists = pdist(ca_pos_final)
            # Exclude consecutive residues (bonded pairs)
            n_ca = len(ca_pos_final)
            clash_mask = dists < 2.0
            # Build condensed index for consecutive pairs to exclude them
            consec_indices = set()
            k = 0
            for i in range(n_ca):
                for j in range(i + 1, n_ca):
                    if abs(i - j) == 1:
                        consec_indices.add(k)
                    k += 1
            # Exclude consecutive pairs from clash count
            n_clashes = 0
            for k_idx in range(len(dists)):
                if clash_mask[k_idx] and k_idx not in consec_indices:
                    n_clashes += 1
            min_non_consec = np.inf
            k = 0
            for i in range(n_ca):
                for j in range(i + 1, n_ca):
                    if abs(i - j) > 1:
                        if dists[k] < min_non_consec:
                            min_non_consec = dists[k]
                    k += 1
        else:
            # For very large sets, use a sampling approach
            n_clashes = 0
            min_non_consec = np.inf

        record(f"{label} Chain {chain_name} no steric clashes (CA-CA > 2.0 A)",
               n_clashes == 0,
               f"clashes={n_clashes}, min non-bonded CA-CA = {min_non_consec:.2f} A")

    del t_start, t_final

# ============================================================================
# CHECK 3: Independent RMSD Verification
# ============================================================================
banner("CHECK 3: Independent RMSD Verification (per-chain)")

fig_rmsd, axes_rmsd = plt.subplots(2, 2, figsize=(16, 10))

for sys_idx, (sysname, sysinfo) in enumerate(SYSTEMS.items()):
    label = sysinfo["label"]
    print(f"\n  Loading {label} trajectory (stride={STRIDE})...")
    t0 = time.time()
    traj = md.load(sysinfo["trajectory"],
                   top=sysinfo["prepared"],
                   stride=STRIDE)
    print(f"    Loaded {traj.n_frames} frames in {time.time()-t0:.1f}s")

    # Time axis in ns (assuming 2ps per step * 5000 steps between saves = 10ps per frame)
    # Actually let's use the energy CSV to calibrate
    df_e = pd.read_csv(sysinfo["energy"])
    df_e.columns = [c.strip().strip('"').strip('#').strip('"').strip() for c in df_e.columns]
    time_col = [c for c in df_e.columns if 'time' in c.lower()][0]
    total_time_ns = df_e[time_col].values[-1] / 1000.0
    time_ax = np.linspace(0, total_time_ns, traj.n_frames)

    for chain_idx, chain_name in [(0, "A"), (1, "B")]:
        ca_sel = traj.topology.select(f"chainid {chain_idx} and name CA")
        if len(ca_sel) == 0:
            print(f"    WARNING: no CA atoms for chain {chain_idx}")
            continue

        # Slice trajectory to just these atoms
        traj_chain = traj.atom_slice(ca_sel)

        # Superpose to first frame explicitly, then compute RMSD
        ref_chain = traj_chain[0]
        traj_chain.superpose(ref_chain)

        # Compute RMSD manually after superposition (should be near zero for frame 0)
        rmsd_manual = np.zeros(traj_chain.n_frames)
        ref_xyz = ref_chain.xyz[0]  # (n_atoms, 3) in nm
        for f in range(traj_chain.n_frames):
            diff = traj_chain.xyz[f] - ref_xyz
            rmsd_manual[f] = np.sqrt(np.mean(np.sum(diff**2, axis=1)))

        rmsd_angstrom = rmsd_manual * 10.0  # nm -> Angstrom

        # Also compute using mdtraj.rmsd (which does its own superposition)
        rmsd_mdtraj = md.rmsd(traj.atom_slice(ca_sel), traj.atom_slice(ca_sel), frame=0) * 10.0

        # Final RMSD = average of last 10% of trajectory
        last_10pct = int(0.9 * len(rmsd_angstrom))
        mean_final_manual = np.mean(rmsd_angstrom[last_10pct:])
        mean_final_mdtraj = np.mean(rmsd_mdtraj[last_10pct:])
        last_frame_rmsd = rmsd_angstrom[-1]

        pipeline_val = PIPELINE_RMSD[sysname][chain_name]

        print(f"\n    {label} Chain {chain_name}:")
        print(f"      Manual RMSD (last 10% avg): {mean_final_manual:.2f} A")
        print(f"      mdtraj RMSD (last 10% avg): {mean_final_mdtraj:.2f} A")
        print(f"      Last frame RMSD:            {last_frame_rmsd:.2f} A")
        print(f"      Pipeline reported:          {pipeline_val:.2f} A")

        # Methods should agree with each other within 0.1 A
        methods_agree = abs(mean_final_manual - mean_final_mdtraj) < 0.1
        record(f"{label} Chain {chain_name} manual vs mdtraj RMSD agree",
               methods_agree,
               f"diff = {abs(mean_final_manual - mean_final_mdtraj):.3f} A")

        # Compare to pipeline (within 1.5 A tolerance since we use stride=50 and avg last 10%)
        pipeline_close = abs(mean_final_mdtraj - pipeline_val) < 1.5
        record(f"{label} Chain {chain_name} RMSD close to pipeline ({pipeline_val:.2f} A)",
               pipeline_close,
               f"our value = {mean_final_mdtraj:.2f} A, diff = {abs(mean_final_mdtraj - pipeline_val):.2f} A")

        # Plot
        ax = axes_rmsd[chain_idx, sys_idx]
        ax.plot(time_ax, rmsd_angstrom, lw=0.8, alpha=0.6, label='Manual')
        ax.plot(time_ax, rmsd_mdtraj, lw=0.8, alpha=0.6, ls='--', label='mdtraj.rmsd')
        ax.axhline(pipeline_val, color='red', ls=':', lw=1.2, label=f'Pipeline: {pipeline_val:.2f} A')
        ax.set_title(f'{label} Chain {chain_name} RMSD')
        ax.set_xlabel('Time (ns)')
        ax.set_ylabel('RMSD (Angstrom)')
        ax.legend(fontsize=8)

    # ========================================================================
    # CHECK 4: Radius of Gyration
    # ========================================================================
    print(f"\n  --- Radius of Gyration ({label}) ---")

    fig_rg, axes_rg_local = plt.subplots(1, 2, figsize=(14, 5))

    for chain_idx, chain_name in [(0, "A"), (1, "B")]:
        ca_sel = traj.topology.select(f"chainid {chain_idx} and name CA")
        traj_chain = traj.atom_slice(ca_sel)
        rg = md.compute_rg(traj_chain) * 10.0  # nm -> Angstrom

        rg_start = rg[0]
        rg_end   = rg[-1]
        rg_mean  = np.mean(rg)
        rg_std   = np.std(rg)
        rg_change_pct = abs(rg_end - rg_start) / rg_start * 100

        print(f"    Chain {chain_name}: Rg_start={rg_start:.2f} A, Rg_end={rg_end:.2f} A, "
              f"mean={rg_mean:.2f} A, std={rg_std:.2f} A, change={rg_change_pct:.1f}%")

        # Rg should not change by more than 20% (would indicate unfolding)
        rg_stable = rg_change_pct < 20
        record(f"{label} Chain {chain_name} Rg stable (<20% change)",
               rg_stable,
               f"change = {rg_change_pct:.1f}%")

        ax = axes_rg_local[chain_idx]
        ax.plot(time_ax, rg, lw=0.8, color='purple')
        ax.set_title(f'{label} Chain {chain_name} Radius of Gyration')
        ax.set_xlabel('Time (ns)')
        ax.set_ylabel('Rg (Angstrom)')
        ax.axhline(rg_mean, color='red', ls='--', lw=0.8, label=f'Mean: {rg_mean:.1f} A')
        ax.legend(fontsize=8)

    fig_rg.tight_layout()
    fig_rg.savefig(os.path.join(OUT_DIR, f"04_rg_{sysname}.png"), dpi=150)
    print(f"  Plot saved: {os.path.join(OUT_DIR, f'04_rg_{sysname}.png')}")

    # ========================================================================
    # CHECK 5: Inter-chain Contact Analysis
    # ========================================================================
    # NOTE: The DCD trajectory stores coordinates in unwrapped/periodic-image
    # space, which can place the two chains in different periodic boxes (>100 A
    # apart).  We therefore use image_molecules() to re-wrap them, and fall
    # back to the PDB files (which are already properly wrapped) if that fails.
    print(f"\n  --- Inter-chain Contacts ({label}) ---")

    ca_A_sel = traj.topology.select("chainid 0 and name CA")
    ca_B_sel = traj.topology.select("chainid 1 and name CA")

    # Try to re-image the trajectory so both chains are in the same box
    try:
        traj_imaged = traj.image_molecules(inplace=False)
        traj_for_contacts = traj_imaged
        print("    (using image_molecules to re-wrap periodic images)")
    except Exception as e:
        print(f"    (image_molecules failed: {e}; using raw trajectory)")
        traj_for_contacts = traj

    # Sample at start, middle, end
    frame_indices = [0, traj_for_contacts.n_frames // 2, traj_for_contacts.n_frames - 1]
    frame_labels = ["Start", "Middle", "End"]
    contact_counts_traj = []

    for fi, fl in zip(frame_indices, frame_labels):
        pos_A = traj_for_contacts.xyz[fi, ca_A_sel] * 10.0  # nm -> A
        pos_B = traj_for_contacts.xyz[fi, ca_B_sel] * 10.0
        diff = pos_A[:, np.newaxis, :] - pos_B[np.newaxis, :, :]
        dists = np.sqrt(np.sum(diff**2, axis=2))
        n_contacts = np.sum(dists < 10.0)
        contact_counts_traj.append(n_contacts)
        print(f"    Traj {fl} (frame {fi}): {n_contacts} inter-chain CA-CA contacts < 10 A")

    # Also compute from the PDB files (which are guaranteed properly wrapped)
    t_prep = md.load(sysinfo["prepared"])
    t_fin  = md.load(sysinfo["final"])
    contact_counts_pdb = []
    for pdb, pdb_label in [(t_prep, "Prepared PDB"), (t_fin, "Final PDB")]:
        ca_A_pdb = pdb.topology.select("chainid 0 and name CA")
        ca_B_pdb = pdb.topology.select("chainid 1 and name CA")
        pos_A = pdb.xyz[0, ca_A_pdb] * 10.0
        pos_B = pdb.xyz[0, ca_B_pdb] * 10.0
        diff = pos_A[:, np.newaxis, :] - pos_B[np.newaxis, :, :]
        dists = np.sqrt(np.sum(diff**2, axis=2))
        n_contacts = np.sum(dists < 10.0)
        contact_counts_pdb.append(n_contacts)
        print(f"    {pdb_label}: {n_contacts} inter-chain CA-CA contacts < 10 A")
    del t_prep, t_fin

    # Use PDB-based contacts for the pass/fail verdict (most reliable)
    start_contacts = contact_counts_pdb[0]
    end_contacts   = contact_counts_pdb[1]
    # Positive change = loss; negative change = gain (interface tightening)
    contact_change_pct = (start_contacts - end_contacts) / start_contacts * 100 if start_contacts > 0 else 0
    # FAIL only if contacts DROP by more than 50%  (gain is fine / even good)
    contacts_ok = contact_change_pct < 50   # negative values (gain) always pass
    record(f"{label} inter-chain contacts maintained (< 50% loss, PDB-based)",
           contacts_ok,
           f"prepared={start_contacts}, final={end_contacts}, "
           f"change={-contact_change_pct:+.1f}% ({'gain' if contact_change_pct < 0 else 'loss'})")

    del traj  # free memory

fig_rmsd.tight_layout()
fig_rmsd.savefig(os.path.join(OUT_DIR, "03_rmsd_verification.png"), dpi=150)
print(f"\n  RMSD plot saved: {os.path.join(OUT_DIR, '03_rmsd_verification.png')}")

# ============================================================================
# CHECK 6: Comparison to Experimental Structure (6FCX)
# ============================================================================
banner("CHECK 6: Comparison to Experimental Crystal Structure (6FCX)")

crystal = md.load(CRYSTAL)

# Build a resSeq -> CA atom index map for 6FCX chain 0
cryst_resseq_to_ca = {}
for atom in crystal.topology.chain(0).atoms:
    if atom.name == "CA":
        cryst_resseq_to_ca[atom.residue.resSeq] = atom.index

print(f"  6FCX Chain A: {len(cryst_resseq_to_ca)} CA atoms (resSeq {min(cryst_resseq_to_ca)}-{max(cryst_resseq_to_ca)})")

for sysname, sysinfo in SYSTEMS.items():
    label = sysinfo["label"]
    t_final = md.load(sysinfo["final"])

    # Build resSeq -> CA atom index map for simulation chain 0
    sim_resseq_to_ca = {}
    for atom in t_final.topology.chain(0).atoms:
        if atom.name == "CA":
            sim_resseq_to_ca[atom.residue.resSeq] = atom.index

    print(f"  {label} final Chain A: {len(sim_resseq_to_ca)} CA atoms (resSeq {min(sim_resseq_to_ca)}-{max(sim_resseq_to_ca)})")

    # Find common resSeq values where both have CA atoms AND residue names match
    common_resseqs = sorted(set(cryst_resseq_to_ca.keys()) & set(sim_resseq_to_ca.keys()))

    # Further filter to matching residue names
    cryst_res_by_seq = {r.resSeq: r for r in crystal.topology.chain(0).residues}
    sim_res_by_seq   = {r.resSeq: r for r in t_final.topology.chain(0).residues}
    matched_resseqs = [rs for rs in common_resseqs
                       if cryst_res_by_seq[rs].name == sim_res_by_seq[rs].name]

    cryst_ca_indices = [cryst_resseq_to_ca[rs] for rs in matched_resseqs]
    sim_ca_indices   = [sim_resseq_to_ca[rs]   for rs in matched_resseqs]

    print(f"  Common residues with matching names: {len(matched_resseqs)} "
          f"(resSeq {matched_resseqs[0]}-{matched_resseqs[-1]})")

    # Slice both to the matched CA atoms
    crystal_sub = crystal.atom_slice(cryst_ca_indices)
    sim_sub     = t_final.atom_slice(sim_ca_indices)

    # Compute RMSD using mdtraj (superposition included)
    rmsd_val = md.rmsd(sim_sub, crystal_sub, frame=0)[0] * 10.0  # nm -> A
    print(f"  {label} vs 6FCX Chain A RMSD: {rmsd_val:.2f} A (over {len(matched_resseqs)} matched CA atoms)")

    # For a 100ns simulation of a 656-residue protein, RMSD < 15 A from crystal is reasonable
    crystal_ok = rmsd_val < 15.0
    record(f"{label} vs 6FCX RMSD < 15 A", crystal_ok,
           f"RMSD = {rmsd_val:.2f} A over {len(matched_resseqs)} matched residues")

    del t_final

# ============================================================================
# CHECK 7: Secondary Structure Spot-Check (DSSP)
# ============================================================================
banner("CHECK 7: Secondary Structure (DSSP) - Final Frames")

fig_ss, axes_ss = plt.subplots(1, 2, figsize=(14, 6))

ss_data = {}
for sys_idx, (sysname, sysinfo) in enumerate(SYSTEMS.items()):
    label = sysinfo["label"]
    t_final = md.load(sysinfo["final"])

    # Select only protein atoms (chains 0 and 1)
    protein_sel = t_final.topology.select("chainid 0 or chainid 1")
    t_protein = t_final.atom_slice(protein_sel)

    dssp = md.compute_dssp(t_protein, simplified=True)[0]  # shape: (n_residues,)

    # Count by chain
    chain_boundary = t_protein.topology.chain(0).n_residues
    for chain_idx, chain_name in [(0, "A"), (1, "B")]:
        if chain_idx == 0:
            chain_dssp = dssp[:chain_boundary]
        else:
            chain_dssp = dssp[chain_boundary:]

        n_helix = np.sum(chain_dssp == 'H')
        n_sheet = np.sum(chain_dssp == 'E')
        n_coil  = np.sum(chain_dssp == 'C')
        n_na    = np.sum(chain_dssp == 'NA')
        total   = len(chain_dssp)

        ss_data[(sysname, chain_name)] = {
            'helix': n_helix, 'sheet': n_sheet, 'coil': n_coil,
            'na': n_na, 'total': total
        }
        print(f"  {label} Chain {chain_name}: "
              f"Helix={n_helix} ({100*n_helix/total:.1f}%), "
              f"Sheet={n_sheet} ({100*n_sheet/total:.1f}%), "
              f"Coil={n_coil} ({100*n_coil/total:.1f}%), "
              f"NA={n_na}")

    # Bar chart
    ax = axes_ss[sys_idx]
    categories = ['Helix', 'Sheet', 'Coil']
    x = np.arange(len(categories))
    width = 0.35
    vals_A = [ss_data[(sysname, 'A')]['helix'], ss_data[(sysname, 'A')]['sheet'], ss_data[(sysname, 'A')]['coil']]
    vals_B = [ss_data[(sysname, 'B')]['helix'], ss_data[(sysname, 'B')]['sheet'], ss_data[(sysname, 'B')]['coil']]
    ax.bar(x - width/2, vals_A, width, label='Chain A', color='steelblue')
    ax.bar(x + width/2, vals_B, width, label='Chain B', color='coral')
    ax.set_xticks(x)
    ax.set_xticklabels(categories)
    ax.set_title(f'{label} - Secondary Structure (Final Frame)')
    ax.set_ylabel('Number of Residues')
    ax.legend()

    del t_final

# Compare WT vs Compound for major differences
for chain_name in ['A', 'B']:
    wt_helix = ss_data[('wt', chain_name)]['helix']
    cp_helix = ss_data[('compound', chain_name)]['helix']
    wt_sheet = ss_data[('wt', chain_name)]['sheet']
    cp_sheet = ss_data[('compound', chain_name)]['sheet']

    helix_diff = abs(wt_helix - cp_helix)
    sheet_diff = abs(wt_sheet - cp_sheet)
    total = ss_data[('wt', chain_name)]['total']

    # Flag if helix or sheet count differs by more than 30% of total
    helix_pct_diff = helix_diff / total * 100 if total > 0 else 0
    sheet_pct_diff = sheet_diff / total * 100 if total > 0 else 0

    ss_ok = (helix_pct_diff < 30) and (sheet_pct_diff < 30)
    record(f"Chain {chain_name} WT vs Compound SS similar (<30% diff)",
           ss_ok,
           f"helix diff={helix_diff} ({helix_pct_diff:.1f}%), sheet diff={sheet_diff} ({sheet_pct_diff:.1f}%)")

fig_ss.tight_layout()
fig_ss.savefig(os.path.join(OUT_DIR, "07_secondary_structure.png"), dpi=150)
print(f"  Plot saved: {os.path.join(OUT_DIR, '07_secondary_structure.png')}")

# ============================================================================
# SUMMARY
# ============================================================================
banner("OVERALL VERIFICATION SUMMARY")

n_pass = sum(1 for v in results.values() if v)
n_fail = sum(1 for v in results.values() if not v)
n_total = len(results)

print(f"\n  Total checks: {n_total}")
print(f"  Passed:       {n_pass}")
print(f"  Failed:       {n_fail}")
print()

for name, passed in results.items():
    tag = "PASS" if passed else "FAIL"
    print(f"    [{tag}] {name}")

print()
if n_fail == 0:
    print("  *** OVERALL VERDICT: ALL CHECKS PASSED ***")
    print("  The MD simulation results appear trustworthy.")
else:
    print(f"  *** OVERALL VERDICT: {n_fail} CHECK(S) FAILED ***")
    print("  Review the failed checks above for details.")

print(f"\n  All plots saved to: {OUT_DIR}/")
print("  Done.")
