#!/usr/bin/env python3
"""
MTHFR Molecular Dynamics Simulation Pipeline
==============================================
Runs OpenMM MD simulations on WT and compound heterozygous MTHFR dimers.
Compares structural dynamics, FAD stability, and dimer interface behavior.

Requirements: conda activate md (OpenMM, PDBFixer, MDTraj, matplotlib, numpy)
Hardware: Apple M1 Max with OpenCL acceleration

Usage:
    python run_md.py                    # Run both simulations (default 10ns)
    python run_md.py --length 100       # Run 100ns simulations
    python run_md.py --analyze-only     # Just analyze existing trajectories

Author: Igor Mihaljko / DSM.Promo | License: CC BY-NC-SA 4.0
"""
import argparse, os, sys, time
from pathlib import Path

try:
    import openmm
    from openmm import app, unit
    from pdbfixer import PDBFixer
    import mdtraj as md
    import numpy as np
    import matplotlib; matplotlib.use('Agg')
    import matplotlib.pyplot as plt
except ImportError as e:
    print(f"Missing dependency: {e}")
    print("Run: conda activate md")
    sys.exit(1)

RESULTS = Path("alphafold/results_all")
MD_OUTPUT = Path("analysis/md_results")
FIGURES = Path("analysis/outputs/figures")

SIMULATIONS = {
    "wt_dimer": {
        "cif": "wt_dimer_run1",
        "label": "Wild-Type Dimer",
        "color": "#2196F3",
    },
    "compound_dimer": {
        "cif": "compound_dimer_run1",
        "label": "Compound Het Dimer",
        "color": "#F44336",
    },
}

def find_best_cif(job_dir):
    """Find the rank-0 (best) CIF model."""
    cifs = sorted(Path(job_dir).rglob("*model_0.cif"))
    return str(cifs[0]) if cifs else None

def prepare_system(cif_path, name):
    """Prepare protein system for MD simulation."""
    print(f"  Preparing {name}...")

    # Fix the structure (add missing atoms, hydrogens, solvent)
    fixer = PDBFixer(filename=cif_path)
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.4)  # pH 7.4

    # Add solvent (water box with 1.0 nm padding)
    fixer.addSolvent(padding=1.0 * unit.nanometers, ionicStrength=0.15 * unit.molar)

    # Save prepared structure
    prep_path = MD_OUTPUT / f"{name}_prepared.pdb"
    with open(prep_path, 'w') as f:
        app.PDBFile.writeFile(fixer.topology, fixer.positions, f)

    print(f"  Saved prepared structure: {prep_path}")
    print(f"  Atoms: {fixer.topology.getNumAtoms()}")

    return fixer.topology, fixer.positions

def run_simulation(topology, positions, name, length_ns=10, platform_name="OpenCL"):
    """Run MD simulation with OpenMM."""
    print(f"\n  Running {length_ns}ns MD for {name} on {platform_name}...")

    # Force field
    forcefield = app.ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')

    # Create system
    system = forcefield.createSystem(
        topology,
        nonbondedMethod=app.PME,
        nonbondedCutoff=1.0 * unit.nanometers,
        constraints=app.HBonds,
    )

    # Integrator (2fs timestep)
    integrator = openmm.LangevinMiddleIntegrator(
        300 * unit.kelvin,     # Temperature
        1.0 / unit.picoseconds, # Friction
        0.002 * unit.picoseconds # Timestep
    )

    # Platform
    try:
        platform = openmm.Platform.getPlatformByName(platform_name)
    except Exception:
        platform = openmm.Platform.getPlatformByName("CPU")
        print(f"  Falling back to CPU platform")

    # Create simulation
    simulation = app.Simulation(topology, system, integrator, platform)
    simulation.context.setPositions(positions)

    # Energy minimization
    print(f"  Energy minimization...")
    simulation.minimizeEnergy(maxIterations=1000)

    # Equilibration (100ps NVT)
    print(f"  Equilibrating (100ps)...")
    simulation.context.setVelocitiesToTemperature(300 * unit.kelvin)
    simulation.step(50000)  # 100ps at 2fs

    # Production MD
    traj_path = MD_OUTPUT / f"{name}_trajectory.dcd"
    log_path = MD_OUTPUT / f"{name}_energy.csv"

    n_steps = int(length_ns * 500000)  # ns to steps at 2fs
    report_interval = 5000  # Save every 10ps

    simulation.reporters.append(app.DCDReporter(str(traj_path), report_interval))
    simulation.reporters.append(app.StateDataReporter(
        str(log_path), report_interval,
        step=True, time=True, potentialEnergy=True,
        temperature=True, speed=True
    ))

    print(f"  Production MD: {length_ns}ns ({n_steps} steps)...")
    start_time = time.time()

    # Run in chunks to show progress
    chunk_size = 50000  # 100ps chunks
    n_chunks = n_steps // chunk_size
    for i in range(n_chunks):
        simulation.step(chunk_size)
        elapsed = time.time() - start_time
        progress = (i + 1) / n_chunks * 100
        speed = (i + 1) * chunk_size * 0.002 / elapsed  # ps/s
        eta = (n_chunks - i - 1) * elapsed / (i + 1)
        print(f"    {progress:.0f}% ({(i+1)*0.1:.1f}ns / {length_ns}ns) "
              f"speed={speed:.1f} ps/s, ETA={eta/60:.0f}min", end='\r')

    remaining = n_steps % chunk_size
    if remaining:
        simulation.step(remaining)

    total_time = time.time() - start_time
    print(f"\n  Completed in {total_time/60:.1f} minutes")
    print(f"  Trajectory: {traj_path}")

    # Save final state
    state = simulation.context.getState(getPositions=True)
    final_path = MD_OUTPUT / f"{name}_final.pdb"
    with open(final_path, 'w') as f:
        app.PDBFile.writeFile(topology, state.getPositions(), f)

    return traj_path

def analyze_trajectory(name, label, color):
    """Analyze MD trajectory and generate plots."""
    print(f"\n  Analyzing {name}...")

    traj_path = MD_OUTPUT / f"{name}_trajectory.dcd"
    top_path = MD_OUTPUT / f"{name}_prepared.pdb"

    if not traj_path.exists():
        print(f"  No trajectory found for {name}")
        return None

    # Load trajectory
    traj = md.load(str(traj_path), top=str(top_path))
    print(f"  Loaded {traj.n_frames} frames, {traj.n_atoms} atoms")

    # Select protein atoms only
    protein = traj.atom_slice(traj.topology.select('protein'))

    # RMSD
    rmsd = md.rmsd(protein, protein, 0) * 10  # Convert to Angstroms

    # RMSF per residue
    rmsf = md.rmsf(protein, protein, 0) * 10  # Angstroms

    # Get CA indices for per-residue RMSF
    ca_indices = protein.topology.select('name CA')
    ca_rmsf = rmsf[ca_indices] if len(ca_indices) > 0 else rmsf

    results = {
        "name": name,
        "label": label,
        "color": color,
        "rmsd": rmsd,
        "rmsf": ca_rmsf,
        "n_frames": traj.n_frames,
        "time_ns": traj.time[-1] / 1000 if len(traj.time) > 0 else 0,
    }

    # Save per-simulation RMSD
    np.savetxt(MD_OUTPUT / f"{name}_rmsd.csv",
               np.column_stack([np.arange(len(rmsd)) * 0.01, rmsd]),
               header="time_ns,rmsd_angstrom", delimiter=",", comments="")

    return results

def plot_comparison(all_results):
    """Generate comparison plots."""
    print("\n  Generating comparison plots...")

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # 1. RMSD over time
    for r in all_results:
        time_ns = np.arange(len(r['rmsd'])) * 0.01
        axes[0,0].plot(time_ns, r['rmsd'], label=r['label'], color=r['color'], alpha=0.8, linewidth=0.5)
    axes[0,0].set_xlabel('Time (ns)')
    axes[0,0].set_ylabel('RMSD (A)')
    axes[0,0].set_title('Backbone RMSD Over Time')
    axes[0,0].legend()

    # 2. RMSD distribution
    for r in all_results:
        axes[0,1].hist(r['rmsd'], bins=50, alpha=0.5, label=f"{r['label']} (mean={np.mean(r['rmsd']):.2f}A)", color=r['color'])
    axes[0,1].set_xlabel('RMSD (A)')
    axes[0,1].set_ylabel('Count')
    axes[0,1].set_title('RMSD Distribution')
    axes[0,1].legend()

    # 3. Per-residue RMSF
    for r in all_results:
        residues = np.arange(1, len(r['rmsf']) + 1)
        axes[1,0].plot(residues, r['rmsf'], label=r['label'], color=r['color'], alpha=0.8, linewidth=0.8)
    axes[1,0].axvline(x=222, color='red', linestyle='--', alpha=0.5, label='pos 222')
    axes[1,0].axvline(x=429, color='purple', linestyle='--', alpha=0.5, label='pos 429')
    axes[1,0].set_xlabel('Residue Number')
    axes[1,0].set_ylabel('RMSF (A)')
    axes[1,0].set_title('Per-Residue Flexibility (RMSF)')
    axes[1,0].legend(fontsize=8)

    # 4. RMSF difference (Compound - WT)
    if len(all_results) >= 2:
        wt = all_results[0]
        cp = all_results[1]
        n = min(len(wt['rmsf']), len(cp['rmsf']))
        diff = cp['rmsf'][:n] - wt['rmsf'][:n]
        residues = np.arange(1, n + 1)
        colors = ['#F44336' if d > 0 else '#4CAF50' for d in diff]
        axes[1,1].bar(residues, diff, width=1, color=colors, alpha=0.7)
        axes[1,1].axhline(y=0, color='black', linewidth=0.5)
        axes[1,1].axvline(x=222, color='red', linestyle='--', alpha=0.5)
        axes[1,1].axvline(x=429, color='purple', linestyle='--', alpha=0.5)
        axes[1,1].set_xlabel('Residue Number')
        axes[1,1].set_ylabel('RMSF Difference (A)\n(Compound - WT)')
        axes[1,1].set_title('Flexibility Difference (red = compound more flexible)')

    plt.suptitle('Molecular Dynamics Comparison: WT vs Compound Heterozygous MTHFR Dimer', fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(FIGURES / 'md_comparison.png', dpi=200, bbox_inches='tight')
    plt.close()
    print(f"  Saved md_comparison.png")

    # Summary stats
    print("\n" + "=" * 60)
    print("MD SIMULATION SUMMARY")
    print("=" * 60)
    for r in all_results:
        print(f"\n  {r['label']}:")
        print(f"    Frames: {r['n_frames']}")
        print(f"    Simulation time: {r['time_ns']:.1f} ns")
        print(f"    Mean RMSD: {np.mean(r['rmsd']):.2f} +/- {np.std(r['rmsd']):.2f} A")
        print(f"    RMSF at pos 222: {r['rmsf'][221]:.2f} A" if len(r['rmsf']) > 221 else "")
        print(f"    RMSF at pos 429: {r['rmsf'][428]:.2f} A" if len(r['rmsf']) > 428 else "")

def main():
    parser = argparse.ArgumentParser(description='MTHFR MD Simulation Pipeline')
    parser.add_argument('--length', type=float, default=10, help='Simulation length in ns (default: 10)')
    parser.add_argument('--analyze-only', action='store_true', help='Only analyze existing trajectories')
    parser.add_argument('--platform', default='OpenCL', help='OpenMM platform (OpenCL, CPU)')
    args = parser.parse_args()

    print("=" * 60)
    print("MTHFR Molecular Dynamics Pipeline")
    print(f"Simulation length: {args.length} ns")
    print("=" * 60)

    MD_OUTPUT.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)

    if not args.analyze_only:
        for name, config in SIMULATIONS.items():
            cif_dir = RESULTS / config["cif"]
            cif_path = find_best_cif(cif_dir)

            if not cif_path:
                print(f"  ERROR: No CIF found for {name} in {cif_dir}")
                continue

            print(f"\n{'='*60}")
            print(f"Simulation: {config['label']}")
            print(f"{'='*60}")

            topology, positions = prepare_system(cif_path, name)
            run_simulation(topology, positions, name,
                         length_ns=args.length, platform_name=args.platform)

    # Analyze all trajectories
    print(f"\n{'='*60}")
    print("Analysis")
    print(f"{'='*60}")

    all_results = []
    for name, config in SIMULATIONS.items():
        result = analyze_trajectory(name, config["label"], config["color"])
        if result:
            all_results.append(result)

    if len(all_results) >= 2:
        plot_comparison(all_results)

    print("\n" + "=" * 60)
    print("DONE!")
    print("=" * 60)

if __name__ == "__main__":
    main()
