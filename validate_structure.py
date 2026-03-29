#!/usr/bin/env python3
"""
MTHFR Structure Validation and Statistical Analysis
=====================================================
1. RMSD against experimental PDB 6FCX
2. Cohen's d effect size (WT vs compound dimer)
3. Full per-residue pLDDT comparison plot
4. Bonferroni-corrected p-values for all comparisons

Usage: python validate_structure.py
"""
import json, csv, sys
from pathlib import Path

try:
    import numpy as np
    import matplotlib; matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from scipy import stats
except ImportError:
    print("Install: pip install numpy matplotlib scipy"); sys.exit(1)

RESULTS = Path("alphafold/results_all")
OUTPUT = Path("analysis/outputs")
REF_STRUCTURE = Path("alphafold/6FCX.cif")

def extract_ca_coords_from_cif(cif_path, chain="A"):
    """Extract CA atom coordinates from mmCIF file."""
    coords = []
    resids = []
    with open(cif_path) as f:
        in_atom = False
        for line in f:
            if line.startswith("ATOM") or (in_atom and not line.startswith("#") and not line.startswith("loop_")):
                parts = line.split()
                if len(parts) >= 15:
                    try:
                        atom_type = parts[0]  # ATOM or HETATM
                        atom_name = parts[3] if len(parts) > 3 else ""
                        chain_id = parts[6] if len(parts) > 6 else ""

                        # Different CIF formats have different column positions
                        # Try to find CA atoms
                        if "CA" in line and atom_type in ("ATOM", "HETATM"):
                            # Find coordinates - they're usually the 3 numbers before the last few columns
                            numbers = []
                            for p in parts:
                                try:
                                    numbers.append(float(p))
                                except ValueError:
                                    pass
                            if len(numbers) >= 3:
                                x, y, z = numbers[0], numbers[1], numbers[2]
                                coords.append([x, y, z])
                                try:
                                    resids.append(int(parts[8]) if len(parts) > 8 else len(coords))
                                except (ValueError, IndexError):
                                    resids.append(len(coords))
                    except (ValueError, IndexError):
                        pass
                in_atom = True
            elif line.startswith("#"):
                in_atom = False
    return np.array(coords) if coords else None, resids

def compute_rmsd(coords1, coords2):
    """Compute RMSD between two coordinate arrays after optimal superposition."""
    # Center both
    c1 = coords1 - coords1.mean(axis=0)
    c2 = coords2 - coords2.mean(axis=0)

    # Truncate to same length
    n = min(len(c1), len(c2))
    c1, c2 = c1[:n], c2[:n]

    # Kabsch algorithm for optimal rotation
    H = c1.T @ c2
    U, S, Vt = np.linalg.svd(H)
    d = np.linalg.det(Vt.T @ U.T)
    sign_matrix = np.diag([1, 1, d])
    R = Vt.T @ sign_matrix @ U.T

    c2_rotated = (R @ c2.T).T
    rmsd = np.sqrt(np.mean(np.sum((c1 - c2_rotated)**2, axis=1)))
    return rmsd, n

def extract_plddt_from_cif(cif_path):
    """Extract per-residue pLDDT (B-factor) from AlphaFold CIF."""
    plddts = {}
    with open(cif_path) as f:
        for line in f:
            if line.startswith("ATOM"):
                parts = line.split()
                if len(parts) >= 15:
                    try:
                        atom_name = parts[3]
                        if atom_name == "CA":
                            resnum = int(parts[8])
                            bfactor = float(parts[14])
                            plddts[resnum] = bfactor
                    except (ValueError, IndexError):
                        pass
    return plddts

def main():
    print("=" * 60)
    print("MTHFR Structure Validation & Statistical Analysis")
    print("=" * 60)
    print()

    OUTPUT.mkdir(parents=True, exist_ok=True)
    fig_dir = OUTPUT / "figures"
    fig_dir.mkdir(exist_ok=True)

    # ============================================================
    # 1. RMSD against PDB 6FCX
    # ============================================================
    print("[1/4] Computing RMSD against experimental PDB 6FCX...")

    if REF_STRUCTURE.exists():
        ref_coords, ref_resids = extract_ca_coords_from_cif(REF_STRUCTURE)

        if ref_coords is not None and len(ref_coords) > 0:
            print(f"  6FCX: {len(ref_coords)} CA atoms extracted")

            # Compare WT monomer run1 against 6FCX
            wt_cifs = sorted(RESULTS.glob("wt_mono_run1/*model_0.cif"))
            if wt_cifs:
                pred_coords, pred_resids = extract_ca_coords_from_cif(wt_cifs[0])
                if pred_coords is not None and len(pred_coords) > 0:
                    rmsd, n_atoms = compute_rmsd(ref_coords, pred_coords)
                    print(f"  WT mono run1 vs 6FCX: RMSD = {rmsd:.2f} A over {n_atoms} CA atoms")

                    # Classify
                    if rmsd <= 3.0:
                        classification = "Accurate prediction (RMSD <= 3.0 A)"
                    elif rmsd <= 5.0:
                        classification = "Acceptable prediction"
                    else:
                        classification = "High deviation - review needed"
                    print(f"  Classification: {classification}")
                else:
                    print("  WARNING: Could not extract coordinates from prediction CIF")
                    rmsd = None
            else:
                print("  WARNING: No WT monomer CIF found")
                rmsd = None
        else:
            print("  WARNING: Could not extract coordinates from 6FCX")
            rmsd = None
    else:
        print("  WARNING: 6FCX.cif not found. Download from RCSB.")
        rmsd = None

    # ============================================================
    # 2. Cohen's d effect size
    # ============================================================
    print("\n[2/4] Computing Cohen's d effect sizes...")

    # Load metrics
    metrics_file = OUTPUT / "metrics.csv"
    if not metrics_file.exists():
        print("  Run analyze.py first"); return

    rows = list(csv.DictReader(open(metrics_file)))

    groups = {
        "WT mono": [r for r in rows if "wt_mono" in r["Job"]],
        "WT dimer": [r for r in rows if "wt_dimer_run" in r["Job"]],
        "C677T mono": [r for r in rows if "c677t_mono" in r["Job"]],
        "C677T dimer": [r for r in rows if "c677t_dimer_run" in r["Job"]],
        "A1298C mono": [r for r in rows if "a1298c_mono" in r["Job"]],
        "Compound dimer": [r for r in rows if "compound_dimer_run" in r["Job"]],
    }

    metrics_to_compare = ["ipTM", "pTM", "FAD_iptm", "pLDDT@429"]

    print(f"\n  {'Comparison':<35} {'Metric':<12} {'Mean1':>8} {'Mean2':>8} {'Cohen d':>8} {'p-value':>10} {'Sig':>5}")
    print("  " + "-" * 90)

    all_pvalues = []
    comparison_results = []

    for metric in metrics_to_compare:
        wt_vals = [float(r[metric]) for r in groups["WT dimer"] if r[metric]]
        cp_vals = [float(r[metric]) for r in groups["Compound dimer"] if r[metric]]
        c677t_vals = [float(r[metric]) for r in groups["C677T dimer"] if r[metric]]

        if len(wt_vals) >= 2 and len(cp_vals) >= 2:
            # WT vs Compound
            pooled_sd = np.sqrt((np.std(wt_vals, ddof=1)**2 + np.std(cp_vals, ddof=1)**2) / 2)
            cohens_d = (np.mean(wt_vals) - np.mean(cp_vals)) / pooled_sd if pooled_sd > 0 else 0
            t, p = stats.ttest_ind(wt_vals, cp_vals)
            all_pvalues.append(p)
            sig = "YES" if p < 0.05 else "no"
            print(f"  {'WT dimer vs Compound dimer':<35} {metric:<12} {np.mean(wt_vals):>8.4f} {np.mean(cp_vals):>8.4f} {cohens_d:>8.2f} {p:>10.4f} {sig:>5}")
            comparison_results.append({
                "comparison": "WT vs Compound",
                "metric": metric,
                "mean1": np.mean(wt_vals),
                "mean2": np.mean(cp_vals),
                "cohens_d": cohens_d,
                "p_raw": p,
                "significant": p < 0.05
            })

            # WT vs C677T
            if len(c677t_vals) >= 2:
                pooled_sd2 = np.sqrt((np.std(wt_vals, ddof=1)**2 + np.std(c677t_vals, ddof=1)**2) / 2)
                cohens_d2 = (np.mean(wt_vals) - np.mean(c677t_vals)) / pooled_sd2 if pooled_sd2 > 0 else 0
                t2, p2 = stats.ttest_ind(wt_vals, c677t_vals)
                all_pvalues.append(p2)
                sig2 = "YES" if p2 < 0.05 else "no"
                print(f"  {'WT dimer vs C677T dimer':<35} {metric:<12} {np.mean(wt_vals):>8.4f} {np.mean(c677t_vals):>8.4f} {cohens_d2:>8.2f} {p2:>10.4f} {sig2:>5}")

    # ============================================================
    # 3. Bonferroni correction
    # ============================================================
    print(f"\n[3/4] Applying Bonferroni correction ({len(all_pvalues)} comparisons)...")

    n_comparisons = len(all_pvalues)
    bonferroni_threshold = 0.05 / n_comparisons
    print(f"  Raw threshold: p < 0.05")
    print(f"  Bonferroni threshold: p < {bonferroni_threshold:.4f} (0.05 / {n_comparisons})")

    for i, p in enumerate(all_pvalues):
        adjusted_p = min(p * n_comparisons, 1.0)
        sig = "YES" if adjusted_p < 0.05 else "no"
        print(f"  Comparison {i+1}: raw p={p:.4f}, adjusted p={adjusted_p:.4f} ({sig})")

    # ============================================================
    # 4. Full per-residue pLDDT plot
    # ============================================================
    print("\n[4/4] Generating per-residue pLDDT comparison plot...")

    # Extract pLDDT from WT and compound dimer best models
    wt_cifs = sorted(RESULTS.glob("wt_mono_run1/*model_0.cif"))
    cp_cifs = sorted(RESULTS.glob("compound_dimer_run1/*model_0.cif"))
    c677t_cifs = sorted(RESULTS.glob("c677t_mono_run1/*model_0.cif"))

    fig, axes = plt.subplots(2, 1, figsize=(14, 8), sharex=True)

    for cif_list, label, color in [
        (wt_cifs, "WT", "#2196F3"),
        (c677t_cifs, "C677T (A222V)", "#FF9800"),
        (cp_cifs, "Compound (A222V+E429A)", "#F44336"),
    ]:
        if cif_list:
            plddts = extract_plddt_from_cif(cif_list[0])
            if plddts:
                residues = sorted(plddts.keys())
                values = [plddts[r] for r in residues]
                axes[0].plot(residues, values, label=label, alpha=0.8, linewidth=0.8, color=color)

    axes[0].axvline(x=222, color='red', linestyle='--', alpha=0.5, label='pos 222 (C677T)')
    axes[0].axvline(x=429, color='purple', linestyle='--', alpha=0.5, label='pos 429 (A1298C)')
    axes[0].set_ylabel('pLDDT', fontsize=12)
    axes[0].set_title('Per-Residue pLDDT Confidence: WT vs C677T vs Compound', fontsize=14)
    axes[0].legend(fontsize=9, loc='lower left')
    axes[0].set_ylim(40, 100)
    axes[0].axhspan(90, 100, alpha=0.1, color='blue', label='Very high confidence')
    axes[0].axhspan(70, 90, alpha=0.1, color='cyan')
    axes[0].axhspan(50, 70, alpha=0.1, color='yellow')

    # Difference plot (WT - Compound)
    if wt_cifs and cp_cifs:
        wt_plddts = extract_plddt_from_cif(wt_cifs[0])
        cp_plddts = extract_plddt_from_cif(cp_cifs[0])
        if wt_plddts and cp_plddts:
            common = sorted(set(wt_plddts.keys()) & set(cp_plddts.keys()))
            diff = [wt_plddts[r] - cp_plddts[r] for r in common]
            axes[1].bar(common, diff, width=1, color=['#F44336' if d < 0 else '#4CAF50' for d in diff], alpha=0.7)
            axes[1].axhline(y=0, color='black', linewidth=0.5)
            axes[1].axvline(x=222, color='red', linestyle='--', alpha=0.5)
            axes[1].axvline(x=429, color='purple', linestyle='--', alpha=0.5)
            axes[1].set_ylabel('pLDDT difference\n(WT - Compound)', fontsize=11)
            axes[1].set_xlabel('Residue Number', fontsize=12)
            axes[1].set_title('pLDDT Difference: Positive = WT higher confidence', fontsize=12)

    plt.tight_layout()
    plt.savefig(fig_dir / "plddt_full_chain.png", dpi=200, bbox_inches='tight')
    plt.close()
    print("  Saved plddt_full_chain.png")

    # ============================================================
    # Summary
    # ============================================================
    print("\n" + "=" * 60)
    print("VALIDATION SUMMARY")
    print("=" * 60)
    if rmsd is not None:
        print(f"  RMSD vs 6FCX:    {rmsd:.2f} A ({classification})")
    print(f"  Cohen's d (ipTM WT vs Compound): {comparison_results[0]['cohens_d']:.2f}" if comparison_results else "")
    print(f"  Raw p-value:     {comparison_results[0]['p_raw']:.4f}" if comparison_results else "")
    print(f"  Bonferroni adj:  {min(comparison_results[0]['p_raw'] * n_comparisons, 1.0):.4f}" if comparison_results else "")
    print(f"  Comparisons:     {n_comparisons}")
    print(f"  Figures saved:   {fig_dir}/plddt_full_chain.png")
    print("=" * 60)

if __name__ == "__main__":
    main()
