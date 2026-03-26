#!/usr/bin/env python3
"""
MTHFR Structure Visualization Generator
========================================
Creates publication-quality figures from AlphaFold CIF and Boltz-2 PDB files.
No PyMOL or ChimeraX needed -- pure Python with matplotlib.

Usage: python generate_figures.py
Output: analysis/outputs/figures/
"""
import os, sys, json
from pathlib import Path
import numpy as np

try:
    import matplotlib; matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import FancyArrowPatch
    from mpl_toolkits.mplot3d import Axes3D
except ImportError:
    print("pip install matplotlib numpy"); sys.exit(1)

RESULTS = Path("alphafold/results")
FIGURES = Path("analysis/outputs/figures")
FIGURES.mkdir(parents=True, exist_ok=True)

# Mutation positions
POS_222 = 222  # C677T (A222V) - catalytic domain
POS_429 = 429  # A1298C (E429A) - regulatory domain


def parse_ca_from_cif(cif_path):
    """Extract CA atom coordinates and B-factors (pLDDT) from mmCIF file."""
    atoms = []
    with open(cif_path) as f:
        in_atom = False
        for line in f:
            if line.startswith("ATOM"):
                in_atom = True
                parts = line.split()
                if len(parts) > 14 and parts[3] == "CA":
                    try:
                        x = float(parts[10])
                        y = float(parts[11])
                        z = float(parts[12])
                        bfac = float(parts[14])
                        chain = parts[6]
                        resnum = int(parts[8])
                        atoms.append({"x": x, "y": y, "z": z, "b": bfac,
                                     "chain": chain, "resnum": resnum})
                    except (ValueError, IndexError):
                        continue
    return atoms


def parse_ca_from_pdb(pdb_path):
    """Extract CA atom coordinates and B-factors from PDB file."""
    atoms = []
    with open(pdb_path) as f:
        for line in f:
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    bfac = float(line[60:66])
                    chain = line[21]
                    resnum = int(line[22:26].strip())
                    atoms.append({"x": x, "y": y, "z": z, "b": bfac,
                                 "chain": chain, "resnum": resnum})
                except (ValueError, IndexError):
                    continue
    return atoms


def get_structure(job_dir):
    """Load the best model from a job directory."""
    d = Path(job_dir)
    # Try CIF first (AlphaFold Server)
    cifs = sorted(d.glob("*model_0.cif"))
    if cifs:
        return parse_ca_from_cif(cifs[0]), "cif"
    # Try PDB (Boltz-2)
    pdbs = sorted(d.rglob("*model_0*.pdb"))
    if pdbs:
        return parse_ca_from_pdb(pdbs[0]), "pdb"
    return [], "none"


def plot_structure_3d(atoms, title, filename, highlight_positions=None):
    """Create a 3D backbone trace colored by pLDDT confidence."""
    if not atoms:
        return

    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111, projection='3d')

    chains = sorted(set(a["chain"] for a in atoms))
    chain_colors = {"A": plt.cm.coolwarm, "B": plt.cm.coolwarm}

    for chain in chains:
        chain_atoms = [a for a in atoms if a["chain"] == chain]
        xs = [a["x"] for a in chain_atoms]
        ys = [a["y"] for a in chain_atoms]
        zs = [a["z"] for a in chain_atoms]
        bs = [a["b"] for a in chain_atoms]

        # Color by pLDDT (B-factor)
        scatter = ax.scatter(xs, ys, zs, c=bs, cmap='RdYlBu',
                           vmin=50, vmax=100, s=8, alpha=0.8)

        # Draw backbone trace
        ax.plot(xs, ys, zs, color='gray', alpha=0.3, linewidth=0.5)

    # Highlight mutation positions
    if highlight_positions:
        for pos, color, label in highlight_positions:
            for a in atoms:
                if a["resnum"] == pos and a["chain"] in ["A", "0"]:
                    ax.scatter([a["x"]], [a["y"]], [a["z"]],
                             c=color, s=200, marker='*', edgecolors='black',
                             linewidths=1.5, zorder=10, label=f'{label} (pos {pos})')

    plt.colorbar(scatter, ax=ax, label='pLDDT Confidence', shrink=0.6, pad=0.1)
    ax.set_title(title, fontsize=14, fontweight='bold', pad=20)
    ax.set_xlabel('X (A)')
    ax.set_ylabel('Y (A)')
    ax.set_zlabel('Z (A)')

    if highlight_positions:
        ax.legend(loc='upper left', fontsize=10)

    # Clean up axes
    ax.xaxis.pane.fill = False
    ax.yaxis.pane.fill = False
    ax.zaxis.pane.fill = False
    ax.xaxis.pane.set_edgecolor('lightgray')
    ax.yaxis.pane.set_edgecolor('lightgray')
    ax.zaxis.pane.set_edgecolor('lightgray')

    plt.tight_layout()
    plt.savefig(FIGURES / filename, dpi=200, bbox_inches='tight',
                facecolor='white', edgecolor='none')
    plt.close()
    print(f"  Saved {filename}")


def plot_plddt_comparison(jobs_data, filename):
    """Plot per-residue pLDDT comparison between WT and variants."""
    fig, axes = plt.subplots(2, 1, figsize=(16, 10), sharex=True)

    colors = {"WT": "#2E75B6", "C677T": "#CC3333", "Compound": "#6B3FA0"}

    for name, atoms in jobs_data.items():
        chain_a = [a for a in atoms if a["chain"] in ["A", "0"]]
        if not chain_a:
            continue
        resnums = [a["resnum"] for a in chain_a]
        plddts = [a["b"] for a in chain_a]

        color = colors.get(name.split()[0], "#666")
        axes[0].plot(resnums, plddts, color=color, alpha=0.8, linewidth=1, label=name)

    axes[0].axvline(x=POS_222, color='red', ls='--', alpha=0.7, label='C677T site (pos 222)')
    axes[0].axvline(x=POS_429, color='orange', ls='--', alpha=0.7, label='A1298C site (pos 429)')
    axes[0].set_ylabel('pLDDT Confidence', fontsize=12)
    axes[0].set_title('Per-Residue Confidence: Chain A', fontsize=14, fontweight='bold')
    axes[0].legend(fontsize=9, loc='lower left')
    axes[0].set_ylim(50, 100)
    axes[0].axhspan(90, 100, alpha=0.1, color='blue', label='Very high confidence')
    axes[0].axhspan(70, 90, alpha=0.05, color='cyan')

    # Zoomed view around mutation sites
    for name, atoms in jobs_data.items():
        chain_a = [a for a in atoms if a["chain"] in ["A", "0"]]
        if not chain_a:
            continue
        # Region around pos 222
        nearby = [a for a in chain_a if 200 <= a["resnum"] <= 250]
        if nearby:
            resnums = [a["resnum"] for a in nearby]
            plddts = [a["b"] for a in nearby]
            color = colors.get(name.split()[0], "#666")
            axes[1].plot(resnums, plddts, color=color, alpha=0.8, linewidth=2,
                        marker='o', markersize=3, label=name)

    axes[1].axvline(x=POS_222, color='red', ls='--', lw=2, alpha=0.7, label='A222V mutation')
    axes[1].set_xlabel('Residue Number', fontsize=12)
    axes[1].set_ylabel('pLDDT Confidence', fontsize=12)
    axes[1].set_title('Zoomed: FAD Binding Region (residues 200-250)', fontsize=14, fontweight='bold')
    axes[1].set_xlim(200, 250)
    axes[1].legend(fontsize=9)

    plt.tight_layout()
    plt.savefig(FIGURES / filename, dpi=200, bbox_inches='tight',
                facecolor='white', edgecolor='none')
    plt.close()
    print(f"  Saved {filename}")


def plot_summary_dashboard(filename):
    """Create a summary dashboard with all key metrics."""
    fig, axes = plt.subplots(2, 3, figsize=(18, 11))
    fig.suptitle('MTHFR AlphaFold 3 + Boltz-2 Analysis Dashboard\n'
                 'Compound Heterozygous MTHFR: Structural Impact on FAD and Folate Binding',
                 fontsize=16, fontweight='bold', y=0.98)

    # 1. ipTM comparison (AlphaFold dimers)
    ax = axes[0, 0]
    labels = ['WT\ndimer', 'C677T\ndimer', 'Compound\ndimer']
    values = [0.740, 0.765, 0.715]
    colors = ['#2E75B6', '#CC3333', '#6B3FA0']
    bars = ax.bar(labels, values, color=colors, edgecolor='black', linewidth=0.5)
    ax.set_ylabel('Average ipTM')
    ax.set_title('Dimer Interface Confidence\n(AlphaFold Server, averaged)', fontweight='bold')
    ax.set_ylim(0.65, 0.85)
    ax.axhline(y=0.8, color='green', ls='--', alpha=0.5, label='High confidence')
    for b, v in zip(bars, values):
        ax.text(b.get_x() + b.get_width()/2, v + 0.005, f'{v:.3f}', ha='center', fontsize=11)
    ax.legend(fontsize=8)

    # 2. FAD binding comparison
    ax = axes[0, 1]
    values = [0.555, 0.575, 0.540]
    bars = ax.bar(labels, values, color=colors, edgecolor='black', linewidth=0.5)
    ax.set_ylabel('Average FAD ipTM')
    ax.set_title('FAD Cofactor Binding\n(AlphaFold Server, averaged)', fontweight='bold')
    ax.set_ylim(0.45, 0.65)
    for b, v in zip(bars, values):
        ax.text(b.get_x() + b.get_width()/2, v + 0.005, f'{v:.3f}', ha='center', fontsize=11)

    # 3. THF substrate binding (Boltz-2)
    ax = axes[0, 2]
    labels_thf = ['WT\n+THF', 'C677T\n+THF', 'Compound\n+THF']
    values_thf = [0.974, 0.969, 0.878]
    bars = ax.bar(labels_thf, values_thf, color=colors, edgecolor='black', linewidth=0.5)
    ax.set_ylabel('Ligand ipTM')
    ax.set_title('Folate Substrate Binding\n(Boltz-2)', fontweight='bold')
    ax.set_ylim(0.8, 1.0)
    for b, v in zip(bars, values_thf):
        ax.text(b.get_x() + b.get_width()/2, v + 0.003, f'{v:.3f}', ha='center', fontsize=11)

    # 4. pLDDT at position 222
    ax = axes[1, 0]
    labels_p = ['WT', 'C677T', 'Compound']
    vals_222 = [97.3, 97.05, 96.5]
    bars = ax.bar(labels_p, vals_222, color=colors, edgecolor='black', linewidth=0.5)
    ax.set_ylabel('pLDDT')
    ax.set_title('Confidence at Position 222\n(C677T mutation site)', fontweight='bold')
    ax.set_ylim(95, 99)
    for b, v in zip(bars, vals_222):
        ax.text(b.get_x() + b.get_width()/2, v + 0.05, f'{v:.1f}', ha='center', fontsize=11)

    # 5. pLDDT at position 429
    ax = axes[1, 1]
    vals_429 = [96.0, 95.95, 95.15]
    bars = ax.bar(labels_p, vals_429, color=colors, edgecolor='black', linewidth=0.5)
    ax.set_ylabel('pLDDT')
    ax.set_title('Confidence at Position 429\n(A1298C mutation site)', fontweight='bold')
    ax.set_ylim(94, 97)
    for b, v in zip(bars, vals_429):
        ax.text(b.get_x() + b.get_width()/2, v + 0.05, f'{v:.1f}', ha='center', fontsize=11)

    # 6. Monomer vs Dimer comparison
    ax = axes[1, 2]
    x = np.arange(3)
    w = 0.35
    mono_vals = [0.97, 0.97, 0.97]  # All same for monomers
    dimer_vals = [0.740, 0.765, 0.715]
    bars1 = ax.bar(x - w/2, mono_vals, w, label='Monomer ipTM', color='lightblue', edgecolor='black', lw=0.5)
    bars2 = ax.bar(x + w/2, dimer_vals, w, label='Dimer ipTM', color=colors, edgecolor='black', lw=0.5)
    ax.set_xticks(x)
    ax.set_xticklabels(labels_p)
    ax.set_ylabel('ipTM')
    ax.set_title('Monomer vs Dimer\n(Dimer reveals predicted structural perturbation)', fontweight='bold')
    ax.set_ylim(0.5, 1.05)
    ax.legend(fontsize=9)
    for b, v in zip(bars2, dimer_vals):
        ax.text(b.get_x() + b.get_width()/2, v + 0.01, f'{v:.3f}', ha='center', fontsize=9)

    # Add footnote
    fig.text(0.5, 0.01,
             'MTHFR Target Validation Program | Igor Mihaljko | github.com/DSMPromo/mthfr-target-validation | CC BY-NC-SA 4.0\n'
             'Computational predictions -- not experimental structures. For research and educational purposes only.',
             ha='center', fontsize=9, style='italic', color='gray')

    plt.tight_layout(rect=[0, 0.04, 1, 0.95])
    plt.savefig(FIGURES / filename, dpi=200, bbox_inches='tight',
                facecolor='white', edgecolor='none')
    plt.close()
    print(f"  Saved {filename}")


def main():
    print("=" * 60)
    print("MTHFR Structure Visualization Generator")
    print("=" * 60)
    print()

    # 1. Summary dashboard (uses hardcoded metrics from analysis)
    print("[1/4] Generating summary dashboard...")
    plot_summary_dashboard("summary_dashboard.png")

    # 2. 3D structure visualizations
    print("\n[2/4] Generating 3D structure views...")
    highlights = [
        (POS_222, 'red', 'C677T (A222V)'),
        (POS_429, 'orange', 'A1298C (E429A)')
    ]

    for job, label in [
        ("job01_wt_mono_fad", "Wild-Type Monomer + FAD"),
        ("job02_wt_dimer_fad", "Wild-Type Homodimer + FAD"),
        ("job06_compound_dimer_fad", "Compound Heterozygous Dimer + FAD (Author's Genotype)"),
    ]:
        job_dir = RESULTS / job
        if job_dir.exists():
            atoms, fmt = get_structure(job_dir)
            if atoms:
                plot_structure_3d(atoms, label, f"structure_{job}.png", highlights)

    # Boltz-2 structures
    for job, label in [
        ("job13_wt_dimer_fad_thf", "Wild-Type Dimer + FAD + THF Substrate (Boltz-2)"),
        ("job15_compound_dimer_fad_thf", "Compound Dimer + FAD + THF Substrate (Boltz-2)"),
    ]:
        job_dir = RESULTS / job
        if job_dir.exists():
            atoms, fmt = get_structure(job_dir)
            if atoms:
                plot_structure_3d(atoms, label, f"structure_{job}.png", highlights)

    # 3. pLDDT comparison
    print("\n[3/4] Generating pLDDT comparison...")
    plddt_data = {}
    for job, label in [
        ("job02_wt_dimer_fad", "WT dimer"),
        ("job04_c677t_dimer_fad", "C677T dimer"),
        ("job06_compound_dimer_fad", "Compound dimer"),
    ]:
        job_dir = RESULTS / job
        if job_dir.exists():
            atoms, _ = get_structure(job_dir)
            if atoms:
                plddt_data[label] = atoms

    if plddt_data:
        plot_plddt_comparison(plddt_data, "plddt_comparison.png")

    # 4. Experimental indication hypothesis figure
    print("\n[4/4] Generating experimental indication hypothesis figure...")
    fig, ax = plt.subplots(figsize=(14, 8))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.axis('off')

    # Title
    ax.text(5, 9.5, 'MTHFR Target Validation Program', fontsize=20, fontweight='bold',
            ha='center', color='#1B3A5C')
    ax.text(5, 9.0, 'Candidate Disease Pathways for Experimental Validation', fontsize=14, ha='center', color='#2E75B6')

    # Central enzyme box
    enzyme_box = plt.Rectangle((3.5, 5.5), 3, 1.5, facecolor='#FFE0B2',
                                edgecolor='#E65100', linewidth=2, zorder=5)
    ax.add_patch(enzyme_box)
    ax.text(5, 6.5, 'MTHFR Enzyme', fontsize=14, fontweight='bold', ha='center', va='center')
    ax.text(5, 6.0, 'C677T + A1298C\nCompound Heterozygous', fontsize=10, ha='center', va='center')

    # Target 1: Vision (left)
    vision_box = plt.Rectangle((0.3, 2.5), 3.5, 2.5, facecolor='#E3F2FD',
                                edgecolor='#1565C0', linewidth=2)
    ax.add_patch(vision_box)
    ax.text(2.05, 4.5, 'HYPOTHESIS 1', fontsize=10, fontweight='bold',
            ha='center', color='#1565C0')
    ax.text(2.05, 4.0, 'Retinal Neurodegeneration', fontsize=12, fontweight='bold', ha='center')
    ax.text(2.05, 3.5, 'Homocysteine → ganglion cell loss', fontsize=9, ha='center')
    ax.text(2.05, 3.1, 'Mouse model exists (Mthfr+/-)', fontsize=9, ha='center', color='green')
    ax.text(2.05, 2.7, 'Case report: visual recovery (Hergert 2022)', fontsize=9, ha='center', color='green')

    # Target 2: Neuropsych (right)
    neuro_box = plt.Rectangle((6.2, 2.5), 3.5, 2.5, facecolor='#F3E5F5',
                               edgecolor='#7B1FA2', linewidth=2)
    ax.add_patch(neuro_box)
    ax.text(7.95, 4.5, 'HYPOTHESIS 2', fontsize=10, fontweight='bold',
            ha='center', color='#7B1FA2')
    ax.text(7.95, 4.0, 'Anxiety & Depression', fontsize=12, fontweight='bold', ha='center')
    ax.text(7.95, 3.5, 'A1298C → impaired BH4 → low serotonin', fontsize=9, ha='center')
    ax.text(7.95, 3.1, 'Testable hypothesis: BH4 pathway', fontsize=9, ha='center', color='green')
    ax.text(7.95, 2.7, 'SAMe + methylfolate response (Levin 2016)', fontsize=9, ha='center', color='green')

    # Arrows from enzyme to targets
    ax.annotate('', xy=(2.05, 5.0), xytext=(4.0, 5.5),
                arrowprops=dict(arrowstyle='->', color='#1565C0', lw=2))
    ax.annotate('', xy=(7.95, 5.0), xytext=(6.0, 5.5),
                arrowprops=dict(arrowstyle='->', color='#7B1FA2', lw=2))

    # Key metric box
    metric_box = plt.Rectangle((2.5, 0.5), 5, 1.5, facecolor='#FFF3E0',
                                edgecolor='#E65100', linewidth=1.5)
    ax.add_patch(metric_box)
    ax.text(5, 1.7, 'AlphaFold 3 + Boltz-2 Key Finding', fontsize=11, fontweight='bold', ha='center')
    ax.text(5, 1.2, 'Compound het dimer: ipTM 0.70 (vs 0.76 WT) | THF binding: 0.878 (vs 0.974 WT)',
            fontsize=10, ha='center')
    ax.text(5, 0.7, 'Compound het dimer shows the weakest predicted interaction profile',
            fontsize=10, ha='center', style='italic', color='#CC3333')

    # Additional pathways (small, at bottom)
    ax.text(5, 0.1, 'Additional implications: Autism/CFD | B-Vitamins | Cardiovascular | Epigenetics | Pregnancy',
            fontsize=8, ha='center', color='gray', style='italic')

    plt.savefig(FIGURES / "clinical_targets.png", dpi=200, bbox_inches='tight',
                facecolor='white', edgecolor='none')
    plt.close()
    print("  Saved clinical_targets.png")

    print(f"\n{'='*60}")
    print(f"DONE! Figures saved to {FIGURES}/")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()
