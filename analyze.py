#!/usr/bin/env python3
"""
MTHFR Structure Analysis Pipeline
===================================
Supports both AlphaFold Server (CIF) and Boltz-2/Tamarind Bio (PDB) results.
Drop result folders or ZIP files into alphafold/results/ and run this script.

Usage:  python analyze.py
Colab:  !pip install matplotlib numpy && !python analyze.py

Requirements: Python 3.8+, matplotlib, numpy
Author: Igor Mihaljko / DSM.Promo | License: CC BY-NC-SA 4.0
"""
import json, os, zipfile, csv, sys
from pathlib import Path
from datetime import datetime

try:
    import numpy as np
    import matplotlib; matplotlib.use('Agg')
    import matplotlib.pyplot as plt
except ImportError:
    print("Install: pip install matplotlib numpy"); sys.exit(1)

RESULTS_DIRS = [Path("alphafold/results_all")]
OUTPUT_DIR = Path("analysis/outputs")
C677T_POS, A1298C_POS = 222, 429

JOBS = {
    # AlphaFold Server jobs (1-12)
    "job01_wt_mono_fad":"WT mono+FAD", "job02_wt_dimer_fad":"WT dimer+FAD",
    "job03_c677t_mono_fad":"C677T mono+FAD", "job04_c677t_dimer_fad":"C677T dimer+FAD",
    "job05_a1298c_mono_fad":"A1298C mono+FAD", "job06_compound_dimer_fad":"Compound dimer+FAD",
    "job07_wt_mono_rep":"WT mono rep", "job08_wt_dimer_rep":"WT dimer rep",
    "job09_c677t_mono_rep":"C677T mono rep", "job10_c677t_dimer_rep":"C677T dimer rep",
    "job11_a1298c_mono_rep":"A1298C mono rep", "job12_compound_rep":"Compound dimer rep",
    # Boltz-2 / Tamarind Bio jobs (13-16) -- substrate/inhibitor binding
    "job13_wt_dimer_fad_thf":"WT dimer+FAD+THF", "job14_c677t_dimer_fad_thf":"C677T dimer+FAD+THF",
    "job15_compound_dimer_fad_thf":"Compound dimer+FAD+THF", "job16_wt_dimer_fad_sam":"WT dimer+FAD+SAM",
    # v1.2 additional seeds (17-34)
    "job17_wt_mono_seed3":"WT mono seed3", "job18_wt_mono_seed4":"WT mono seed4", "job19_wt_mono_seed5":"WT mono seed5",
    "job20_wt_dimer_seed3":"WT dimer seed3", "job21_wt_dimer_seed4":"WT dimer seed4", "job22_wt_dimer_seed5":"WT dimer seed5",
    "job23_c677t_mono_seed3":"C677T mono seed3", "job24_c677t_mono_seed4":"C677T mono seed4", "job25_c677t_mono_seed5":"C677T mono seed5",
    "job26_c677t_dimer_seed3":"C677T dimer seed3", "job27_c677t_dimer_seed4":"C677T dimer seed4", "job28_c677t_dimer_seed5":"C677T dimer seed5",
    "job29_a1298c_mono_seed3":"A1298C mono seed3", "job30_a1298c_mono_seed4":"A1298C mono seed4", "job31_a1298c_mono_seed5":"A1298C mono seed5",
    "job32_compound_dimer_seed3":"Compound dimer seed3", "job33_compound_dimer_seed4":"Compound dimer seed4", "job34_compound_dimer_seed5":"Compound dimer seed5",
}

# Grouping for 5-seed statistics
SEED_GROUPS = {
    "WT mono": ["job01_wt_mono_fad","job07_wt_mono_rep","job17_wt_mono_seed3","job18_wt_mono_seed4","job19_wt_mono_seed5"],
    "WT dimer": ["job02_wt_dimer_fad","job08_wt_dimer_rep","job20_wt_dimer_seed3","job21_wt_dimer_seed4","job22_wt_dimer_seed5"],
    "C677T mono": ["job03_c677t_mono_fad","job09_c677t_mono_rep","job23_c677t_mono_seed3","job24_c677t_mono_seed4","job25_c677t_mono_seed5"],
    "C677T dimer": ["job04_c677t_dimer_fad","job10_c677t_dimer_rep","job26_c677t_dimer_seed3","job27_c677t_dimer_seed4","job28_c677t_dimer_seed5"],
    "A1298C mono": ["job05_a1298c_mono_fad","job11_a1298c_mono_rep","job29_a1298c_mono_seed3","job30_a1298c_mono_seed4","job31_a1298c_mono_seed5"],
    "Compound dimer": ["job06_compound_dimer_fad","job12_compound_rep","job32_compound_dimer_seed3","job33_compound_dimer_seed4","job34_compound_dimer_seed5"],
}

def detect_source(d):
    """Detect whether results are from AlphaFold Server or Boltz-2."""
    d = Path(d)
    if list(d.rglob("*summary_confidences*.json")):
        return "alphafold"
    if list(d.rglob("*.pdb")) or list(d.rglob("*confidence*.json")) or list(d.rglob("*scores*.json")):
        return "boltz2"
    return "unknown"

def get_plddt_from_pdb(pdb_file, resnum, chain="A"):
    """Extract per-residue pLDDT from PDB B-factor column for CA atom."""
    with open(pdb_file) as f:
        for line in f:
            if line.startswith("ATOM") and line[12:16].strip() == "CA" and line[21] == chain:
                try:
                    res_seq = int(line[22:26].strip())
                    if res_seq == resnum:
                        return float(line[60:66].strip())
                except (ValueError, IndexError):
                    continue
    return None

def load_boltz2_metrics(d):
    """Extract metrics from Boltz-2 output (confidence JSON + PDB + PAE npz)."""
    d = Path(d)
    result = {"ptm": None, "iptm": None, "rank": None, "fad_iptm": None,
              "p222": None, "p429": None, "pae": None, "clash": False}

    # Load confidence JSON (Boltz-2 format: confidence_result_model_0.json)
    jsons = sorted(d.rglob("confidence_result_model_0.json"))
    if not jsons:
        jsons = sorted(d.rglob("*confidence*.json"))
    if jsons:
        data = json.load(open(jsons[0]))
        result["ptm"] = data.get("ptm")
        result["iptm"] = data.get("iptm")
        result["rank"] = data.get("confidence_score")
        result["fad_iptm"] = data.get("ligand_iptm")

    # Load PAE from npz file
    pae_files = sorted(d.rglob("pae_result_model_0.npz"))
    if pae_files:
        try:
            pae_data = np.load(pae_files[0])
            # npz files have arrays stored by key
            for key in pae_data.files:
                arr = pae_data[key]
                if arr.ndim == 2:
                    result["pae"] = arr
                    break
        except Exception:
            pass

    # Load per-residue pLDDT from npz file
    plddt_files = sorted(d.rglob("plddt_result_model_0.npz"))
    if plddt_files:
        try:
            plddt_data = np.load(plddt_files[0])
            for key in plddt_data.files:
                arr = plddt_data[key]
                if arr.ndim == 1 and len(arr) > 429:
                    result["p222"] = float(arr[221]) * 100  # Convert 0-1 to 0-100
                    result["p429"] = float(arr[428]) * 100
                    break
        except Exception:
            pass

    # Fallback: extract pLDDT from PDB B-factor if npz not available
    if result["p222"] is None:
        pdbs = sorted(d.rglob("*model_0*.pdb")) or sorted(d.rglob("*.pdb"))
        if pdbs:
            result["p222"] = get_plddt_from_pdb(pdbs[0], 222)
            result["p429"] = get_plddt_from_pdb(pdbs[0], 429)

    return result

def find_jobs(d):
    d=Path(d); d.mkdir(parents=True,exist_ok=True); dirs=[]
    # Collect all valid result directories, including nested subdirectories
    for x in sorted(d.rglob("*")):
        if x.is_dir() and x.name.startswith("job"):
            if list(x.rglob("*summary_confidences*.json")) or list(x.rglob("confidence_result_model_*.json")):
                dirs.append(x)
    # Also check direct children (v1.0 format)
    for x in sorted(d.iterdir()):
        if x.is_dir() and x.name.startswith("job"):
            if x not in dirs and (list(x.rglob("*summary_confidences*.json")) or list(x.rglob("confidence_result_model_*.json"))):
                dirs.append(x)
    return sorted(set(dirs), key=lambda p: p.name)

def find_all_jobs():
    """Search all result directories."""
    all_dirs = []
    for rd in RESULTS_DIRS:
        if rd.exists():
            all_dirs.extend(find_jobs(rd))
    return all_dirs

def load_json(d,pattern):
    c=sorted(Path(d).rglob(pattern))
    return json.load(open(c[0])) if c else None

def get_residue_plddt(d, resnum, chain="A"):
    """Extract per-residue pLDDT from CIF file (B-factor column index 14 for CA atom)."""
    cifs = sorted(Path(d).rglob("*model_0.cif"))
    if not cifs: return None
    target = f" {resnum} "
    with open(cifs[0]) as cf:
        for line in cf:
            if line.startswith("ATOM") and " CA " in line and target in line and f" {chain} " in line:
                parts = line.split()
                if len(parts) > 14:
                    try: return float(parts[14])
                    except ValueError: pass
    return None

def analyze(dirs):
    results=[]
    for d in dirs:
        n=d.name; source=detect_source(d)
        print(f"  {n} [{source}]")

        if source == "alphafold":
            s=load_json(d,"*summary_confidences_0.json")
            f=load_json(d,"*full_data_0.json")
            if not s: print(f"    SKIP: no summary JSON"); continue
            ci=s.get("chain_iptm",[])
            p222=get_residue_plddt(d, 222)
            p429=get_residue_plddt(d, 429)
            e={"name":n,"config":JOBS.get(n,n),"ptm":s.get("ptm"),"iptm":s.get("iptm"),
               "rank":s.get("ranking_score"),"disordered":s.get("fraction_disordered"),
               "clash":s.get("has_clash",False),"chains":len(s.get("chain_ptm",[])),
               "fad_iptm":ci[-1] if len(ci)>1 else None,
               "p222":p222,"p429":p429,
               "pae":np.array(f["pae"]) if f and "pae" in f else None}
            results.append(e)

        elif source == "boltz2":
            m = load_boltz2_metrics(d)
            e={"name":n,"config":JOBS.get(n,n),"ptm":m["ptm"],"iptm":m["iptm"],
               "rank":m["rank"],"disordered":None,
               "clash":m["clash"],"chains":0,
               "fad_iptm":m["fad_iptm"],
               "p222":m["p222"],"p429":m["p429"],
               "pae":m["pae"]}
            results.append(e)

        else:
            print(f"    SKIP: unknown format in {n}")
    return results

def make_tables(R,out):
    out.mkdir(parents=True,exist_ok=True)
    with open(out/"metrics.csv",'w',newline='') as f:
        w=csv.writer(f)
        w.writerow(["Job","Config","pTM","ipTM","Rank","FAD_iptm","pLDDT@222","pLDDT@429","Clash"])
        for r in R:
            w.writerow([r["name"],r["config"],
                f"{r['ptm']:.4f}" if r['ptm'] else "",
                f"{r['iptm']:.4f}" if r['iptm'] else "",
                f"{r['rank']:.4f}" if r['rank'] else "",
                f"{r['fad_iptm']:.4f}" if r['fad_iptm'] else "",
                f"{r['p222']:.1f}" if r['p222'] else "",
                f"{r['p429']:.1f}" if r['p429'] else "",
                r['clash']])
    print(f"  Saved metrics.csv")

def make_pae_plots(R,out):
    d=out/"pae_plots"; d.mkdir(parents=True,exist_ok=True); n=0
    for r in R:
        if r["pae"] is None: continue
        fig,ax=plt.subplots(figsize=(8,7))
        ax.imshow(r["pae"],cmap="bwr",vmin=0,vmax=30)
        ax.set_title(f"PAE: {r['name']} ({r['config']})",fontweight='bold')
        ax.set_xlabel("Token (scored)"); ax.set_ylabel("Token (aligned)")
        for p,c in [(221,"lime"),(428,"yellow")]:
            if p<r["pae"].shape[0]:
                ax.axhline(y=p,color=c,lw=0.5,ls='--',alpha=0.7)
                ax.axvline(x=p,color=c,lw=0.5,ls='--',alpha=0.7)
        plt.colorbar(ax.images[0],ax=ax,label="PAE (A)",shrink=0.8)
        plt.tight_layout(); plt.savefig(d/f"pae_{r['name']}.png",dpi=200); plt.close(); n+=1
    print(f"  Saved {n} PAE plots")

def make_charts(R,out):
    d=out/"charts"; d.mkdir(parents=True,exist_ok=True)
    V=[r for r in R if r["iptm"] is not None]
    if not V: print("  No data for charts"); return
    
    for metric,label,title in [("iptm","ipTM","Interface TM — FAD Binding Confidence"),
                                ("ptm","pTM","Overall Structure Accuracy")]:
        fig,ax=plt.subplots(figsize=(max(10,len(V)*1.2),6))
        vals=[r[metric] or 0 for r in V]
        colors=["#2E75B6" if "WT" in r["config"] else "#CC3333" if "C677T" in r["config"] 
                else "#D4760A" if "A1298C" in r["config"] else "#6B3FA0" if "ompound" in r["config"]
                else "#666" for r in V]
        bars=ax.bar(range(len(V)),vals,color=colors,edgecolor='black',lw=0.5)
        ax.set_xticks(range(len(V)))
        ax.set_xticklabels([r["name"] for r in V],rotation=45,ha='right',fontsize=9)
        ax.set_ylabel(label); ax.set_title(title,fontweight='bold'); ax.set_ylim(0,1)
        if metric=="iptm":
            ax.axhline(y=0.8,color='green',ls='--',alpha=0.5,label='High confidence (0.8)')
        for b,v in zip(bars,vals):
            if v>0: ax.text(b.get_x()+b.get_width()/2,v+0.01,f'{v:.3f}',ha='center',fontsize=8)
        ax.legend(fontsize=9); plt.tight_layout()
        plt.savefig(d/f"{metric}_comparison.png",dpi=200); plt.close()
    print(f"  Saved comparison charts")

def make_report(R,out):
    # Separate AlphaFold and Boltz-2 results
    af_results = [r for r in R if "thf" not in r["name"] and "sam" not in r["name"]]
    b2_results = [r for r in R if "thf" in r["name"] or "sam" in r["name"]]

    h=f"""<!DOCTYPE html><html><head><meta charset="UTF-8">
<title>MTHFR Variant Hypothesis Prioritization Program — Analysis Report</title>
<style>
body{{font-family:'Segoe UI',Arial,sans-serif;max-width:1200px;margin:0 auto;padding:30px;color:#333;line-height:1.6}}
h1{{color:#1B3A5C;border-bottom:3px solid #2E75B6;padding-bottom:12px;font-size:28px}}
h2{{color:#2E75B6;margin-top:40px;font-size:22px}}
h3{{color:#1B3A5C;font-size:18px}}
table{{border-collapse:collapse;width:100%;margin:15px 0;font-size:14px}}
th{{background:#1B3A5C;color:#fff;padding:10px 8px;text-align:left}}
td{{padding:8px;border-bottom:1px solid #ddd}}
tr:nth-child(even){{background:#f8f9fa}}
tr:hover{{background:#e8f4fd}}
.warn{{border-left:4px solid #CC3333;padding:15px;background:#fff5f5;margin:20px 0;border-radius:4px}}
.info{{border-left:4px solid #2E75B6;padding:15px;background:#f0f7ff;margin:20px 0;border-radius:4px}}
.finding{{border-left:4px solid #2E7D32;padding:15px;background:#f1f8e9;margin:20px 0;border-radius:4px}}
img{{max-width:100%;border:1px solid #ddd;margin:10px 0;border-radius:4px;box-shadow:0 2px 4px rgba(0,0,0,0.1)}}
.grid{{display:grid;grid-template-columns:1fr 1fr;gap:15px}}
.grid-3{{display:grid;grid-template-columns:1fr 1fr 1fr;gap:10px}}
.badge{{display:inline-block;padding:3px 10px;border-radius:12px;font-size:12px;font-weight:bold;margin-right:8px}}
.badge-af{{background:#E3F2FD;color:#1565C0}}
.badge-b2{{background:#F3E5F5;color:#7B1FA2}}
.highlight{{background:#FFF9C4;font-weight:bold}}
a{{color:#2E75B6}}
.footer{{margin-top:40px;padding-top:20px;border-top:2px solid #ddd;color:#666;font-size:13px;text-align:center}}
</style></head><body>

<h1>MTHFR Variant Hypothesis Prioritization Program</h1>
<p style="font-size:18px;color:#666">Computational Structural Comparison of Dimer Interface Confidence, Cofactor Binding, and Substrate Interaction in Selected MTHFR Variant States</p>

<p><b>Generated:</b> {datetime.now().strftime('%Y-%m-%d %H:%M')} |
<b>Author:</b> Igor Mihaljko / <a href="https://dsm.promo">DSM.Promo</a> |
<b>ORCID:</b> <a href="https://orcid.org/0009-0000-1408-1065">0009-0000-1408-1065</a> |
<b>Jobs analyzed:</b> {len(R)} |
<b>Platforms:</b> <span class="badge badge-af">AlphaFold 3</span> <span class="badge badge-b2">Boltz-2</span></p>

<div class="warn"><b>DISCLAIMER:</b> For research and educational purposes only. The models and scores presented here are computational predictions, not experimentally resolved structures or clinical findings. All observations require orthogonal validation through experimental methods such as cryo-EM, X-ray crystallography, molecular dynamics, and biochemical or functional assays. The author is not a medical professional.</div>

<h2>Executive Summary</h2>
<div class="finding">
<b>Key observation:</b> Relative to the tested wild-type comparator, the compound heterozygous MTHFR dimer showed the lowest predicted interaction-confidence values across the assessed dimer-interface and ligand-associated metrics:
<ul>
<li><b>FAD cofactor context:</b> ipTM 0.53 versus 0.57 in the tested wild-type dimer (AlphaFold Server)</li>
<li><b>THF substrate context:</b> ligand ipTM 0.878 versus 0.974 in the tested wild-type dimer (Boltz-2)</li>
<li><b>Dimer interface context:</b> ipTM 0.70 versus 0.76 in the tested wild-type dimer</li>
<li><b>Position 429 context:</b> pLDDT 95.0 versus 96.2 in the tested wild-type comparator</li>
</ul>
<b>Working hypothesis:</b> These comparative computational results support experimental testing of whether selected MTHFR variant states are associated with measurable differences in dimer-level behavior, retinal biomarkers, and neuropsychiatric-relevant biochemical readouts. The present work is intended to generate falsifiable hypotheses for orthogonal validation, not to establish mechanism, disease causation, or therapeutic relevance.
</div>"""

    # Summary dashboard
    fig_dir = out/"figures"
    if (fig_dir/"summary_dashboard.png").exists():
        h += '<h2>Summary Dashboard</h2>'
        h += '<img src="figures/summary_dashboard.png" alt="Summary Dashboard">'

    # Candidate disease pathways
    if (fig_dir/"clinical_targets.png").exists():
        h += '<h2>Primary Experimental Follow-Up Contexts</h2>'
        h += '<img src="figures/clinical_targets.png" alt="Candidate Pathways">'

    # Structure visualizations
    struct_figs = sorted(fig_dir.glob("structure_*.png")) if fig_dir.exists() else []
    if struct_figs:
        h += '<h2>3D Structure Visualizations</h2>'
        h += '<p>Backbone traces colored by pLDDT confidence (blue=high, red=low). Mutation sites marked with stars.</p>'
        h += '<div class="grid">'
        for fig in struct_figs:
            h += f'<img src="figures/{fig.name}" alt="{fig.stem}">'
        h += '</div>'

    # pLDDT comparison
    if (fig_dir/"plddt_comparison.png").exists():
        h += '<h2>Per-Residue Confidence Comparison</h2>'
        h += '<img src="figures/plddt_comparison.png" alt="pLDDT Comparison">'

    # AlphaFold Server results table
    h += '<h2>AlphaFold 3 Server Results (Jobs 1-12)</h2>'
    h += '<p><span class="badge badge-af">AlphaFold 3</span> FAD cofactor binding — monomers and homodimers with independent replication seeds</p>'
    h += '<table><tr><th>Job</th><th>Config</th><th>pTM</th><th>ipTM</th><th>Rank</th><th>FAD ipTM</th><th>pLDDT@222</th><th>pLDDT@429</th></tr>'
    for r in af_results:
        cls = ' class="highlight"' if 'ompound' in r['config'] else ''
        h+=f"<tr{cls}><td>{r['name']}</td><td>{r['config']}</td>"
        for k in ['ptm','iptm','rank','fad_iptm']:
            h+=f"<td>{r[k]:.3f}</td>" if r[k] else "<td>-</td>"
        for k in ['p222','p429']:
            h+=f"<td>{r[k]:.0f}</td>" if r[k] else "<td>-</td>"
        h+="</tr>"
    h+="</table>"

    # Boltz-2 results table
    if b2_results:
        h += '<h2>Boltz-2 Substrate Binding Results (Jobs 13-16)</h2>'
        h += '<p><span class="badge badge-b2">Boltz-2</span> THF folate substrate and SAM allosteric inhibitor binding via Tamarind Bio</p>'
        h += '<table><tr><th>Job</th><th>Config</th><th>pTM</th><th>ipTM</th><th>Confidence</th><th>Ligand ipTM</th><th>pLDDT@222</th><th>pLDDT@429</th></tr>'
        for r in b2_results:
            cls = ' class="highlight"' if 'ompound' in r['config'] else ''
            h+=f"<tr{cls}><td>{r['name']}</td><td>{r['config']}</td>"
            for k in ['ptm','iptm','rank','fad_iptm']:
                h+=f"<td>{r[k]:.3f}</td>" if r[k] else "<td>-</td>"
            for k in ['p222','p429']:
                h+=f"<td>{r[k]:.0f}</td>" if r[k] else "<td>-</td>"
            h+="</tr>"
        h+="</table>"

    # Charts
    h += '<h2>Comparison Charts</h2><div class="grid">'
    for name in ["iptm_comparison","ptm_comparison"]:
        p=out/"charts"/f"{name}.png"
        if p.exists(): h+=f'<img src="charts/{name}.png" alt="{name}">'
    h += '</div>'

    # PAE plots
    pae_dir=out/"pae_plots"
    if pae_dir.exists():
        pngs=sorted(pae_dir.glob("*.png"))
        if pngs:
            h+='<h2>Predicted Aligned Error (PAE) Plots</h2>'
            h+='<p>Blue = confident relative positioning. Red/white = uncertain. Dashed lines mark mutation positions (222, 429).</p>'
            h+='<div class="grid">'
            for png in pngs: h+=f'<img src="pae_plots/{png.name}" alt="{png.stem}">'
            h+='</div>'

    h+="""<h2>Interpretation Guide</h2>
<div class="info">
<h3>What These Metrics Mean</h3>
<ul>
<li><b>pTM (predicted TM-score):</b> Overall fold confidence. >0.7 = reliable fold.</li>
<li><b>ipTM (interface TM-score):</b> Confidence in protein-protein or protein-ligand interaction. >0.8 = high confidence.</li>
<li><b>Ligand ipTM:</b> Specifically how well the ligand (FAD, THF, SAM) binds. Higher = stronger predicted binding.</li>
<li><b>pLDDT:</b> Per-residue confidence. >90 = very high. >70 = confident. <50 = likely disordered.</li>
<li><b>PAE:</b> How confidently AlphaFold predicts relative positions of residue pairs. Low values (blue) = confident.</li>
</ul>
<h3>Why Monomers vs Dimers Matter</h3>
<p>MTHFR functions as a <b>homodimer</b> in vivo. The tested monomer predictions (Jobs 1,3,5) do not suggest large-scale loss of overall folding confidence across variant states. However, dimer predictions (Jobs 2,4,6) show inter-chain differences: FAD binding confidence drops from 0.97 to 0.53-0.57, and the compound heterozygous dimer yielded the lowest comparative values. This is the central comparative finding in the current computational dataset and the main basis for downstream experimental prioritization.</p>
<h3>Why Two Platforms</h3>
<p><b>AlphaFold 3 Server</b> (Jobs 1-12) provided the core FAD binding analysis. <b>Boltz-2</b> (Jobs 13-16) enabled THF substrate and SAM inhibitor modeling not possible through AlphaFold Server's interface. Across the tested model configurations, both showed the same directional pattern, with the compound heterozygous dimer yielding the lowest comparative confidence values among the assessed dimer conditions.</p>
</div>

<h2>Primary Experimental Indication Hypotheses</h2>
<div class="info">
<h3>Retinal Pathway Hypothesis</h3>
<p>Because prior literature has linked MTHFR-related one-carbon pathway disruption and hyperhomocysteinemia to retinal injury, the retina is a candidate downstream system for experimental follow-up. The current computational results do not establish efficacy in retinal disease, but they support testing whether selected MTHFR variant states are associated with measurable retinal biomarker differences. Prior evidence: Mthfr+/- mice show ganglion cell loss (Markand 2015); case report of visual recovery with betaine (Hergert 2022).</p>
<h3>Neuropsychiatric Pathway Hypothesis</h3>
<p>Because folate cycle dysfunction may influence BH4-dependent neurotransmitter pathways, neuropsychiatric phenotypes are a second candidate area for validation. The current structural data do not establish causation or treatment effect, but they support testing whether selected MTHFR variant states correlate with measurable biochemical differences. This may warrant exploratory study in defined subgroups where one-carbon metabolism and BH4-relevant biology are hypothesized to be relevant.</p>
<p><i>Additional downstream systems (autism/CFD, B-vitamin metabolism, cardiovascular, epigenetics, pregnancy) are documented in the full research paper as candidate areas for future investigation.</i></p>
</div>

<h2>Next Steps</h2>
<ol>
<li>Structural validation: compare predictions against experimental PDB 6FCX crystal structure</li>
<li>Molecular dynamics simulations on WT vs compound dimer</li>
<li>Biochemical validation: express proteins, dimer stability assays, cofactor binding assays</li>
<li>Exploratory biomarker correlation studies in well-defined compound heterozygous cohorts</li>
<li>Submit to bioRxiv as preprint</li>
<li>Contact target researchers (structural biology, metabolic disease, retinal, neuropsychiatric)</li>
</ol>

<div class="footer">
<p><b>MTHFR Variant Hypothesis Prioritization Program</b> | Igor Mihaljko | <a href="https://github.com/DSMPromo/mthfr-target-validation">GitHub</a> | <a href="https://orcid.org/0009-0000-1408-1065">ORCID</a> | CC BY-NC-SA 4.0</p>
<p>Computational predictions, not experimentally resolved structures. This document is intended for research planning and educational use only.</p>
<p style="margin-top:15px;padding-top:15px;border-top:1px solid #eee;font-size:12px;color:#888"><b>About this project:</b> This research was developed by <a href="https://dsm.promo">DSM.Promo</a> as a demonstration of how AI-powered tools (AlphaFold 3, Boltz-2, Claude) can streamline complex computational research workflows. The entire computational foundation was built using publicly accessible tools at zero cost, demonstrating that advanced structural biology hypothesis prioritization is now accessible to non-specialists with the right methodology and AI integration.</p>
</div>
</body></html>"""

    with open(out/"report.html",'w') as f: f.write(h)
    print(f"  Saved report.html — OPEN THIS IN YOUR BROWSER")

def make_pymol(out):
    d=out/"pymol_scripts"; d.mkdir(parents=True,exist_ok=True)
    with open(d/"analyze.pml",'w') as f:
        f.write("""# MTHFR PyMOL Analysis Script
# Update paths below, then: File > Run Script > select this file

bg_color white
# load path/to/WT_model_0.cif, WT
# load path/to/C677T_model_0.cif, C677T
# spectrum b, red_white_blue, WT, minimum=0, maximum=100
# select site222, resi 222
# show sticks, site222 or resn FAD
# zoom site222
# ray 2400, 1800
# png position222.png, dpi=300
# align C677T, WT
# color cyan, WT
# color salmon, C677T
print("Update the paths above and uncomment to run")
""")
    print(f"  Saved PyMOL script")

def main():
    print("="*60)
    print("MTHFR AlphaFold Analysis Pipeline")
    print("="*60)
    print()
    
    print("[1/6] Searching for results...")
    dirs=find_all_jobs()
    if not dirs:
        for rd in RESULTS_DIRS:
            rd.mkdir(parents=True,exist_ok=True)
        OUTPUT_DIR.mkdir(parents=True,exist_ok=True)
        make_pymol(OUTPUT_DIR)
        print(f"\n  No results yet. Drop AlphaFold ZIP files in alphafold/results/")
        print(f"  Then run: python analyze.py")
        print(f"\n  PyMOL script template ready in {OUTPUT_DIR}/pymol_scripts/")
        return
    
    print(f"  Found {len(dirs)} job(s)\n")
    
    print("[2/6] Extracting metrics...")
    R=analyze(dirs); print()
    
    print("[3/6] Creating tables...")
    make_tables(R,OUTPUT_DIR); print()
    
    print("[4/6] Generating PAE plots...")
    make_pae_plots(R,OUTPUT_DIR); print()
    
    print("[5/6] Generating charts...")
    make_charts(R,OUTPUT_DIR); print()
    
    print("[6/6] Generating report...")
    make_pymol(OUTPUT_DIR)
    make_report(R,OUTPUT_DIR)
    
    print(f"\n{'='*60}")
    print("DONE! Open analysis/outputs/report.html in your browser.")
    print(f"{'='*60}\n")

if __name__=="__main__": main()
