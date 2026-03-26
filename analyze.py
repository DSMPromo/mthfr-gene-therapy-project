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

RESULTS_DIR = Path("alphafold/results")
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
    # Collect all valid result directories (skip ZIPs if folder already exists)
    for x in sorted(d.iterdir()):
        if x.is_dir():
            if list(x.rglob("*summary_confidences*.json")) or list(x.rglob("confidence_result_model_*.json")):
                dirs.append(x)
    return dirs

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
    h=f"""<!DOCTYPE html><html><head><meta charset="UTF-8">
<title>MTHFR Analysis Report</title>
<style>
body{{font-family:Arial;max-width:1100px;margin:0 auto;padding:20px;color:#333}}
h1{{color:#1B3A5C;border-bottom:3px solid #2E75B6;padding-bottom:10px}}
h2{{color:#2E75B6}}
table{{border-collapse:collapse;width:100%;margin:15px 0}}
th{{background:#1B3A5C;color:#fff;padding:8px}}
td{{padding:6px 8px;border-bottom:1px solid #ddd}}
tr:nth-child(even){{background:#f2f2f2}}
.warn{{border:2px solid #CC3333;padding:12px;background:#fff5f5;margin:15px 0}}
img{{max-width:100%;border:1px solid #ddd;margin:8px 0}}
.grid{{display:grid;grid-template-columns:1fr 1fr;gap:15px}}
</style></head><body>
<h1>MTHFR AlphaFold 3 Analysis Report</h1>
<p><b>Generated:</b> {datetime.now().strftime('%Y-%m-%d %H:%M')}</p>
<p><b>Author:</b> Igor Mihaljko / DSM.Promo | <b>Jobs analyzed:</b> {len(R)}</p>
<div class="warn"><b>DISCLAIMER:</b> For research/educational purposes only. Computational predictions require experimental validation. Author is not a medical professional.</div>
<h2>Metrics Comparison</h2>
<table><tr><th>Job</th><th>Config</th><th>pTM</th><th>ipTM</th><th>Rank</th><th>FAD iptm</th><th>pLDDT@222</th><th>pLDDT@429</th></tr>"""
    for r in R:
        h+=f"<tr><td>{r['name']}</td><td>{r['config']}</td>"
        for k in ['ptm','iptm','rank','fad_iptm']:
            h+=f"<td>{r[k]:.3f}</td>" if r[k] else "<td>-</td>"
        for k in ['p222','p429']:
            h+=f"<td>{r[k]:.0f}</td>" if r[k] else "<td>-</td>"
        h+="</tr>"
    h+="</table>"
    
    # Add chart images
    for name in ["iptm_comparison","ptm_comparison"]:
        p=out/"charts"/f"{name}.png"
        if p.exists(): h+=f'<h2>{name.replace("_"," ").title()}</h2><img src="charts/{name}.png">'
    
    # Add PAE plots
    pae_dir=out/"pae_plots"
    if pae_dir.exists():
        pngs=sorted(pae_dir.glob("*.png"))
        if pngs:
            h+='<h2>PAE Plots</h2><div class="grid">'
            for png in pngs: h+=f'<img src="pae_plots/{png.name}">'
            h+='</div>'
    
    h+="""<h2>Next Steps</h2><ol>
<li>Open CIF files in PyMOL — use scripts in pymol_scripts/</li>
<li>Superimpose WT and variant structures</li>
<li>Measure FAD distances at active site</li>
<li>Update research paper with findings</li>
<li>Share with research collaborators</li></ol></body></html>"""
    
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
    dirs=find_jobs(RESULTS_DIR)
    if not dirs:
        RESULTS_DIR.mkdir(parents=True,exist_ok=True)
        OUTPUT_DIR.mkdir(parents=True,exist_ok=True)
        make_pymol(OUTPUT_DIR)
        print(f"\n  No results yet. Drop AlphaFold ZIP files in {RESULTS_DIR}/")
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
