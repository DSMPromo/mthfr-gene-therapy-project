#!/usr/bin/env python3
"""
MTHFR AlphaFold Analysis Pipeline
===================================
Drop AlphaFold Server ZIP files into alphafold/results/ and run this script.

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

JOBS = {"Job01":"WT mono+FAD","Job02":"WT dimer+FAD","Job03":"C677T mono+FAD",
        "Job04":"C677T dimer+FAD","Job05":"A1298C mono+FAD","Job06":"Compound dimer+FAD",
        "Job07":"WT mono rep","Job08":"WT dimer rep","Job09":"C677T mono rep",
        "Job10":"C677T dimer rep","Job11":"A1298C mono rep","Job12":"Compound rep",
        "Job13":"WT+FAD+THF","Job14":"C677T+FAD+THF","Job15":"Compound+FAD+THF","Job16":"WT+FAD+SAM"}

def find_jobs(d):
    d=Path(d); d.mkdir(parents=True,exist_ok=True); dirs=[]
    for z in sorted(d.glob("*.zip")):
        ed=d/z.stem
        if not ed.exists():
            print(f"  Extracting: {z.name}")
            with zipfile.ZipFile(z,'r') as zf: zf.extractall(ed)
        dirs.append(ed)
    for x in sorted(d.iterdir()):
        if x.is_dir() and x not in dirs and list(x.rglob("*summary_confidences*.json")): dirs.append(x)
    return dirs

def load_json(d,pattern):
    c=sorted(Path(d).rglob(pattern))
    return json.load(open(c[0])) if c else None

def analyze(dirs):
    results=[]
    for d in dirs:
        n=d.name; print(f"  {n}")
        s=load_json(d,"*summary_confidences_0.json")
        f=load_json(d,"*full_data_0.json")
        if not s: print(f"    SKIP: no summary JSON"); continue
        ci=s.get("chain_iptm",[])
        e={"name":n,"config":JOBS.get(n,n),"ptm":s.get("ptm"),"iptm":s.get("iptm"),
           "rank":s.get("ranking_score"),"disordered":s.get("fraction_disordered"),
           "clash":s.get("has_clash",False),"chains":len(s.get("chain_ptm",[])),
           "fad_iptm":ci[-1] if len(ci)>1 else None,
           "p222":f["atom_plddts"][221] if f and len(f.get("atom_plddts",[]))>221 else None,
           "p429":f["atom_plddts"][428] if f and len(f.get("atom_plddts",[]))>428 else None,
           "pae":np.array(f["pae"]) if f and "pae" in f else None}
        results.append(e)
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
<li>Update white paper with findings</li>
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
