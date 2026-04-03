#!/usr/bin/env python3
"""Convert RESEARCH_PAPER_DRAFT.md to a publication-quality PDF with embedded figures and title page."""
import markdown, base64
from weasyprint import HTML
from pathlib import Path

SRC = Path("docs/RESEARCH_PAPER_DRAFT.md")
OUT = Path("docs/MTHFR_Research_Paper.pdf")
FIG = Path("analysis/outputs")

def img_b64(path, width="100%"):
    if not path.exists(): return ""
    data = base64.b64encode(path.read_bytes()).decode()
    return f'<img src="data:image/png;base64,{data}" style="width:{width};max-width:100%">'

def figure(path, num, caption, width="100%"):
    if not path.exists(): return ""
    return f'''<div class="figure">
        {img_b64(path, width)}
        <p class="caption"><b>Figure {num}.</b> {caption}</p>
    </div>'''

md_text = SRC.read_text()
html_body = markdown.markdown(md_text, extensions=["tables", "fenced_code"])

# Title page
title_page = f"""
<div class="title-page">
    <div class="title-badge">PREPRINT</div>
    <h1 class="main-title">MTHFR Variant Hypothesis Prioritization</h1>
    <p class="subtitle">Computational Prioritization of Selected MTHFR Variant States<br>for Experimental Follow-Up</p>
    <p class="subtitle2">A Hypothesis-Prioritization Study Using AlphaFold 3 and Boltz-2</p>
    <div class="title-divider"></div>
    <p class="author">Igor Mihaljko</p>
    <p class="affil">Independent Researcher | DSM.Promo | Chicago, IL</p>
    <p class="affil">ORCID: 0009-0000-1408-1065</p>
    <p class="affil">igor@dsm.promo</p>
    <div class="title-divider"></div>
    <p class="date">April 2026 | v1.6</p>
    <div class="key-result">
        <p class="kr-label">CORE COMPUTATIONAL OBSERVATION</p>
        <p class="kr-text">Monomer predictions showed preserved overall folding confidence across tested variant states.
        Dimer predictions showed comparative inter-chain differences, with the compound heterozygous
        dimer yielding the lowest reported interaction-confidence values.</p>
        <table class="kr-table">
            <tr><th>Metric</th><th>WT Dimer (n=10)</th><th>Compound Het (n=10)</th><th>Significance</th></tr>
            <tr><td>ipTM (AlphaFold 3)</td><td>0.752 &plusmn;0.023</td><td style="color:#CC3333;font-weight:bold">0.714 &plusmn;0.026</td><td>p=0.003</td></tr>
            <tr><td>FAD binding</td><td>0.563 &plusmn;0.014</td><td style="color:#CC3333;font-weight:bold">0.542 &plusmn;0.020</td><td>p=0.016</td></tr>
            <tr><td>pLDDT @ pos 429</td><td>96.04 &plusmn;0.22</td><td style="color:#CC3333;font-weight:bold">95.29 &plusmn;0.21</td><td>p&lt;0.000001</td></tr>
            <tr><td>Equilibrium RMSD (100ns MD)</td><td>8.17 &plusmn;0.30 &Aring;</td><td style="color:#CC3333;font-weight:bold">6.88 &plusmn;0.30 &Aring;</td><td>Cohen's d=4.31</td></tr>
        </table>
    </div>
    <p class="disclaimer-small">For research and educational purposes only. Computational predictions, not experimentally resolved structures.</p>
</div>
<div style="page-break-after:always"></div>
"""

# Figures section
figures_section = f"""
<div style="page-break-before:always"></div>
<h2>Figures</h2>

{figure(FIG/"figures"/"summary_dashboard.png", 1, "Summary dashboard of all 16 structural predictions across AlphaFold 3 (Jobs 1-12) and Boltz-2 (Jobs 13-16) platforms. Compound heterozygous dimer consistently yielded the lowest comparative confidence values.")}

{figure(FIG/"figures"/"clinical_targets.png", 2, "Primary experimental follow-up context. Retinal biomarkers are identified as the primary validation pathway based on quantifiable OCT endpoints, established animal models, and demonstrated genotype-response correlation in human studies.")}

<div style="page-break-before:always"></div>

{figure(FIG/"charts"/"iptm_comparison.png", 3, "Interface predicted TM-score (ipTM) comparison across variant states. Compound heterozygous dimers averaged lower ipTM than wild-type and single-variant comparators, consistent with a possible combined dimer-level perturbation in this modeling setup.")}

{figure(FIG/"charts"/"ptm_comparison.png", 4, "Predicted TM-score (pTM) comparison across all tested configurations. Monomer predictions remained broadly similar; dimer predictions showed comparative differences.")}

{figure(FIG/"figures"/"plddt_comparison.png", 5, "Per-residue confidence (pLDDT) at mutation sites 222 and 429 across variant states. Position 429 showed the largest comparative decrease in compound dimers (95.29 vs 96.04 in wild-type, n=10), consistent with possible regulatory-domain involvement.")}

{figure(FIG/"figures"/"md_comparison.png", 6, "100ns molecular dynamics comparison: WT vs compound heterozygous MTHFR dimer (OpenMM, Amber14/TIP3P-FB, 300K, RTX 4090). PBC-corrected per-chain RMSD shows compound dimer adopts a more compact ensemble (equilibrium RMSD 6.88 vs 8.17 A, Cohen's d = 4.31, p < 1e-323). Both systems reach equilibrium at ~62ns.")}

<div style="page-break-before:always"></div>
<h3>3D Structure Comparisons</h3>

<div class="grid">
{figure(FIG/"figures"/"structure_job01_wt_mono_fad.png", "7a", "Wild-type monomer + FAD (Job 01). Backbone colored by pLDDT confidence (blue = high).", "48%")}
{figure(FIG/"figures"/"structure_job06_compound_dimer_fad.png", "7b", "Compound heterozygous dimer + FAD (Job 06). Author genotype context.", "48%")}
</div>

<div class="grid">
{figure(FIG/"figures"/"structure_job02_wt_dimer_fad.png", "7c", "Wild-type dimer + FAD (Job 02). Baseline homodimer.", "48%")}
{figure(FIG/"figures"/"structure_job13_wt_dimer_fad_thf.png", "7d", "Wild-type dimer + FAD + THF (Job 13, Boltz-2). Substrate binding context.", "48%")}
</div>

<div style="page-break-before:always"></div>
<h3>Predicted Aligned Error (PAE) Heatmaps</h3>
<p class="caption">Blue regions indicate high positional confidence; yellow/red indicate lower confidence or disorder. Off-diagonal blocks in dimer plots represent inter-chain contact confidence.</p>

<div class="grid">
{figure(FIG/"pae_plots"/"pae_job01_wt_mono_fad.png", "8a", "WT Mono (Job 01)", "48%")}
{figure(FIG/"pae_plots"/"pae_job03_c677t_mono_fad.png", "8b", "C677T Mono (Job 03)", "48%")}
</div>
<div class="grid">
{figure(FIG/"pae_plots"/"pae_job02_wt_dimer_fad.png", "8c", "WT Dimer (Job 02)", "48%")}
{figure(FIG/"pae_plots"/"pae_job06_compound_dimer_fad.png", "8d", "Compound Dimer (Job 06)", "48%")}
</div>
<div class="grid">
{figure(FIG/"pae_plots"/"pae_job13_wt_dimer_fad_thf.png", "8e", "WT + THF, Boltz-2 (Job 13)", "48%")}
{figure(FIG/"pae_plots"/"pae_job15_compound_dimer_fad_thf.png", "8f", "Compound + THF, Boltz-2 (Job 15)", "48%")}
</div>
"""

html_doc = f"""<!DOCTYPE html>
<html><head><meta charset="utf-8">
<style>
@page {{
    size: letter;
    margin: 0.85in 0.85in 0.85in 0.85in;
    @bottom-center {{
        content: "Page " counter(page) " of " counter(pages);
        font-size: 8px; color: #999;
    }}
}}
@page:first {{
    margin: 0;
    @bottom-center {{ content: none; }}
}}
body {{ font-family: 'Helvetica Neue', Helvetica, Arial, sans-serif; font-size: 10.5pt; line-height: 1.5; color: #1a1a1a; }}
/* Title page */
.title-page {{
    width: 100%; min-height: 9in;
    display: flex; flex-direction: column; align-items: center; justify-content: center;
    text-align: center; padding: 0.8in 1.2in;
    background: linear-gradient(180deg, #f8fafc 0%, #e8eef5 50%, #f0f4f8 100%);
    box-sizing: border-box;
}}
.title-badge {{
    background: #2E75B6; color: white; padding: 4px 20px; border-radius: 15px;
    font-size: 11pt; font-weight: bold; letter-spacing: 2px; margin-bottom: 30px;
}}
.main-title {{ font-size: 22pt; color: #1B3A5C; margin: 10px 0; border: none; line-height: 1.2; }}
.subtitle {{ font-size: 13pt; color: #2E75B6; margin: 5px 0; }}
.subtitle2 {{ font-size: 11pt; color: #666; margin: 5px 0 20px 0; font-style: italic; }}
.title-divider {{ width: 200px; height: 2px; background: #2E75B6; margin: 15px auto; }}
.author {{ font-size: 14pt; color: #1B3A5C; font-weight: bold; margin: 5px 0; }}
.affil {{ font-size: 10pt; color: #666; margin: 2px 0; }}
.date {{ font-size: 11pt; color: #1B3A5C; font-weight: bold; margin: 15px 0; }}
.key-result {{
    background: white; border: 2px solid #2E75B6; border-radius: 8px;
    padding: 15px 20px; margin: 20px 0; max-width: 500px; box-shadow: 0 2px 8px rgba(0,0,0,0.1);
}}
.kr-label {{ font-size: 9pt; color: #2E75B6; font-weight: bold; letter-spacing: 1px; margin-bottom: 8px; }}
.kr-text {{ font-size: 9.5pt; color: #333; line-height: 1.4; margin: 0; }}
.kr-table {{ width: 100%; margin-top: 10px; font-size: 9pt; }}
.kr-table th {{ background: #1B3A5C; color: white; padding: 4px 8px; font-size: 8.5pt; }}
.kr-table td {{ padding: 3px 8px; border-bottom: 1px solid #eee; font-size: 8.5pt; }}
.disclaimer-small {{ font-size: 8pt; color: #999; margin-top: 20px; }}
/* Content */
h1 {{ font-size: 20pt; color: #1B3A5C; border-bottom: 3px solid #2E75B6; padding-bottom: 6px; margin-top: 0; }}
h2 {{ font-size: 14pt; color: #1B3A5C; border-bottom: 1px solid #ddd; padding-bottom: 4px; margin-top: 22px; page-break-after: avoid; }}
h3 {{ font-size: 11.5pt; color: #2E75B6; margin-top: 15px; page-break-after: avoid; }}
h4 {{ font-size: 10.5pt; color: #555; margin-top: 10px; }}
p {{ margin: 5px 0; text-align: justify; }}
table {{ width: 100%; border-collapse: collapse; margin: 10px 0; font-size: 9pt; page-break-inside: avoid; }}
th {{ background: #1B3A5C; color: white; padding: 5px 7px; text-align: left; font-weight: bold; }}
td {{ padding: 4px 7px; border-bottom: 1px solid #ddd; }}
tr:nth-child(even) {{ background: #f8f9fa; }}
blockquote {{ border-left: 4px solid #2E75B6; margin: 10px 0; padding: 8px 12px; background: #f0f4f8; font-style: italic; font-size: 9.5pt; }}
strong {{ color: #1B3A5C; }}
ul, ol {{ margin: 5px 0; padding-left: 20px; }}
li {{ margin: 2px 0; }}
hr {{ border: none; border-top: 2px solid #2E75B6; margin: 20px 0; }}
code {{ background: #f4f4f4; padding: 1px 3px; border-radius: 2px; font-size: 9pt; }}
a {{ color: #2E75B6; text-decoration: none; }}
/* Figures */
.figure {{ margin: 12px 0; text-align: center; page-break-inside: avoid; }}
.figure img {{ border: 1px solid #ddd; border-radius: 4px; box-shadow: 0 1px 3px rgba(0,0,0,0.08); }}
.caption {{ font-size: 8.5pt; color: #666; margin-top: 4px; text-align: center; font-style: italic; }}
.grid {{ display: flex; gap: 8px; justify-content: center; margin: 8px 0; flex-wrap: wrap; }}
.grid .figure {{ flex: 1; min-width: 45%; }}
.footer-note {{ margin-top: 25px; padding-top: 10px; border-top: 2px solid #ddd; font-size: 8pt; color: #888; text-align: center; }}
</style>
</head><body>
{title_page}
{html_body}
{figures_section}
<div class="footer-note">
    <p><b>MTHFR Variant Hypothesis Prioritization Program</b> | Igor Mihaljko | DSM.Promo | ORCID: 0009-0000-1408-1065</p>
    <p>GitHub: github.com/DSMPromo/mthfr-target-validation | CC BY-NC-SA 4.0</p>
    <p style="margin-top:6px">This research was developed by DSM.Promo as a demonstration of how AI-powered tools
    (AlphaFold 3, Boltz-2, Claude) can streamline complex computational research workflows.</p>
</div>
</body></html>"""

HTML(string=html_doc).write_pdf(str(OUT))
size_kb = OUT.stat().st_size / 1024
print(f"PDF saved to {OUT}")
print(f"File size: {size_kb:.0f} KB ({size_kb/1024:.1f} MB)")
