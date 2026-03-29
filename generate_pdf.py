#!/usr/bin/env python3
"""Convert RESEARCH_PAPER_DRAFT.md to a publication-quality PDF."""
import markdown
from weasyprint import HTML
from pathlib import Path

SRC = Path("docs/RESEARCH_PAPER_DRAFT.md")
OUT = Path("docs/MTHFR_Research_Paper.pdf")

md_text = SRC.read_text()
html_body = markdown.markdown(md_text, extensions=["tables", "fenced_code"])

html_doc = f"""<!DOCTYPE html>
<html><head><meta charset="utf-8">
<style>
@page {{
    size: letter;
    margin: 1in 1in 1in 1in;
    @bottom-center {{
        content: "Page " counter(page) " of " counter(pages);
        font-size: 9px;
        color: #999;
    }}
}}
body {{
    font-family: 'Helvetica Neue', Helvetica, Arial, sans-serif;
    font-size: 11pt;
    line-height: 1.6;
    color: #1a1a1a;
    max-width: 100%;
}}
h1 {{
    font-size: 22pt;
    color: #1B3A5C;
    border-bottom: 3px solid #2E75B6;
    padding-bottom: 8px;
    margin-top: 0;
}}
h2 {{
    font-size: 16pt;
    color: #1B3A5C;
    border-bottom: 1px solid #ddd;
    padding-bottom: 5px;
    margin-top: 30px;
    page-break-after: avoid;
}}
h3 {{
    font-size: 13pt;
    color: #2E75B6;
    margin-top: 20px;
    page-break-after: avoid;
}}
h4 {{
    font-size: 11pt;
    color: #555;
    margin-top: 15px;
}}
p {{
    margin: 8px 0;
    text-align: justify;
}}
table {{
    width: 100%;
    border-collapse: collapse;
    margin: 15px 0;
    font-size: 10pt;
    page-break-inside: avoid;
}}
th {{
    background: #1B3A5C;
    color: white;
    padding: 8px 10px;
    text-align: left;
    font-weight: bold;
}}
td {{
    padding: 6px 10px;
    border-bottom: 1px solid #ddd;
}}
tr:nth-child(even) {{
    background: #f8f9fa;
}}
blockquote {{
    border-left: 4px solid #2E75B6;
    margin: 15px 0;
    padding: 10px 15px;
    background: #f0f4f8;
    font-style: italic;
}}
strong {{
    color: #1B3A5C;
}}
ul, ol {{
    margin: 8px 0;
    padding-left: 25px;
}}
li {{
    margin: 4px 0;
}}
hr {{
    border: none;
    border-top: 2px solid #2E75B6;
    margin: 30px 0;
}}
code {{
    background: #f4f4f4;
    padding: 2px 5px;
    border-radius: 3px;
    font-size: 10pt;
}}
a {{
    color: #2E75B6;
    text-decoration: none;
}}
.title-block {{
    text-align: center;
    margin-bottom: 30px;
    padding: 20px;
    background: linear-gradient(135deg, #f0f4f8, #e8eef5);
    border-radius: 8px;
}}
.disclaimer {{
    background: #fff3cd;
    border: 1px solid #ffc107;
    padding: 12px;
    border-radius: 4px;
    font-size: 10pt;
    margin: 20px 0;
}}
.footer-note {{
    margin-top: 40px;
    padding-top: 15px;
    border-top: 2px solid #ddd;
    font-size: 9pt;
    color: #888;
    text-align: center;
}}
</style>
</head><body>
{html_body}
<div class="footer-note">
    <p>MTHFR Variant Hypothesis Prioritization Program | Igor Mihaljko | DSM.Promo | ORCID: 0009-0000-1408-1065</p>
    <p>GitHub: github.com/DSMPromo/mthfr-target-validation | CC BY-NC-SA 4.0</p>
    <p style="margin-top:10px">This research was developed by DSM.Promo as a demonstration of how AI-powered tools
    (AlphaFold 3, Boltz-2, Claude) can streamline complex computational research workflows.</p>
</div>
</body></html>"""

HTML(string=html_doc).write_pdf(str(OUT))
print(f"PDF saved to {OUT}")
print(f"File size: {OUT.stat().st_size / 1024:.0f} KB")
