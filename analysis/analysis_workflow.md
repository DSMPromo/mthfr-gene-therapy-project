# Analysis Workflow: Comparing AlphaFold Results

## Step 1: Download and Organize

```
MTHFR_Project/
├── Job01_WT_mono/
│   ├── fold_*_model_0.cif       ← Best ranked structure
│   ├── fold_*_model_1.cif       ← 2nd ranked
│   ├── fold_*_model_2.cif       ← etc.
│   ├── fold_*_summary_confidences_0.json
│   ├── fold_*_full_data_0.json
│   └── fold_*_job_request.json
├── Job02_WT_dimer/
├── Job03_C677T_mono/
├── Job04_C677T_dimer/
├── Job05_A1298C_mono/
├── Job06_compound_dimer/
└── ...
```

## Step 2: Extract Confidence Metrics

Open each `summary_confidences_0.json` file and record these values in `metrics_template.csv`:

| Metric | Where to Find It | What It Means |
|--------|-----------------|---------------|
| `ptm` | Top-level field | Overall structure accuracy (0-1, higher = better) |
| `iptm` | Top-level field | Interface accuracy — KEY METRIC for FAD binding |
| `ranking_score` | Top-level field | Combined quality score for ranking |
| `chain_iptm` | Array field | Per-chain interface scores (find the FAD chain) |
| `chain_ptm` | Array field | Per-chain structural accuracy |
| `has_clash` | Boolean field | Check for atomic clashes (should be false) |

## Step 3: PyMOL Visualization

### Install PyMOL
- Free educational version: https://pymol.org/edu/
- Or install via conda: `conda install -c conda-forge pymol-open-source`

### Load and Color a Structure
```pymol
# Load the best-ranked structure
load fold_Job01_model_0.cif, WT_mono

# Color by pLDDT confidence (stored in B-factor field)
spectrum b, red_white_blue, minimum=0, maximum=100

# Blue = high confidence (>90)
# White = moderate (70-90)
# Red = low confidence (<50)
```

### Zoom to Mutation Sites
```pymol
# Position 222 (C677T site)
select site222, resi 222
zoom site222
show sticks, site222
show sticks, resn FAD

# Take screenshot
ray 2400, 1800
png WT_position222.png, dpi=300

# Position 429 (A1298C site)
select site429, resi 429
zoom site429
show sticks, site429
png WT_position429.png, dpi=300
```

### Superimpose Wild-Type and Variant
```pymol
# Load both structures
load fold_Job01_model_0.cif, WT
load fold_Job03_model_0.cif, C677T

# Align variant onto wild-type
align C677T, WT

# Color differently
color cyan, WT
color salmon, C677T

# Show both at position 222
select both222, resi 222
show sticks, both222
show sticks, resn FAD
zoom both222

# Screenshot
ray 2400, 1800
png WT_vs_C677T_position222.png, dpi=300
```

### Measure FAD-Protein Distances
```pymol
# Use the measurement wizard
# Wizard → Measurement → click two atoms

# Or use the command:
distance dist1, /WT//A/222/CA, /WT//B/FAD/N5
# (adjust chain IDs and atom names based on your structure)
```

## Step 4: PAE Plot Comparison

The PAE (Predicted Aligned Error) plot is in `full_data_0.json` under the `pae` field.

### Using ChimeraX (recommended)
1. Open ChimeraX
2. Load the .cif file
3. Load the full_data JSON: `alphafold pae #1 file full_data_0.json`
4. The PAE plot appears automatically
5. Blue = low error (high confidence), Red = high error

### Using Python
```python
import json
import matplotlib.pyplot as plt
import numpy as np

with open('fold_Job01_full_data_0.json') as f:
    data = json.load(f)

pae = np.array(data['pae'])
plt.figure(figsize=(10, 8))
plt.imshow(pae, cmap='bwr_r', vmin=0, vmax=30)
plt.colorbar(label='Predicted Aligned Error (Å)')
plt.title('PAE Plot - Wild-Type MTHFR')
plt.xlabel('Token index')
plt.ylabel('Token index')
plt.savefig('pae_WT.png', dpi=300, bbox_inches='tight')
plt.show()
```

## Step 5: Create Comparison Table

Use the `metrics_template.csv` to create a publication-quality table:

| Job | Configuration | pTM | ipTM | chain_iptm (FAD) | pLDDT@222 | pLDDT@429 | Notes |
|-----|--------------|-----|------|-------------------|-----------|-----------|-------|
| 1 | WT mono + FAD | ? | ? | ? | ? | ? | Baseline |
| 2 | WT dimer + FAD | ? | ? | ? | ? | ? | Native state |
| 3 | C677T mono + FAD | ? | ? | ? | ? | ? | Catalytic domain |
| 4 | C677T dimer + FAD | ? | ? | ? | ? | ? | TT genotype |
| 5 | A1298C mono + FAD | ? | ? | ? | ? | ? | Regulatory domain |
| 6 | Compound dimer | ? | ? | ? | ? | ? | YOUR genotype |

## Step 6: Write Findings

Template sentence:
> "AlphaFold 3 structural modeling revealed that the C677T variant (A222V) showed a [X%] reduction in ipTM score at the protein-FAD interface compared to wild-type ([ipTM_WT] vs [ipTM_C677T]), consistent with the known destabilization of FAD cofactor binding in the thermolabile enzyme. The compound heterozygous dimer (modeling the author's genotype) showed [describe findings]."

## What to Save for the White Paper
- [ ] Metrics comparison table (filled CSV)
- [ ] Side-by-side PyMOL images: position 222 (WT vs C677T)
- [ ] Side-by-side PyMOL images: position 429 (WT vs A1298C)
- [ ] PAE plot comparisons (WT vs each variant)
- [ ] Superimposed structure overview image
- [ ] FAD distance measurements
- [ ] Written analysis paragraph
