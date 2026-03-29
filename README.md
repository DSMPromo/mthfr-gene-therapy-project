# MTHFR Variant Target Validation Program

### Computational Prioritization of Selected MTHFR Variant States for Experimental Validation

> This project is a computational hypothesis-prioritization study designed to identify experimentally testable differences among wild-type, single-variant, and compound heterozygous MTHFR states. In the tested AlphaFold 3 and Boltz-2 model configurations, the compound heterozygous dimer yielded the lowest reported interaction-confidence values across the assessed dimer-interface and ligand-associated metrics relative to the tested wild-type and single-variant comparators. These observations do not establish mechanism, disease causation, clinical relevance, or therapeutic benefit. They support a bounded experimental agenda focused on dimer-level biology, retinal biomarkers, and neuropsychiatric-relevant biochemical readouts.
>
> **Core observation:** Across the tested model configurations, monomer predictions showed preserved overall folding confidence, whereas dimer predictions showed comparative inter-chain differences. This is the central computational finding of the current study.

[![License: CC BY-NC-SA 4.0](https://img.shields.io/badge/License-CC%20BY--NC--SA%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by-nc-sa/4.0/)
[![Status: Active Research](https://img.shields.io/badge/Status-Active%20Research-green.svg)]()
[![AlphaFold 3](https://img.shields.io/badge/AlphaFold-3%20Server-blue.svg)](https://alphafoldserver.com)
[![Sequences: UniProt P42898](https://img.shields.io/badge/UniProt-P42898-teal.svg)](https://www.uniprot.org/uniprotkb/P42898/entry)
[![Structure: PDB 6FCX](https://img.shields.io/badge/PDB-6FCX-purple.svg)](https://www.rcsb.org/structure/6FCX)

---

## Table of Contents

- [Who I Am](#who-i-am-and-who-i-am-not)
- [The Problem](#the-problem-nobody-connected)
- [The Science](#the-mthfr-enzyme)
- [Primary Targets](#primary-research-targets)
- [Additional Pathways](#additional-disease-pathways)
- [AlphaFold 3 Results](#alphafold-3-structural-predictions)
- [Research Findings](#what-the-research-shows)
- [Safety Architecture](#safety-architecture-what-if-something-goes-wrong)
- [Quick Start: Replicate This](#quick-start-replicate-this-research)
- [Development Roadmap](#development-roadmap)
- [Project Structure](#project-structure)
- [References](#verified-references)
- [How to Help](#how-you-can-help)

---

## About the Author

My name is Igor Mihaljko. I come from outside biomedical research. I'm a cybersecurity specialist and cloud architect with 20+ years in IT. I approached this as a structured computational investigation after discovering I carry compound heterozygous MTHFR variants (C677T + A1298C).

When I started researching what this means, I found literature spanning multiple disease categories and built this project to synthesize that literature alongside structural predictions in one open framework. I built this as an open starting point for researchers who want to evaluate these questions experimentally. The sequences, protocols, and analysis pipeline are included for replication and critique.

See [FOUNDERS_NOTE.md](FOUNDERS_NOTE.md) for the full story behind this project.

---

## Disclaimer

**This project is for research and educational purposes only.** Nothing in this document constitutes medical advice. All therapeutic concepts are speculative and would need years of validation by qualified researchers. The author is not a medical professional. See [DISCLAIMER.md](DISCLAIMER.md) for full details.

---

## A Cross-Specialty Problem That May Be Underintegrated in Current Research

MTHFR dysfunction has been studied extensively within individual medical specialties -- ophthalmology, psychiatry, neurology, cardiology, nutrition, reproductive medicine, and genetics. However, these investigations have largely proceeded in parallel rather than as an integrated body of evidence.

This project attempts to bring the structural biology perspective into that conversation. We do not claim to have connected dots that others missed; we are synthesizing published literature and adding computational structural data to support experimental prioritization.

### Possible Reasons MTHFR Has Not Yet Advanced as a Gene-Editing Research Program

As of March 2026, we did not identify a registered interventional clinical trial focused on MTHFR gene-editing in ClinicalTrials.gov searches. This absence should not be interpreted as evidence of scientific disinterest or target invalidity. Several factors may contribute:

| Observation | Context |
|-------------|---------|
| Common variant, not classic disease gene | High population prevalence has led these variants to be framed as common polymorphisms rather than direct interventional targets |
| Supplementation available | Methylfolate, riboflavin, B12 partially compensate, reducing perceived urgency for genetic correction |
| Cross-specialty fragmentation | Relevant disease associations span 7+ specialties; no single field owns the integrated picture |
| Structural data only recently available | First human MTHFR crystal structure 2018 (PDB 6FCX), cryo-EM 2024 (PDB 8QA5, 8QA6) |
| Gene-editing precedent emerging | Human base-editing programs have entered clinical testing or early clinical use for selected targets, though target biology, delivery, benefit-risk, and regulatory path remain highly target-specific |

**Why this question is now technically more testable:** Recent advances in structure prediction, human base-editing programs for other targets, scalable nucleic acid delivery systems, and new human MTHFR structural data make experimental assessment of MTHFR correction more technically discussable than before. They do not establish MTHFR as a validated therapeutic target. This project asks whether these conditions support a bounded experimental investigation.

---

## The MTHFR Enzyme

**MTHFR** (methylenetetrahydrofolate reductase, EC 1.5.1.20) is a FAD-dependent homodimeric enzyme that converts 5,10-methylenetetrahydrofolate to 5-methyltetrahydrofolate, the primary circulating form of folate and the methyl donor for converting homocysteine to methionine.

### Known Structures

| PDB | Year | Method | Resolution | What It Shows |
|-----|------|--------|------------|---------------|
| [6FCX](https://www.rcsb.org/structure/6FCX) | 2018 | X-ray | 2.5 A | First human MTHFR crystal structure -- homodimer with catalytic TIM-barrel + regulatory SAM-binding domain (Froese et al., *Nat Commun*) |
| [8QA5](https://www.rcsb.org/structure/8QA5) | 2024 | Cryo-EM | -- | MTHFR + SAH dis-inhibited (active) state -- single SAH bound, catalytic domain flexible (Froese et al., *Nat Commun*) |
| [8QA6](https://www.rcsb.org/structure/8QA6) | 2024 | Cryo-EM | -- | MTHFR + SAM inhibited state -- dual SAM binding reorients catalytic domain, blocks substrate access |

### The Two Key Variants

| Variant | rsID | Position | Domain | Mechanism | Population Frequency |
|---------|------|----------|--------|-----------|---------------------|
| **C677T (A222V)** | [rs1801133](https://www.ncbi.nlm.nih.gov/clinvar/?term=rs1801133) | 222 (catalytic) | FAD-binding TIM-barrel | Displaces helix alpha-5, disrupting FAD-interacting residues Asn168, Arg171, Lys172. Creates thermolabile enzyme. ~35% activity loss per allele. | 24-50% (varies by ethnicity) |
| **A1298C (E429A)** | [rs1801131](https://www.ncbi.nlm.nih.gov/clinvar/?term=rs1801131) | 429 (regulatory) | SAM-binding domain | Has been discussed in relation to BH4-related pathways and lies near the regulatory region implicated by 2024 cryo-EM structures. Reported activity reduction in some studies, with interpretation varying across study design and context. | ~25% European |

**Compound heterozygosity** has been reported in prior literature to associate with lower overall enzyme activity relative to wild-type, with possible contributions from both catalytic and regulatory context.

---

## Author Genotype Context

| Gene | Variant | rsID | Zygosity | Classification | Relevance |
|------|---------|------|----------|----------------|-----------|
| **MTHFR** | **C677T (A222V)** | rs1801133 | **Hetero (AG)** | Clinically Significant | Catalytic-domain variant with literature support for altered enzyme stability and reduced activity in some contexts |
| **MTHFR** | **A1298C (E429A)** | rs1801131 | **Hetero (TG)** | Conflicting/Uncertain | Regulatory-domain variant discussed in relation to one-carbon and neurotransmitter-relevant pathways, with mixed clinical interpretation across sources |
| **MTRR** | **I22M (c.66A>G)** | rs1801394 | **Hetero (AG)** | Likely Pathogenic | Variant in downstream B12-related metabolism that could be relevant to one-carbon pathway context, pending phenotype-specific interpretation. |

> **Voluntary Disclosure:** I am sharing my genetic variant data voluntarily because I believe transparency strengthens this research. This is my personal decision. You should never feel pressured to share your genetic information publicly. If you want to contribute your own data to this project, you can do so anonymously.

---

## Primary Experimental Indication Hypotheses

The computational findings in this project generate specific, testable hypotheses rather than therapeutic claims. Two downstream systems are prioritized for experimental validation based on existing literature:

### Hypothesis 1: Retinal Pathway

Because prior literature has linked MTHFR-related one-carbon pathway disruption and hyperhomocysteinemia to retinal injury, the retina is a candidate downstream system for validation. The current computational results do not establish efficacy in retinal disease, but they support testing whether selected MTHFR states are associated with measurable retinal biomarkers.

**Supporting literature:**
- **Mthfr+/- mice** show 2x retinal homocysteine, ~20% ganglion cell loss, thinner nerve fiber layer, and vascular leakage, without elevated IOP (Markand et al., 2015, *IOVS*)
- **Case report:** betaine treatment was associated with visual recovery in a patient with severe MTHFR deficiency (Hergert et al., 2022, *J Neurology*)
- **Retinal perfusion improved** by MTHFR genotype-guided supplementation; compound C677T/A1298C carriers showed significant improvement (Jiang et al., 2023, *Clin Ophthalmol*)
- MTHFR modifies CRB1-related retinopathies, worsening rd8 retinal phenotype (Markand et al., 2015, *Exp Eye Res*)

**Proposed validation:** Correlate MTHFR genotype with retinal OCT measurements, homocysteine levels, and visual acuity in defined patient cohorts.

### Hypothesis 2: Neuropsychiatric Pathway

Because one-carbon pathway disruption may influence BH4-dependent neurotransmitter biology, neuropsychiatric phenotypes are a second candidate area for validation. The current structural data do not establish causation or treatment effect, but they support testing whether selected MTHFR states are associated with measurable biochemical and clinical features.

**Supporting literature:**
- **MTHFR polymorphism** associated with schizophrenia, major depression, autism, bipolar disorder, and ADHD (Wan et al., 2018, *Transl Psychiatry*)
- Compound heterozygous patients show treatment-responsive anxiety with SAMe and methylated B vitamins when MTHFR status is identified (Levin & Varga, 2016)
- A1298C has been discussed in relation to BH4-relevant pathways in the literature, but the mechanism, effect size, and subgroup specificity remain incompletely defined
- MTHFR-related one-carbon pathway alterations may influence BH4-dependent neurotransmitter synthesis, warranting investigation in defined subgroups with biochemical abnormalities

**Proposed validation:** Measure BH4, homocysteine, SAMe/SAH ratio, and neurotransmitter metabolites in compound heterozygous individuals vs. wild-type controls.

---

## Additional Disease Pathways

MTHFR dysfunction has been associated with additional disease-relevant pathways in the literature. If future biochemical and in vivo studies support a causal contribution of selected MTHFR states, sequence correction strategies could be evaluated as one possible upstream intervention concept. These pathways are documented here for completeness.

Published literature has also reported MTHFR associations with autism/cerebral folate deficiency (Pu et al., 2020; Frye et al., 2018), B-vitamin processing failure (McNulty et al., 2006), cardiovascular disease via homocysteine elevation, transgenerational epigenetic effects in mouse models (Karahan et al., 2021), and adverse pregnancy outcomes including preeclampsia (Xia et al., 2013) and recurrent pregnancy loss (Kumar et al., 2024). These are documented in the [full research paper](docs/RESEARCH_PAPER_DRAFT.md) but are not the primary focus of this computational study.

**Working hypothesis:** The compound heterozygous MTHFR dimer shows lower predicted confidence for cofactor-associated and substrate-associated interactions relative to wild-type. These comparative computational results support experimental testing of whether selected MTHFR variant states are associated with measurable differences in retinal and neuropsychiatric-relevant biomarkers or pathway readouts, and whether sequence correction strategies could theoretically alter enzyme behavior relative to wild-type, subject to extensive experimental validation.

---

## AlphaFold 3 Structural Predictions

We used [AlphaFold 3 Server](https://alphafoldserver.com) to predict structures of wild-type, C677T, A1298C, and compound heterozygous MTHFR in complex with FAD cofactor.

### Methodology

- **Sequences:** All derived from canonical [UniProt P42898](https://www.uniprot.org/uniprotkb/P42898/entry) (656 amino acids) with verified mutations at positions 222 and 429
- **Predictions:** 16 jobs total (monomers, homodimers, heterodimers, with FAD, THF, and SAM)
- **Replication:** 10 independent seeds per configuration (64 total predictions) enabling statistical comparison (t-test)
- **Upload-ready JSON files** included in `alphafold/jobs/json/` for one-click replication

### Job Design

| Jobs | Protein | Configuration | Purpose |
|------|---------|--------------|---------|
| 01-02 | Wild-type | Monomer + Dimer + FAD | Baseline |
| 03-04 | C677T (A222V) | Monomer + Dimer + FAD | Catalytic domain variant |
| 05 | A1298C (E429A) | Monomer + FAD | Regulatory domain variant |
| 06 | Compound + WT | Heterodimer + FAD | Compound heterozygous (author's genotype) |
| 07-12 | All variants | Same as 01-06 | Replication with independent seeds |
| 13-15 | WT, C677T, Compound | Dimer + FAD + THF | Substrate binding comparison |
| 16 | Wild-type | Dimer + FAD + SAM | Allosteric inhibitor binding |

### AlphaFold 3 Results (March 2026)

#### Monomer Predictions

| Job | Variant | pTM | ipTM (FAD) | FAD Binding | pLDDT@222 | pLDDT@429 |
|-----|---------|-----|-----------|-------------|-----------|-----------|
| 01 | **WT mono** | 0.80 | 0.97 | 0.97 | 98.5 | 97.5 |
| 03 | **C677T mono** | 0.81 | 0.97 | 0.97 | 98.1 | 97.8 |
| 05 | **A1298C mono** | 0.81 | 0.97 | 0.97 | 98.5 | 97.5 |
| 07 | **WT rep** | 0.81 | 0.98 | 0.98 | 98.5 | 97.5 |
| 09 | **C677T rep** | 0.82 | 0.97 | 0.97 | 98.0 | 97.6 |
| 11 | **A1298C rep** | 0.81 | 0.97 | 0.97 | 98.4 | 97.5 |

#### Dimer Predictions (All 12 Jobs Complete, Replicated)

| Job | Variant | pTM | ipTM | FAD Binding | pLDDT@222 | pLDDT@429 |
|-----|---------|-----|------|-------------|-----------|-----------|
| 02 | **WT dimer** | 0.79 | 0.76 | 0.57 | 97.4 | 96.2 |
| 08 | **WT dimer rep** | 0.76 | 0.72 | 0.54 | 97.2 | 95.8 |
| 04 | **C677T dimer** | 0.77 | 0.75 | 0.57 | 97.0 | 96.0 |
| 10 | **C677T dimer rep** | 0.80 | 0.78 | 0.58 | 97.1 | 95.9 |
| 06 | **Compound dimer (author genotype context)** | **0.73** | **0.70** | **0.53** | **96.6** | **95.0** |
| 12 | **Compound dimer rep** | **0.76** | **0.73** | **0.55** | **96.4** | **95.3** |

#### Summary: Averaged Across Replicates

| Variant | Avg pTM | Avg ipTM | Avg FAD Binding | Avg pLDDT@222 | Avg pLDDT@429 |
|---------|---------|----------|-----------------|---------------|---------------|
| **WT dimer** (n=10) | 0.786 +/-0.016 | 0.752 +/-0.023 | 0.566 +/-0.018 | 97.3 | 96.0 |
| **C677T dimer** | 0.785 | 0.765 | 0.575 | 97.05 | 95.95 |
| **Compound dimer** (n=10) | **0.744 +/-0.027** | **0.714 +/-0.026** | **0.540 +/-0.023** | **96.5** | **95.3** |

**Key observations (replicated across independent seeds):**
- **Monomer predictions showed preserved overall folding confidence** (ipTM 0.97-0.98) -- the tested monomer predictions do not suggest large-scale loss of overall fold confidence
- **Dimer predictions reveal comparative inter-chain differences** -- reported FAD-associated confidence values are lower in the tested dimer models than in the tested monomer models, consistent with the distinct interaction context of the homodimer
- **The compound heterozygous dimer yielded the lowest reported interaction-confidence values** across the reported comparative metrics in both runs: pTM (0.73/0.76), ipTM (0.70/0.73), FAD binding (0.53/0.55), pLDDT@429 (95.0/95.3)
- **Compound heterozygous dimers averaged lower than the tested comparators** -- ipTM 0.714 +/-0.026 vs WT 0.752 +/-0.023 (t-test: t=3.437, **p=0.003**), pTM p=0.003, pLDDT@429 p<0.000001. Three of four metrics now survive Bonferroni correction for 12 comparisons, consistent with a robust combined dimer-level perturbation in this modeling setup
- **10-seed replication** -- independent random seeds produced consistent directional trends across the tested configurations, with 3 of 4 key metrics surviving Bonferroni correction (ipTM p=0.035, pTM p=0.031, pLDDT@429 p=0.000005)
- **Position 429 showed the largest confidence decrease** among the reported local metrics in compound dimers -- pLDDT 95.0-95.3 vs 95.8-96.2 in WT, consistent with possible regulatory-domain involvement at the dimer level

> **The core computational observation:** Monomer predictions showed preserved overall folding confidence, whereas dimer predictions showed comparative inter-chain differences. This is the central computational finding around which the experimental agenda is organized.

### Structural Validation and Statistical Analysis

**RMSD against experimental structure:** The wild-type monomer prediction was compared to the experimentally determined crystal structure (PDB: 6FCX, Froese et al., 2018) using gemmi-based superposition. The resulting RMSD of 1.42 A over 596 matched CA atoms classifies the prediction as high accuracy (RMSD <= 2.0 A), supporting the use of AlphaFold 3 predictions for comparative analysis within this protein system.

**Effect sizes (Cohen's d, WT dimer vs compound dimer, n=10 each):**

| Metric | WT Dimer (n=10) | Compound Dimer (n=10) | Cohen's d | Raw p | Bonferroni p |
|--------|---------------|---------------------|-----------|-------|-------------|
| ipTM | 0.752 +/-0.023 | 0.714 +/-0.026 | 1.54 (large) | 0.003 | **0.035** |
| pTM | 0.786 +/-0.016 | 0.744 +/-0.027 | 1.56 (large) | 0.003 | **0.031** |
| FAD binding | 0.566 +/-0.018 | 0.540 +/-0.023 | 1.19 (large) | 0.016 | 0.188 |
| pLDDT@429 | 96.04 +/-0.22 | 95.29 +/-0.21 | 3.45 (very large) | <0.000001 | **0.000005** |

After Bonferroni correction for 12 comparisons, three metrics remain statistically significant: ipTM (adjusted p=0.035), pTM (adjusted p=0.031), and pLDDT@429 (adjusted p=0.000005). This identifies both the dimer interface and the regulatory domain as robust comparative signals in this modeling setup.

**Per-residue pLDDT comparison:** A full-chain pLDDT comparison across all 656 residues was generated for WT, C677T, and compound heterozygous states (see analysis/outputs/figures/plddt_full_chain.png).

> **Important:** These are computational predictions, not experimental structures. All confidence metrics should be interpreted as hypothesis generators, not proof of mechanism. See the [full research paper draft](docs/RESEARCH_PAPER_DRAFT.md) for complete methodology and limitations.

---

## Structural and Literature Context

The 2024 cryo-EM structures (PDB 8QA5, 8QA6) revealed for the first time how SAM-mediated allosteric inhibition works through **dual SAM binding** that reorients the catalytic domain and blocks substrate access. This provides structural context relevant to interpretation of the A1298C (E429A) variant, which sits in the regulatory domain near the SAM-binding pocket.

The 2006 structural perturbation study (Pejchal et al.) showed that the C677T (A222V) variant displaces **helix alpha-5**, which carries three FAD-interacting residues (Asn168, Arg171, Lys172). This is consistent with the observed thermolability; the variant enzyme may have reduced FAD cofactor stability at physiological temperature.

### The MTRR Compounding Effect

The author also carries MTRR I22M (rs1801394), which is involved in B12 regeneration, the enzyme immediately downstream of MTHFR. This combination could place additive pressure on one-carbon metabolism by affecting upstream folate handling and downstream B12-dependent recycling. Whether this co-occurrence contributes to phenotypic variability among compound heterozygous individuals warrants investigation.

### Translational Context: Human Base-Editing Programs in Other Targets

| Program | Company | Target | Delivery | Status |
|---------|---------|--------|----------|--------|
| VERVE-102 | Verve Therapeutics | PCSK9 (single base edit) | IV LNP | Phase 1b |
| Beam-101 | Beam Therapeutics | Sickle cell (base edit) | Ex vivo | Phase I/II |
| Casgevy | Vertex/CRISPR | Sickle cell / beta-thal | Ex vivo | **FDA Approved (2023)** |
| YOLT-101 | YolTech | PCSK9 (next-gen ABE) | IV LNP | Phase 1 (2025) |

At the sequence level, C677T appears theoretically compatible with **adenine base editing (ABE)** logic, which converts A-T base pairs to G-C base pairs. However, edit window, guide design, bystander edits, tissue targeting, and functional rescue would all need empirical validation.

---

## Safety Architecture: What If Something Goes Wrong?

Any future gene-editing program would need to address safety. The following are examples of safety layers used or discussed in broader editing research that could be relevant to MTHFR correction if it advances to experimental testing. Each would require target-specific validation.

| Layer | Concept | Precedent |
|-------|---------|-----------|
| **1. Inherent design** | Base editing avoids DNA double-strand breaks; LNP-delivered mRNA degrades within 48 hours | Used in VERVE-102, BEAM-101 |
| **2. Off-target screening** | Cas-OFFinder, GUIDE-seq, whole-genome sequencing before any intervention | Standard in current editing programs |
| **3. Anti-CRISPR proteins** | Emergency neutralization of residual editor activity (e.g., AcrIIA4) | Demonstrated in preclinical models |
| **4. Reversible epigenetic editing** | CRISPRoff as a non-permanent first step before committing to base editing | In early translational development |
| **5. Post-intervention monitoring** | Biomarker panels, off-target sequencing at defined timepoints | Standard clinical protocol design |

**Key point:** At the target base, a candidate correction strategy would aim to restore the reference allele sequence rather than introduce a novel coding change.

See the [full research paper](docs/RESEARCH_PAPER_DRAFT.md#47-safety-architecture-built-in-safeguards-and-emergency-off-switches) for complete safety analysis.

---

## Quick Start: Replicate This Research

The workflow is designed to be reproducible with publicly accessible tools.

### Option A: Molecular Dynamics on Google Colab (One Click)

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/DSMPromo/mthfr-target-validation/blob/main/MTHFR_MD_Simulation.ipynb)

Run 100ns molecular dynamics simulations comparing WT and compound heterozygous MTHFR dimers. Requires a Google account and Colab GPU runtime (free T4 or paid A100).

### Option B: AlphaFold Analysis on Google Colab

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/DSMPromo/mthfr-target-validation/blob/main/MTHFR_AlphaFold_Analyzer.ipynb)

Run the automated analysis pipeline on pre-computed AlphaFold 3 and Boltz-2 results. Generates metrics, PAE heatmaps, comparison charts, and HTML report.

### Option C: AlphaFold Server (Run Your Own Predictions)
1. Go to [alphafoldserver.com](https://alphafoldserver.com) and sign in
2. Click **"Upload JSON"**
3. Upload `alphafold/jobs/json_all/ALL_30_JOBS.json` from this repo
4. Submit all 30 jobs as drafts, then confirm each one

### Option D: Local Setup
Follow the detailed step-by-step instructions in [`alphafold/jobs/submission_plan.md`](alphafold/jobs/submission_plan.md)

### After Jobs Complete
```bash
# Clone the repo
git clone https://github.com/DSMPromo/mthfr-target-validation.git
cd mthfr-target-validation

# Set up environment
python3 -m venv .venv && source .venv/bin/activate
pip install matplotlib numpy

# Drop your AlphaFold results into alphafold/results_all/

# Run the analysis
python analyze.py

# Open the report
open analysis/outputs/report.html
```

---

## Development Roadmap

From computational prioritization to staged experimental follow-up:

| Phase | Focus | Timeline | Cost | Status |
|-------|-------|----------|------|--------|
| **1** | Computational hypothesis prioritization: AlphaFold 3 (64 predictions, 10 seeds) + Boltz-2 (4 substrate/inhibitor), RMSD validation vs PDB 6FCX, statistical analysis | **Done** | **$0** | ✅ **Complete** |
| **2a** | Molecular dynamics simulations: OpenMM/Amber14, WT dimer vs compound dimer, RMSD/RMSF/flexibility analysis | **In progress** | **$0** | 🔄 **Running** |
| **2b** | Preprint (bioRxiv submitted), researcher outreach, extended structural benchmarking | 1-3 months | ~$500 | In progress |
| **3** | Biochemical validation: expression of wild-type and selected variant proteins, dimer stability assays, FAD and THF interaction assays | 6-12 months | $50K-150K | Requires lab partner |
| **4** | Exploratory biomarker correlation studies, including homocysteine, methylation panels, retinal OCT, and BH4-related readouts in well-defined compound heterozygous cohorts | 12-18 months | $100K-300K | Requires clinical collaborator |
| **5** | If biochemical and biomarker findings support continued investigation: sequence-level guide assessment, edit-window analysis, bystander-risk assessment, and cell-based feasibility studies | 18-24 months | $200K-500K | Requires gene-editing partner |
| **6** | If cell-based feasibility and safety data support further progression: preclinical animal studies using appropriate Mthfr-relevant models and multi-organ outcome measures | 24-36 months | $500K-2M | Requires institutional partner |
| **7** | If preclinical data support advancement: IND-enabling planning, regulatory strategy development, and early clinical trial concept design | 3-7+ years | $5M+ | Requires biotech or pharma partner |

> Advancement between phases is conditional and depends on supportive data at each prior stage. This roadmap is a planning framework, not a claim of therapeutic readiness.

**Where we are now:** Phase 1 complete. The computational hypothesis prioritization is documented and replicable. Built for $0.

**The next step:** Identify a researcher with biochemical or structural biology capabilities to evaluate whether these computational observations translate to measurable differences in experimental systems.

---

## Project Structure

```
mthfr-target-validation/
|-- README.md                           <-- You are here
|-- analyze.py                          <-- Automated analysis pipeline (AlphaFold 3 + Boltz-2)
|-- generate_figures.py                 <-- Publication figure generator
|-- validate_language.py                <-- Language compliance scanner (21 rules)
|-- requirements.txt                    <-- Python dependencies
|-- MTHFR_AlphaFold_Analyzer.ipynb      <-- Google Colab: AlphaFold analysis
|-- MTHFR_MD_Simulation.ipynb           <-- Google Colab: 100ns MD simulations
|-- run_md.py                           <-- Local MD pipeline (OpenMM/M1 Max)
|-- validate_structure.py               <-- RMSD, Cohen's d, Bonferroni validation
|-- MTHFR-Research-Findings.zip         <-- Downloadable archive of all findings
|-- DISCLAIMER.md / CONTRIBUTING.md / LICENSE
|
|-- docs/
|   |-- MTHFR_Master_Document.md                 <-- Complete research document (v6.0)
|   |-- MTHFR_Research_Paper.pdf                 <-- Publication-ready PDF with figures
|   |-- RESEARCH_PAPER_DRAFT.md                  <-- Publication-quality paper draft
|
|-- sequences/                          <-- Verified UniProt P42898 sequences
|   |-- MTHFR_wildtype.fasta            <-- 656aa, A at 222, E at 429
|   |-- MTHFR_C677T.fasta              <-- A222V (catalytic domain variant)
|   |-- MTHFR_A1298C.fasta             <-- E429A (regulatory domain variant)
|   |-- MTHFR_compound.fasta           <-- Both mutations
|
|-- alphafold/
|   |-- jobs/
|   |   |-- submission_plan.md          <-- 34 jobs, step-by-step
|   |   |-- json/                       <-- Original 12 job JSONs
|   |   |-- json_all/                   <-- All 30 AlphaFold Server JSONs + ALL_34_JOBS.json
|   |-- results_all/                    <-- All 64 results (unified folder)
|       |-- wt_mono_run1..10/            <-- WT monomer, 10 seeds
|       |-- wt_dimer_run1..10/           <-- WT dimer, 10 seeds
|       |-- c677t_mono_run1..10/         <-- C677T monomer, 10 seeds
|       |-- c677t_dimer_run1..10/        <-- C677T dimer, 10 seeds
|       |-- a1298c_mono_run1..10/        <-- A1298C monomer, 10 seeds
|       |-- compound_dimer_run1..10/     <-- Compound dimer, 10 seeds
|       |-- wt_dimer_thf/              <-- WT + THF (Boltz-2)
|       |-- c677t_dimer_thf/           <-- C677T + THF (Boltz-2)
|       |-- compound_dimer_thf/        <-- Compound + THF (Boltz-2)
|       |-- wt_dimer_sam/              <-- WT + SAM (Boltz-2)
|
|-- analysis/
|   |-- analysis_workflow.md            <-- PyMOL commands, analysis steps
|   |-- metrics_template.csv            <-- Recording template
|   |-- outputs/
|       |-- metrics.csv                 <-- All 64 jobs: pTM, ipTM, pLDDT, ligand scores
|       |-- report.html                 <-- Interactive HTML report
|       |-- charts/                     <-- ipTM and pTM comparison bar charts
|       |-- pae_plots/                  <-- 64 PAE heatmaps
|       |-- figures/                    <-- 3D structures, dashboard, pLDDT comparison
|       |-- pymol_scripts/              <-- Ready-to-use PyMOL analysis scripts
|
|-- outreach/
|   |-- target_researchers.md           <-- Collaboration targets
|   |-- email_template.md              <-- Outreach template
|
|-- .github/ISSUE_TEMPLATE/             <-- Collaboration issue template
```

---

## Verified References

All references have been verified against PubMed and primary sources (March 2026).

### Vision and Retinal Neurodegeneration
1. Markand S, et al. (2015). *Invest Ophthalmol Vis Sci*, 56(4):2684-2695. [PMID: 25766590](https://pubmed.ncbi.nlm.nih.gov/25766590/) -- Mthfr+/- mice: ~20% ganglion cell loss, retinal vasculopathy
2. Hergert RM, et al. (2022). *J Neurol*, 269:4571-4577. [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC8940754/) -- Case report: betaine treatment was associated with visual recovery in severe MTHFR deficiency
3. Jiang H, Liu Z, et al. (2023). *Clin Ophthalmol*, 17:1035-1043. [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC10106310/) -- Retinal perfusion improved by MTHFR genotype

### Autism and Cerebral Folate Deficiency
4. Pu D, et al. (2020). *BMC Pediatr*, 20:471. [PMID: 32972375](https://pubmed.ncbi.nlm.nih.gov/32972375/) -- Meta-analysis: C677T significantly associated with ASD
5. Frye RE, et al. (2018). *Mol Psychiatry*, 23:631-636. [Nature](https://www.nature.com/articles/mp2016168) -- Folinic acid RCT: improved verbal communication in ASD

### Neuropsychiatric Disorders
6. Wan L, et al. (2018). *Transl Psychiatry*, 8:242. [PMID: 30397195](https://pubmed.ncbi.nlm.nih.gov/30397195/) -- MTHFR and schizophrenia, depression, bipolar, ADHD
7. Levin BL, Varga E. (2016). *J Genet Couns*, 25(5):901-911. [Wiley](https://onlinelibrary.wiley.com/doi/full/10.1007/s10897-016-9956-7) -- MTHFR genetic counseling and clinical evidence

### Structural Biology
8. Froese DS, et al. (2018). *Nat Commun*, 9:4735. [Nature](https://www.nature.com/articles/s41467-018-04735-2) -- First human MTHFR crystal structure (PDB: [6FCX](https://www.rcsb.org/structure/6FCX))
9. Froese DS, et al. (2024). *Nat Commun*. [Nature](https://www.nature.com/articles/s41467-024-47174-y) -- Cryo-EM: MTHFR + SAM allosteric regulation (PDB: [8QA5](https://www.rcsb.org/structure/8QA5), [8QA6](https://www.rcsb.org/structure/8QA6))
10. Pejchal R, et al. (2006). *Biochemistry*, 45(15):4808-4818. [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC1868400/) -- A222V structural perturbation: helix alpha-5 displacement disrupts FAD binding

### B-Vitamins and Cardiovascular
11. McNulty H, et al. (2006). *Circulation*, 113(1):74-80. [AHA](https://www.ahajournals.org/doi/10.1161/circulationaha.105.580332) -- Riboflavin lowers homocysteine in 677TT homozygotes

### Epigenetics and DNA Methylation
12. Karahan G, Chan D, Shirane K, et al. (2021). *Development*, 148(13):dev199492. [PMID: 34128976](https://pubmed.ncbi.nlm.nih.gov/34128976/) -- Paternal MTHFR deficiency causes transgenerational loss of sperm DNA methylation
13. Pentieva K, et al. (2020). *Biochimie*, 173:17-26. [ScienceDirect](https://www.sciencedirect.com/science/article/abs/pii/S0300908420300778) -- Riboflavin alters global/gene-specific DNA methylation in 677TT adults

### Pregnancy and Reproductive Health
14. Xia XP, et al. (2013). *Hypertension Research*. [Nature](https://www.nature.com/articles/hr2012117) -- Meta-analysis: C677T and preeclampsia risk
15. Kumar A, et al. (2024). *BMC Pregnancy Childbirth*. [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC11066823/) -- A1298C and recurrent pregnancy loss

### Comprehensive Reviews
16. MTHFR Gene Polymorphisms: A Single Gene with Wide-Ranging Clinical Implications. (2025). *Genes*, 16(4):441. [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC12027316/) -- Comprehensive 2025 review

### Computational Methods
17. Abramson J, et al. (2024). *Nature*, 630:493-500. [Nature](https://www.nature.com/articles/s41586-024-07487-w) -- AlphaFold 3: accurate structure prediction of biomolecular interactions

### Gene Therapy and Base Editing
18. YOLT-101 Phase 1 Clinical Data. (2025). [YolTech](https://www.yoltx.com/news/press-release/95) -- Next-gen adenine base editor, Phase 1 clinical data

### Key Database Links
[UniProt P42898](https://www.uniprot.org/uniprotkb/P42898/entry) | [PDB 6FCX](https://www.rcsb.org/structure/6FCX) | [PDB 8QA5](https://www.rcsb.org/structure/8QA5) | [PDB 8QA6](https://www.rcsb.org/structure/8QA6) | [ClinVar rs1801133](https://www.ncbi.nlm.nih.gov/clinvar/?term=rs1801133) | [ClinVar rs1801131](https://www.ncbi.nlm.nih.gov/clinvar/?term=rs1801131)

---

## How You Can Help

See [CONTRIBUTING.md](CONTRIBUTING.md) for detailed instructions.

**If you're a scientist:** Our structural and literature synthesis identifies a bounded experimental agenda centered on dimer-level effects, retinal biomarkers, and neuropsychiatric biochemistry. We are choosing to focus where the existing literature is most developed. [Open an issue](../../issues/new?template=collaboration.md) or email me.

**What we're looking for:**
- **Ophthalmologists / Retinal researchers:** Existing Mthfr+/- mouse models provide one possible system for testing whether correction strategies alter retinal biomarkers. Delivery route, tissue targeting, and edit feasibility would need separate evaluation
- **Psychiatrists / Neuroscientists:** Could A1298C-related one-carbon pathway alterations be a contributing factor in defined subgroups with biochemical abnormalities? The mechanism and subgroup specificity remain incompletely defined
- **Structural biologists:** Are the AlphaFold/Boltz-2 predictions biologically sensible? Cryo-EM of C677T/A1298C variants with FAD would provide orthogonal validation
- **Gene editing researchers:** Evaluate whether C677T is technically tractable for guide design, edit-window fit, bystander risk assessment, and cell-based rescue studies
- **Replication:** Independent researchers reproducing and extending the computational analysis

**If you're not a scientist:** Share this project with researchers who might find it useful. If you are a patient or family member, discuss any testing or interpretation questions with a qualified clinician.

---

## License

[Creative Commons Attribution-NonCommercial-ShareAlike 4.0](https://creativecommons.org/licenses/by-nc-sa/4.0/)

### Reuse and Licensing Summary

| You Want To... | Allowed? | Details |
|----------------|----------|---------|
| Read, download, and share this research | ✅ Yes | No permission needed |
| Replicate the AlphaFold analysis in your lab | ✅ Yes | That's why we built it |
| Cite this work in your own papers | ✅ Yes | Credit Igor Mihaljko as original author |
| Build on this research for academic/educational purposes | ✅ Yes | Your work must also be shared openly (ShareAlike) |
| Use this work in commercial therapeutic or product development | ❌ Not without permission | Contact igor@dsm.promo for commercial licensing |
| Take this work and make it closed/proprietary | ❌ No | ShareAlike requires derivative work to remain open |

For academic and noncommercial use, this work can be shared and built upon under the terms of the CC BY-NC-SA 4.0 license. For any use outside that license, please contact the author.

---

## Contact

**Igor Mihaljko** | CEO and Founder, DSM.Promo | Chicago, IL | igor@dsm.promo | [ORCID: 0009-0000-1408-1065](https://orcid.org/0009-0000-1408-1065)

---

*I'm not a doctor. I'm not a biohacker. I'm a cybersecurity guy and cloud architect who looked at his DNA results, started asking questions, and used AI tools to organize what I found. I don't have a lab. I don't have a PhD. I built this project to hand it to the real scientists who can take it further.*

*The structural workflow is reproducible within the documented pipeline, the cited references have been checked against primary or database sources, and the analysis pipeline is open source. I am not qualified to determine whether these computational differences are biologically meaningful, but I have tried to make the question easier for qualified researchers to test.*

*If you are a clinician or researcher, I hope this project helps clarify a tractable experimental path. If you are a patient or family member, this project should be read as a research document, not as medical guidance. Variant status needs clinical interpretation in context with symptoms, labs, and the broader medical picture.*

*My goal is to help move a testable research question into the hands of qualified scientists who can evaluate it rigorously.*
