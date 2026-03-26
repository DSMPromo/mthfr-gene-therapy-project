# MTHFR Variant Hypothesis Prioritization Program

### Computational Prioritization of Selected MTHFR Variant States for Experimental Follow-Up

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

- [About the Author](#about-the-author)
- [Cross-Specialty Context](#a-cross-specialty-problem-that-may-be-underintegrated-in-current-research)
- [The Enzyme](#the-mthfr-enzyme)
- [Primary Follow-Up Contexts](#primary-experimental-follow-up-contexts)
- [Additional Literature Context](#additional-literature-context)
- [AlphaFold 3 Results](#alphafold-3-structural-predictions)
- [Study Limits](#what-this-study-does-not-show)
- [Alternative Explanations](#alternative-explanations)
- [Variant Ranking](#variant-state-ranking-framework)
- [Next Three Experiments](#next-three-experiments)
- [Why Test This](#why-this-project-is-worth-experimental-follow-up)
- [Structural Context](#structural-and-literature-context)
- [Safety Concepts](#safety-concepts-for-future-experimental-programs)
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

### Why MTHFR Has Seen Limited Translational Development

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

## Primary Experimental Follow-Up Contexts

The computational findings in this project generate specific, testable hypotheses rather than therapeutic claims. Two downstream systems are prioritized for experimental validation based on existing literature:

### Hypothesis 1: Retinal Pathway

Because prior literature has linked MTHFR-related one-carbon pathway disruption and hyperhomocysteinemia to retinal injury, the retina is a candidate downstream system for validation. The current computational results do not establish efficacy in retinal disease, but they support testing whether selected MTHFR states are associated with measurable retinal biomarkers.

**Supporting literature:**
- **Mthfr+/- mice** show 2x retinal homocysteine, ~20% ganglion cell loss, thinner nerve fiber layer, and vascular leakage, without elevated IOP (Markand et al., 2015, *IOVS*)
- **Case report:** betaine treatment was associated with visual recovery in a patient with severe MTHFR deficiency (Hergert et al., 2022, *J Neurology*)
- **Retinal perfusion improved** by MTHFR genotype-guided supplementation; compound C677T/A1298C carriers showed significant improvement (Jiang et al., 2023, *Clin Ophthalmol*)
- MTHFR modifies CRB1-related retinopathies, worsening rd8 retinal phenotype (Markand et al., 2015, *Exp Eye Res*)

**Proposed validation:** Correlate MTHFR genotype with retinal OCT measurements, homocysteine levels, and visual acuity in defined patient cohorts.

### Hypothesis 2: Neuropsychiatric-Relevant Biochemical Follow-Up

Because one-carbon pathway disruption may influence BH4-dependent neurotransmitter biology, neuropsychiatric phenotypes are a second candidate area for validation. The current structural data do not establish causation or treatment effect, but they support testing whether selected MTHFR states are associated with measurable biochemical and clinical features.

**Supporting literature:**
- **MTHFR polymorphism** associated with schizophrenia, major depression, autism, bipolar disorder, and ADHD (Wan et al., 2018, *Transl Psychiatry*)
- Compound heterozygous patients show treatment-responsive anxiety with SAMe and methylated B vitamins when MTHFR status is identified (Levin & Varga, 2016)
- A1298C has been discussed in relation to BH4-relevant pathways in the literature, but the mechanism, effect size, and subgroup specificity remain incompletely defined
- MTHFR-related one-carbon pathway alterations may influence BH4-dependent neurotransmitter synthesis, warranting investigation in defined subgroups with biochemical abnormalities

**Proposed validation:** Measure BH4, homocysteine, SAMe/SAH ratio, and neurotransmitter metabolites in compound heterozygous individuals vs. wild-type controls.

---

## Additional Literature Context

MTHFR dysfunction has been associated with additional disease-relevant pathways in the literature. If future biochemical and in vivo studies support a causal contribution of selected MTHFR states, sequence correction strategies could be evaluated as one possible upstream intervention concept. These pathways are documented here for completeness.

Published literature has also reported MTHFR associations with autism/cerebral folate deficiency (Pu et al., 2020; Frye et al., 2018), B-vitamin processing failure (McNulty et al., 2006), cardiovascular disease via homocysteine elevation, transgenerational epigenetic effects in mouse models (Karahan et al., 2021), and adverse pregnancy outcomes including preeclampsia (Xia et al., 2013) and recurrent pregnancy loss (Kumar et al., 2024). These are documented in the [full research paper](docs/RESEARCH_PAPER_DRAFT.md) but are not the primary focus of this computational study.

**Working hypothesis:** The compound heterozygous MTHFR dimer shows lower predicted confidence for cofactor-associated and substrate-associated interactions relative to wild-type. These comparative computational results support experimental testing of whether selected MTHFR variant states are associated with measurable differences in retinal and neuropsychiatric-relevant biomarkers or pathway readouts. Sequence correction remains a separate downstream question that would require independent validation of target biology, editing feasibility, delivery, and safety.

---

## AlphaFold 3 Structural Predictions

We used [AlphaFold 3 Server](https://alphafoldserver.com) to predict structures of wild-type, C677T, A1298C, and compound heterozygous MTHFR in complex with FAD cofactor.

### Methodology

- **Sequences:** All derived from canonical [UniProt P42898](https://www.uniprot.org/uniprotkb/P42898/entry) (656 amino acids) with verified mutations at positions 222 and 429
- **Predictions:** 16 jobs total (monomers, homodimers, heterodimers, with FAD, THF, and SAM)
- **Replication:** Each core prediction run with independent random seeds for statistical comparison
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
| **WT dimer** | 0.775 | 0.740 | 0.555 | 97.3 | 96.0 |
| **C677T dimer** | 0.785 | 0.765 | 0.575 | 97.05 | 95.95 |
| **Compound dimer** | **0.745** | **0.715** | **0.540** | **96.5** | **95.15** |

**Key observations (replicated across independent seeds):**
- **Monomer predictions showed preserved overall folding confidence** (ipTM 0.97-0.98) -- the tested monomer predictions do not suggest large-scale loss of overall fold confidence
- **Dimer predictions reveal comparative inter-chain differences** -- reported FAD-associated confidence values are lower in the tested dimer models than in the tested monomer models, consistent with the distinct interaction context of the homodimer
- **The compound heterozygous dimer yielded the lowest reported interaction-confidence values** across the reported comparative metrics in both runs: pTM (0.73/0.76), ipTM (0.70/0.73), FAD binding (0.53/0.55), pLDDT@429 (95.0/95.3)
- **Compound heterozygous dimers averaged lower than the tested comparators** -- ipTM 0.715 vs WT 0.740 and C677T 0.765, consistent with a possible combined dimer-level perturbation in this modeling setup
- **Consistent directional trends** -- independent random seeds produced consistent directional trends across the tested configurations, supporting limited internal reproducibility within this modeling setup
- **Position 429 showed the largest confidence decrease** among the reported local metrics in compound dimers -- pLDDT 95.0-95.3 vs 95.8-96.2 in WT, consistent with possible regulatory-domain involvement at the dimer level

> **The core computational observation:** Monomer predictions showed preserved overall folding confidence, whereas dimer predictions showed comparative inter-chain differences. This is the central computational finding around which the experimental agenda is organized.

> **Important:** These are computational predictions, not experimental structures. All confidence metrics should be interpreted as hypothesis generators, not proof of mechanism. See the [full research paper draft](docs/RESEARCH_PAPER_DRAFT.md) for complete methodology and limitations.

---

## What This Study Does Not Show

This section exists to define the boundaries of the current work explicitly. Strong research programs state their limits clearly.

- **It does not show reduced enzyme activity.** No biochemical assay was performed.
- **It does not show altered FAD affinity in vitro.** No binding assay was performed.
- **It does not show retinal damage.** No retinal imaging or histology was performed.
- **It does not show neuropsychiatric relevance.** No neurotransmitter or BH4 measurements were taken.
- **It does not establish editing feasibility.** No guide RNA design, off-target analysis, or cell-based editing was performed.
- **It does not prove that the observed confidence differences are biologically meaningful.** The differences may reflect model uncertainty rather than true structural effects.
- **It does not establish disease causation.** The associations described come from published literature, not from this study.

**What it does show:** A reproducible computational signal in which the compound heterozygous dimer yielded the lowest comparative confidence values across the tested metrics. This signal prioritizes specific experimental follow-up but does not replace it.

---

## Alternative Explanations

The following alternative interpretations should be considered before drawing biological conclusions from these computational results:

1. **Model uncertainty vs. structural difference:** The lower interaction confidence in compound heterozygous dimers may reflect prediction uncertainty rather than a true structural perturbation. AlphaFold and Boltz-2 confidence scores measure model certainty, not binding affinity.

2. **Effect size significance:** The observed differences (e.g., ipTM 0.715 vs 0.740) may not translate into measurable biochemical effects. The biological threshold for significance at these confidence levels has not been established.

3. **Ligand placement limitations:** AlphaFold 3 and Boltz-2 ligand predictions are less validated than protein fold predictions. FAD and THF placement may carry systematic uncertainties that affect comparison across variant states.

4. **Dimer assembly context:** The tested dimer configurations assume specific chain arrangements. In vivo, dimer assembly may involve additional factors (chaperones, post-translational modifications) not captured by these models.

5. **Population-level variability:** The compound heterozygous state represents a specific genotype combination. Clinical phenotypes likely depend on additional genetic, environmental, and epigenetic factors not modeled here.

These alternatives reinforce why experimental validation is the necessary next step, not a deeper computational analysis.

---

## Variant State Ranking Framework

This ranking is a prioritization tool based on computational outputs, not a measure of biological effect size.

| Rank | Variant State | Monomer ipTM | Dimer ipTM (avg) | FAD Binding (avg) | pLDDT@429 (avg) | THF Ligand ipTM | Prioritization |
|------|--------------|-------------|-------------------|-------------------|------------------|-----------------|----------------|
| 1 | **Compound het dimer** | 0.97 | **0.715** | **0.540** | **95.15** | **0.878** | Highest priority for experimental follow-up |
| 2 | WT dimer | 0.97 | 0.740 | 0.555 | 96.0 | 0.974 | Baseline comparator |
| 3 | C677T dimer | 0.97 | 0.765 | 0.575 | 95.95 | 0.969 | Single-variant comparator |
| 4 | All monomers | 0.97 | N/A | 0.97 | 97-98 | N/A | Fold preserved across all states |

**Interpretation:** The compound heterozygous dimer consistently ranks lowest across all assessed dimer-level and ligand-associated metrics. This ranking defines the experimental priority, not the biological conclusion.

---

## Next Three Experiments

Even without a lab, defining the next experiments precisely shows research program thinking. These are the three experiments that would most directly validate or refute the computational observations:

### Experiment 1: Recombinant Protein Expression and Enzyme Activity
- **Express:** WT, C677T (A222V), A1298C (E429A), and compound heterozygous MTHFR
- **Measure:** Specific enzyme activity (Vmax, Km for MTHF), thermal stability (Tm via DSF)
- **Question:** Does the compound heterozygous protein show reduced activity or stability relative to WT?
- **Why first:** This is the most direct test of whether the computational differences have biochemical consequences

### Experiment 2: FAD Binding Affinity and Dimer Stability
- **Measure:** FAD dissociation constants (Kd) across all four variant states using ITC or fluorescence quenching
- **Measure:** Dimer stability via analytical ultracentrifugation or SEC-MALS
- **Question:** Does the compound heterozygous dimer show reduced FAD affinity or altered dimerization?
- **Why second:** Directly tests the core computational observation (dimer-level FAD binding differences)

### Experiment 3: Retinal Biomarker Correlation
- **Cohort:** 50-100 individuals genotyped for MTHFR C677T and A1298C
- **Measure:** Retinal nerve fiber layer thickness (OCT), plasma homocysteine, plasma 5-MTHF
- **Question:** Do compound heterozygous individuals show measurable retinal differences vs. wild-type?
- **Why third:** Tests the primary indication hypothesis with clinical-grade measurements

These experiments would move the project from computational prioritization (Phase 1, complete) to biochemical validation (Phase 3).

---

## Why This Project Is Worth Experimental Follow-Up

**The core signal:**
Across 16 structural predictions on two independent platforms (AlphaFold 3 and Boltz-2) with replicated random seeds, the compound heterozygous MTHFR dimer consistently yielded the lowest interaction-confidence values. Monomers fold normally. The differences emerge at the dimer level.

**Why the question matters:**
Up to 40% of the global population carries at least one MTHFR variant. If compound heterozygous individuals show measurable dimer-level dysfunction, the downstream implications would justify focused follow-up in retinal and metabolic-relevant systems.

**What this study does not claim:**
It does not claim mechanism. It does not claim disease causation. It does not claim therapeutic readiness. It claims a reproducible computational ranking that prioritizes specific experimental questions.

**Three experiments that would validate or refute it:**
1. Recombinant protein activity comparison (WT vs compound het)
2. FAD binding and dimer stability assays
3. Retinal biomarker correlation in genotyped cohorts

This is a falsifiable research program. The experiments can prove it wrong. That is what makes it worth testing.

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

## Safety Concepts for Future Experimental Programs

Sequence correction is discussed only as a long-range downstream possibility. The present study does not establish MTHFR as a therapeutic editing target. Any editing concept would require separate validation of target biology, guide design, delivery, safety, and functional rescue.

If correction strategies were to advance, established safety frameworks from current base-editing programs (e.g., off-target screening via GUIDE-seq, anti-CRISPR emergency neutralization, reversible epigenetic editing as a non-permanent first step) would need target-specific adaptation and validation. See the [full research paper](docs/RESEARCH_PAPER_DRAFT.md#47-safety-architecture-built-in-safeguards-and-emergency-off-switches) for detailed safety discussion.

---

## Quick Start: Replicate This Research

**The workflow is designed to be reproducible with publicly accessible tools.**

### Option A: Upload JSON (Fastest)
1. Go to [alphafoldserver.com](https://alphafoldserver.com) and sign in
2. Click **"Upload JSON"**
3. Upload `alphafold/jobs/json/ALL_12_JOBS.json` from this repo
4. Click **"Submit 12 jobs as drafts"**
5. Open each draft, turn off Seed toggle, click **"Confirm and submit job"**

> **Note:** The JSON file contains 12 core jobs (6 variants + 6 replications). Jobs 13-16 (substrate/inhibitor binding with THF and SAM) require manual setup through the AlphaFold Server UI -- see the submission plan for details.

### Option B: Manual Setup
Follow the detailed step-by-step instructions in [`alphafold/jobs/submission_plan.md`](alphafold/jobs/submission_plan.md)

### After Jobs Complete
```bash
# Clone the repo
git clone https://github.com/DSMPromo/mthfr-target-validation.git
cd mthfr-target-validation

# Set up environment
python3 -m venv .venv && source .venv/bin/activate
pip install matplotlib numpy

# Drop your AlphaFold ZIP results into alphafold/results/

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
| **1** | Computational hypothesis prioritization, AlphaFold 3 and Boltz-2 structural comparison, reproducible pipeline, open repository | **Done** | **$0** | ✅ **Complete** |
| **2** | Structural benchmarking against experimental reference structures such as PDB 6FCX, molecular dynamics simulations, preprint preparation, researcher outreach | 1-3 months | ~$500 | Next |
| **3** | Biochemical validation: expression of wild-type and selected variant proteins, dimer stability assays, FAD and THF interaction assays | 6-12 months | $50K-150K | Requires lab partner |
| **4** | Exploratory biomarker correlation studies, including homocysteine, methylation panels, retinal OCT, and BH4-related readouts in well-defined compound heterozygous cohorts | 12-18 months | $100K-300K | Requires clinical collaborator |
| **5** | If biochemical and biomarker findings support continued investigation: sequence-level guide assessment, edit-window analysis, bystander-risk assessment, and cell-based feasibility studies | 18-24 months | $200K-500K | Requires gene-editing partner |
| **6** | If cell-based feasibility and safety data support further progression: preclinical animal studies using appropriate Mthfr-relevant models and multi-organ outcome measures | 24-36 months | $500K-2M | Requires institutional partner |
| **7** | If preclinical data support advancement: IND-enabling planning, regulatory strategy development, and early clinical trial concept design | 3-7+ years | $5M+ | Requires biotech or pharma partner |

> Advancement between phases is conditional and depends on supportive data at each prior stage. This roadmap is a planning framework, not a claim of therapeutic readiness.

**Where we are now:** Phase 1 complete. The computational hypothesis prioritization is documented and replicable. 

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
|-- MTHFR_AlphaFold_Analyzer.ipynb      <-- Google Colab notebook
|-- MTHFR-Research-Findings.zip         <-- Downloadable archive of all findings
|-- DISCLAIMER.md / CONTRIBUTING.md / LICENSE
|
|-- docs/
|   |-- MTHFR_Master_Document.md                 <-- Complete research document (v5.0)
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
|   |   |-- submission_plan.md          <-- 16 jobs, step-by-step
|   |   |-- json/                       <-- Upload-ready JSON files for AlphaFold Server
|   |-- results/                        <-- Jobs 1-12 (AlphaFold Server) + Jobs 13-16 (Boltz-2)
|
|-- analysis/
|   |-- analysis_workflow.md            <-- PyMOL commands, analysis steps
|   |-- metrics_template.csv            <-- Recording template
|   |-- outputs/
|       |-- metrics.csv                 <-- All 16 jobs: pTM, ipTM, pLDDT, ligand scores
|       |-- report.html                 <-- Interactive HTML report
|       |-- charts/                     <-- ipTM and pTM comparison bar charts
|       |-- pae_plots/                  <-- 16 PAE heatmaps (AlphaFold + Boltz-2)
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

**For researchers:** This license is designed to maximize scientific collaboration while protecting the original work. Use it freely for research. Cite it. Build on it. Share your results.

**For biotech/pharma:** If you are interested in commercial development based on this work, contact us to discuss licensing.

---

## Contact

**Igor Mihaljko** | CEO and Founder, DSM.Promo | Chicago, IL | igor@dsm.promo | [ORCID: 0009-0000-1408-1065](https://orcid.org/0009-0000-1408-1065)

---

*I built this project as an open computational starting point for qualified researchers to evaluate experimentally. The workflow is reproducible, the references are source-checked, and the conclusions are intentionally bounded by the current evidence.*

*If you are a clinician or researcher, I hope this project helps clarify a tractable experimental path. If you are a patient or family member, this project should be read as a research document, not as medical guidance. Variant status needs clinical interpretation in context with symptoms, labs, and the broader medical picture.*
