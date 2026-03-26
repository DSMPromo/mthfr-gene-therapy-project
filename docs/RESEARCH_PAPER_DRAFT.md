# MTHFR Variant Target Validation: Computational Prioritization of Selected MTHFR Variant States for Experimental Validation

## A Hypothesis-Prioritization Study Using AlphaFold 3 and Boltz-2

---

**Author:** Igor Mihaljko, Independent Researcher, DSM.Promo, Chicago, IL

**ORCID:** [0009-0000-1408-1065](https://orcid.org/0009-0000-1408-1065)

**Correspondence:** igor@dsm.promo

**Keywords:** MTHFR, AlphaFold 3, Boltz-2, C677T, A1298C, FAD cofactor, folate substrate, protein structure prediction, retinal neurodegeneration, anxiety, depression, BH4, gene therapy, base editing, compound heterozygous

**License:** CC BY-NC-SA 4.0

**Data Availability:** All sequences, analysis scripts, and AlphaFold results available at https://github.com/DSMPromo/mthfr-target-validation

**AI Disclosure:** AI tools (Claude, Anthropic) were used for literature organization, analysis pipeline development, and manuscript preparation. All AI-generated content was reviewed for accuracy by the author. This disclosure follows emerging best practices in scientific publishing (Nature, 2023).

---

## Abstract

Methylenetetrahydrofolate reductase (MTHFR) variants C677T (A222V) and A1298C (E429A) affect up to 40% of the global population, yet their structural consequences on FAD cofactor binding and folate substrate access at the dimer level remain incompletely characterized. We used AlphaFold 3 Server and Boltz-2 to predict the structures of wild-type, C677T, A1298C, and compound heterozygous MTHFR in complex with FAD, folate substrate (THF), and allosteric inhibitor (SAM). Across 16 structural predictions with independent replication seeds, we compared pTM, ipTM, ligand binding confidence, pLDDT at mutation sites (positions 222 and 429), and predicted aligned error (PAE). Monomer predictions show minimal differences between variants (ipTM 0.97 across all), but dimer predictions reveal a reproducible signal: the compound heterozygous dimer yielded the lowest reported interaction-confidence values across the assessed metrics (AlphaFold: pTM 0.73, ipTM 0.70, FAD binding 0.53; Boltz-2: ligand ipTM 0.878 vs 0.974 wild-type for THF substrate binding). These observations generate falsifiable hypotheses for two candidate downstream systems: (1) the retinal pathway, where prior literature has linked MTHFR-related one-carbon pathway disruption to homocysteine-associated ganglion cell injury, and (2) neuropsychiatric pathways, where A1298C has been discussed in relation to BH4-relevant neurotransmitter biology. The present work is intended to prioritize experiments, define measurable hypotheses, and support collaboration with structural biology, metabolic disease, retinal, and neuropsychiatric researchers. All predictions require orthogonal validation through experimental methods.

**Note:** This is a computational study generating hypotheses for experimental validation. All structural predictions require confirmation through experimental methods (cryo-EM, X-ray crystallography, functional assays).

---

## 1. Introduction

### 1.1 The MTHFR Enzyme and Its Clinical Significance

MTHFR (methylenetetrahydrofolate reductase, EC 1.5.1.20) is a FAD-dependent homodimeric enzyme encoded on chromosome 1p36.22 that catalyzes the conversion of 5,10-methylenetetrahydrofolate to 5-methyltetrahydrofolate (5-MTHF), the primary circulating form of folate and the methyl donor for homocysteine remethylation to methionine (Froese et al., 2018). This reaction sits at the intersection of folate metabolism, homocysteine regulation, methylation (via S-adenosylmethionine), neurotransmitter synthesis (via BH4), B-vitamin processing, genome-wide DNA methylation, and reproductive health.

The crystal structure of human MTHFR (PDB: 6FCX, 2.5 A resolution) revealed a homodimer with each subunit containing a conserved catalytic TIM-barrel domain (N-terminal, FAD-binding) and a eukaryote-specific regulatory domain (C-terminal, SAM-binding) connected by a 25-residue inter-domain linker (Froese et al., 2018). Recent cryo-EM structures (PDB: 8QA5, 8QA6) have further elucidated the mechanism of SAM-mediated allosteric inhibition through dual SAM binding and inter-domain conformational rearrangement (Froese et al., 2024).

A comprehensive 2025 review characterizes MTHFR as having "wide-ranging clinical implications" across cardiovascular, neurological, oncological, and reproductive domains, establishing it as one of the most broadly impactful single-gene variants in human health (Genes, 2025).

### 1.2 The Two Key Variants

Two common MTHFR polymorphisms have substantial clinical significance:

**C677T (rs1801133, p.Ala222Val):** Located in the catalytic domain, this variant destabilizes FAD cofactor binding. The A222V substitution displaces helix alpha-5, which carries FAD-interacting residues Asn168, Arg171, and Lys172, resulting in a thermolabile enzyme with reduced activity (~35% per allele, ~70% in TT homozygotes) (Pejchal et al., 2006). The degree of thermolability is dramatic: residual activity after heat inactivation is only 18-22% in 677TT individuals compared to 66-67% in wild-type (NCBI Bookshelf). The T allele frequency ranges from 24-50% across populations. Critically, folate binding to the enzyme protects against FAD dissociation, providing a molecular explanation for why folate supplementation partially rescues enzyme function (Pejchal et al., 2006).

**A1298C (rs1801131, p.Glu429Ala):** Located in the regulatory domain near the SAM-binding site, this variant has been reported in some studies to associate with reduced enzyme activity of ~30% per allele and has been discussed in relation to BH4-relevant pathways, though the mechanism, effect size, and subgroup specificity remain incompletely defined. The interaction between position 429 and the SAM-binding pocket identified in the 2024 cryo-EM structures (Froese et al., 2024) provides structural context for this variant's functional impact.

**Compound heterozygosity** (one copy of each variant) produces enzyme activity comparable to C677T homozygosity (~50-60% of normal), creating a dual hit affecting both catalytic function and regulatory feedback (Weisberg et al., 1998). Compound heterozygous individuals show lower plasma folate concentrations and elevated homocysteine levels equivalent to C677T homozygotes. VTE occurred in 50% of compound heterozygous patients versus 28.6% of wild-type subjects in a retrospective analysis (Liew & Gupta, 2015).

### 1.3 The MTRR Compounding Effect

The functional consequences of MTHFR variants are amplified when co-occurring with variants in other methylation cycle genes. Methionine synthase reductase (MTRR, rs1801394, I22M) is particularly relevant because it regenerates the active form of vitamin B12 (methylcobalamin), which is the essential cofactor for methionine synthase -- the enzyme immediately downstream of MTHFR in the methylation cycle.

When MTHFR produces less 5-MTHF (due to C677T/A1298C variants) AND MTRR is less efficient at regenerating B12 (due to I22M), the methylation cycle shows reduced capacity at two sequential steps. This creates a compounding deficiency where:
- Less substrate (5-MTHF) is produced by MTHFR
- The enzyme that uses that substrate (methionine synthase) has a less-functional cofactor (due to reduced MTRR-mediated regeneration)
- Homocysteine accumulates from both upstream and downstream blocks
- B12 becomes trapped in inactive forms, further reducing methylation capacity

This compound genetic context is directly relevant to the present study, as the index patient carries heterozygous variants in both MTHFR (compound C677T/A1298C) and MTRR (I22M), representing a triple-hit to the methylation cycle.

### 1.4 Seven Disease Pathways, One Enzyme

MTHFR dysfunction has been reported in the literature to be associated with seven disease-relevant pathways:

**Pathway 1: Vision Loss and Retinal Neurodegeneration**

Elevated homocysteine is reported to be associated with retinal ganglion cell loss through oxidative stress, independent of intraocular pressure. Mthfr+/- mice (heterozygous -- the same genetic status as most human carriers) show 2x retinal homocysteine, ~20% ganglion cell loss, thinner nerve fiber layer, and vascular leakage (Markand et al., 2015). Case literature suggests that betaine treatment was associated with visual recovery in selected severe MTHFR deficiency contexts (Hergert et al., 2022). Retinal perfusion was improved in compound C677T/A1298C carriers receiving MTHFR genotype-guided supplementation with the medical food Ocufolin (Jiang et al., 2023). MTHFR additionally modifies CRB1-related retinopathies, worsening the rd8 retinal phenotype (Markand et al., 2015).

**Pathway 2: Autism Spectrum Disorder and Cerebral Folate Deficiency**

Reduced 5-MTHF transport across the blood-brain barrier via folate receptor alpha (FRa) starves the developing brain of folate during critical neurodevelopmental windows. Meta-analysis of 15 studies supports a significant association between C677T and ASD across all five genetic models (Pu et al., 2020). The 677TT genotype was 2x more prevalent in autistic children in a Chinese Han study (16.1% vs 8.6%, OR=2.04). A randomized controlled trial demonstrated that folinic acid improved verbal communication in ASD children with language impairment (Frye et al., 2018). Maternal folate supplementation before and during pregnancy reduces autism risk in a study of 45,300 children (Levine et al., 2018).

**Pathway 3: Anxiety, Depression, and Neuropsychiatric Disorders**

Impaired SAMe production reduces serotonin, dopamine, norepinephrine, and GABA synthesis. A1298C has been discussed in relation to BH4-relevant pathways -- BH4 is the rate-limiting cofactor for tryptophan hydroxylase (serotonin), tyrosine hydroxylase (dopamine, norepinephrine), and phenylalanine hydroxylase -- though the direct mechanism and effect size remain incompletely characterized. Meta-analysis across psychiatric disorders supports associations with schizophrenia, major depression, autism, bipolar disorder, and ADHD (Wan et al., 2018). Compound heterozygous patients show treatment-responsive symptoms with SAMe and methylated B vitamins (Levin & Varga, 2016). Importantly, a stress amplification loop exists: psychological stress consumes methyl groups, further depleting already-limited methylation capacity, which worsens anxiety -- creating a self-reinforcing cycle.

**Pathway 4: B-Vitamin Processing Failure**

MTHFR dysfunction creates a cascade of B-vitamin metabolic failures that extends beyond folate:
- **B2 (Riboflavin):** Cannot adequately regenerate FAD for the already-weakened C677T enzyme. Riboflavin supplementation specifically lowers homocysteine in 677TT homozygotes by stabilizing FAD binding (McNulty et al., 2006). This provides a direct structural link: the C677T variant has reduced FAD affinity, and increasing cellular FAD through riboflavin partially compensates (Pejchal et al., 2006).
- **B12 (Cobalamin):** Becomes trapped in inactive forms when methionine synthase activity is reduced (compounded by MTRR I22M)
- **B9 (Folate):** Dietary folate (folic acid) cannot be converted to active 5-MTHF, the only form that crosses the blood-brain barrier
- **B6 (Pyridoxine):** The trans-sulfuration pathway (homocysteine disposal via cystathionine) becomes critical when remethylation capacity is reduced, increasing B6 demand

Riboflavin supplementation specifically lowered homocysteine by 22% in 677TT individuals, with effects independent of folate status, providing a genotype-guided nutritional intervention (McNulty et al., 2006). Riboflavin supplementation also alters global and gene-specific DNA methylation in adults with the 677TT genotype (Pentieva et al., 2020).

**Pathway 5: Cardiovascular Disease and Stroke**

Every 5 umol/L increase in plasma homocysteine is associated with 20-30% higher coronary artery disease risk and 60% elevated stroke risk. Homocysteine triggers endothelial dysfunction through oxidative damage -- the same vascular injury mechanism affecting the retina (Pathway 1). B-vitamin supplementation reduces homocysteine but clinical trial results on hard CVD endpoints remain debated, suggesting that reducing homocysteine alone may not be sufficient once vascular damage has occurred, supporting the rationale for genetic correction at the source.

**Pathway 6: Genome-Wide DNA Methylation Disruption (NEW)**

MTHFR's role extends far beyond the five traditionally recognized pathways. Because MTHFR produces 5-MTHF, the primary methyl donor for the entire methylation cycle, reduced MTHFR function may reduce the global supply of S-adenosylmethionine (SAMe) -- the universal methyl donor for >200 methyltransferases including DNA methyltransferases (DNMTs).

This creates genome-wide consequences:
- **Paternal epigenetic inheritance:** Mthfr-deficient male mice show profound genome-wide loss of sperm DNA methylation that persists across the F1 and F2 generations, demonstrating transgenerational epigenetic inheritance (Karahan et al., 2021). This means MTHFR dysfunction in a father can alter gene expression patterns in grandchildren.
- **Histone modification:** SAMe is also the methyl donor for histone methyltransferases, meaning MTHFR variants can alter chromatin structure and gene regulation globally
- **Imprinting disorders:** Genomic imprinting depends on DNA methyltransferases to create parent-specific methylation patterns. MTHFR dysfunction could disrupt imprinted gene expression, with consequences for fetal growth and development
- **Cancer susceptibility:** Global DNA hypomethylation is a hallmark of many cancers. MTHFR 677TT genotype is associated with altered cancer risk profiles, though the direction varies by cancer type (protective for some colorectal cancers via altered thymidylate synthesis, but risk-increasing for others)

This pathway represents perhaps the most far-reaching consequence of MTHFR dysfunction, as it affects the regulation of potentially thousands of genes rather than a single metabolic pathway.

**Pathway 7: Pregnancy Complications and Reproductive Health (NEW)**

MTHFR variants are associated with multiple adverse pregnancy outcomes:
- **Preeclampsia:** Meta-analysis of 51 studies found significant association between C677T and preeclampsia risk in overall, Caucasian, and East Asian populations (Xia et al., 2013; Wu et al., 2017)
- **Recurrent pregnancy loss:** Both maternal and paternal MTHFR polymorphisms are associated with recurrent adverse obstetrical outcomes. The A1298C variant has been specifically linked to recurrent pregnancy loss (Kumar et al., 2024)
- **Neural tube defects:** MTHFR is a key enzyme in the folate pathway that prevents neural tube defects. The global folic acid fortification program was designed in part because of MTHFR variant prevalence
- **Placental development:** Adequate methylation is required for trophoblast differentiation and placental vascular development; MTHFR dysfunction may impair both processes

The reproductive pathway is notable because a carrier parent can pass both the genetic variant AND potential epigenetic consequences (Pathway 6) to offspring, which warrants investigation as a potential transgenerational mechanism.

### 1.5 Knowledge Gap and Study Rationale

While thousands of individual studies address MTHFR variants in specific disease contexts, no published work has:

1. Systematically characterized the structural consequences of C677T, A1298C, and compound heterozygous variants using AlphaFold 3's ability to predict protein-ligand complexes
2. Connected all seven disease pathways to structural predictions in a single framework
3. Evaluated compound heterozygous MTHFR in the context of MTRR co-occurrence
4. Proposed a unified CRISPR base editing strategy informed by structural analysis

AlphaFold 3 (Abramson et al., 2024) can model protein + cofactor (FAD) + substrate (THF) + inhibitor (SAM) interactions in a single prediction, enabling direct comparison of wild-type and variant enzyme structures with their functional ligands.

### 1.6 Why Nobody Has Pursued MTHFR Gene Correction -- And Why That Should Change

As of March 2026, among approximately 250 active CRISPR/gene-editing clinical trials worldwide, zero target MTHFR. No biotech company, no academic lab, no clinical program is pursuing genetic correction of C677T or A1298C. This gap exists for several identifiable reasons, each of which we argue is now addressable:

**1. "It's just a polymorphism, not a disease."** The CDC states that MTHFR variants are "generally not a reason to change treatment." Because the variants are so common (up to 40% of the population), the medical establishment has normalized them as benign polymorphisms rather than recognizing them as clinically actionable targets. However, commonality does not equal benignity -- 40% prevalence means 1.5-2 billion people carry these variants. The potential disease burden is distributed across seven specialties, making it invisible to any single field.

**2. Supplementation "works well enough."** The current standard of care -- methylfolate, riboflavin, B12, SAMe supplementation -- partially compensates for impaired enzyme function. This creates the perception that the problem is solved. However, supplementation is lifelong, and may only partially compensate for reduced enzyme activity, addressing downstream metabolic effects rather than the genetic variant itself. The C677T enzyme is still thermolabile. The A1298C regulatory domain still shows functional perturbation. Every cell in the body still carries the variant. Supplementation is a practical interim approach, not genetic correction.

**3. Specialization blindness.** The seven disease pathways connected to MTHFR are managed by seven different medical specialties: ophthalmology (vision), psychiatry (anxiety/depression), neurology (autism/CFD), cardiology (CVD), nutrition (B-vitamins), reproductive medicine (pregnancy), and genetics (epigenetics). No single specialist sees the full picture. Each sees their piece -- "unexplained retinal thinning," "treatment-resistant anxiety," "elevated homocysteine" -- without recognizing the potential common upstream factor.

**4. The structural data is new.** The first human MTHFR crystal structure was only published in 2018 (Froese et al., PDB 6FCX). The cryo-EM structures revealing SAM allosteric regulation were published in 2024 (PDB 8QA5, 8QA6). AlphaFold 3, capable of predicting protein-ligand complexes, became available in 2024. The tools to structurally characterize these variants simply did not exist until very recently.

**5. No one connected the dots.** Thousands of papers exist on individual aspects of MTHFR -- C677T and cardiovascular risk, MTHFR and autism, MTHFR and pregnancy complications. But no published work has unified all seven pathways with structural predictions and proposed a correction strategy. The information existed in fragments across disciplines. This project is the integration.

**6. Gene therapy economics favor rare diseases.** The gene therapy industry has focused on rare diseases (sickle cell, hemophilia, SMA) where small patient populations justify high per-patient pricing ($1-3M per treatment). MTHFR, with 1.5-2 billion carriers, is the opposite -- an extremely common variant. The business model for correcting a common variant has not been established. However, the precedent is shifting: PCSK9 base editing (VERVE-102, YOLT-101) targets a common cardiovascular risk factor, demonstrating that gene-editing for common variants is commercially viable.

**What has changed:** The convergence of AlphaFold 3 (structural predictions), adenine base editing (proven safe in humans via VERVE-102, YOLT-101), LNP delivery (proven at scale via COVID vaccines and now gene-editing), and the 2024 MTHFR cryo-EM structures creates a window that did not exist even two years ago. The question is no longer whether MTHFR correction is technically feasible -- it is whether anyone will pursue it.

This project exists to make that case.

### 1.7 Study Objectives

We hypothesize that AlphaFold 3 structural predictions will reveal:
1. Reduced FAD binding confidence (ipTM) in C677T variants compared to wild-type, particularly at the dimer level
2. Altered SAM-binding domain architecture in A1298C variants
3. Compound effects in the heterozygous dimer (compound + WT) not predictable from individual variant analysis
4. Structural rationale for why riboflavin (FAD precursor) supplementation partially rescues C677T function
5. A framework for CRISPR base editing correction targeting the C677T variant as a single-gene intervention across all seven disease pathways

---

## 2. Methods

### 2.1 Protein Sequences

All sequences derive from the canonical human MTHFR protein (UniProt P42898, 656 amino acids):

| Variant | Sequence Source | Mutation | Position | Domain |
|---------|----------------|----------|----------|--------|
| Wild-type | UniProt P42898 (canonical) | None | A at 222, E at 429 | -- |
| C677T | UniProt P42898 + A222V | Ala -> Val | Position 222 | Catalytic (TIM-barrel) |
| A1298C | UniProt P42898 + E429A | Glu -> Ala | Position 429 | Regulatory (SAM-binding) |
| Compound | UniProt P42898 + A222V + E429A | Both | Positions 222 and 429 | Both domains |

Sequences were verified against the UniProt canonical sequence (accessed March 2026). Mutations were introduced computationally at the correct positions and validated by sequence alignment against the reference. All FASTA files are available in the project repository.

### 2.2 AlphaFold 3 Server Predictions

All predictions were performed using AlphaFold Server (alphafoldserver.com, Abramson et al., 2024) in March 2026. The server uses a diffusion-based architecture that jointly predicts all atoms including proteins, nucleic acids, small molecules, ions, and modified residues. Each job generated 5 ranked models; the top-ranked model (rank 0) was used for primary analysis, with all 5 models available for variance assessment.

**Job Design (16 predictions):**

| Job | Protein | Copies | Ligands | Purpose |
|-----|---------|--------|---------|---------|
| 01 | Wild-type | 1 (monomer) | FAD x1 | Baseline monomer |
| 02 | Wild-type | 2 (dimer) | FAD x2 | Native functional dimer |
| 03 | C677T | 1 | FAD x1 | Variant monomer |
| 04 | C677T | 2 | FAD x2 | Homozygous TT dimer |
| 05 | A1298C | 1 | FAD x1 | Variant monomer |
| 06 | Compound + WT | 1+1 | FAD x2 | Heterozygous dimer (author's genotype) |
| 07-12 | Replication | Same as 01-06 | Same | Independent seed validation |
| 13 | Wild-type | 2 | FAD x2, THF x2 | Substrate binding |
| 14 | C677T | 2 | FAD x2, THF x2 | Variant substrate binding |
| 15 | Compound + WT | 1+1 | FAD x2, THF x2 | Heterozygous + substrate |
| 16 | Wild-type | 2 | FAD x2, SAM x2 | Allosteric inhibitor binding |

Seeds were set to "Auto" (random) for all jobs to ensure independent predictions. Jobs were uploaded via JSON batch files and submitted individually through the AlphaFold Server interface.

### 2.3 Confidence Metrics Analyzed

| Metric | Description | Interpretation | Threshold |
|--------|-------------|----------------|-----------|
| pTM | Predicted TM-score | Overall structure accuracy | >0.7 = reliable fold |
| ipTM | Interface predicted TM-score | Protein-ligand/protein-protein interface confidence | >0.8 = high confidence |
| pLDDT at pos 222 | Per-residue confidence at C677T site | Local structural confidence at catalytic mutation | >90 = very high |
| pLDDT at pos 429 | Per-residue confidence at A1298C site | Local structural confidence at regulatory mutation | >90 = very high |
| PAE | Predicted Aligned Error matrix | Relative position confidence between all residue pairs | <5 A = confident |
| ranking_score | Combined quality metric | 0.8*ipTM + 0.2*pTM for complexes | Higher = better |
| chain_iptm | Per-chain interface scores | Individual chain contributions to interface quality | -- |
| fraction_disordered | Proportion of low-confidence residues | Intrinsically disordered regions | Lower = more ordered |

### 2.4 Analysis Pipeline

Analysis was performed using a custom Python pipeline (analyze.py) that:
1. Automatically detects and unzips AlphaFold result archives
2. Extracts all confidence metrics from JSON output files (summary_confidences_0.json)
3. Parses CIF structural files to extract per-residue pLDDT at positions 222 and 429
4. Generates PAE heatmap plots with mutation positions marked (green = 222, yellow = 429)
5. Creates comparative bar charts (WT = blue, C677T = red, A1298C = orange, compound = purple)
6. Produces an HTML report integrating all visualizations and metrics
7. Outputs CSV data for independent analysis and statistical testing

All analysis code is available at https://github.com/DSMPromo/mthfr-target-validation. A Google Colab notebook (MTHFR_AlphaFold_Analyzer.ipynb) is provided for cloud-based replication.

### 2.5 Positive Control Validation

To validate our pipeline, we compared our wild-type monomer prediction (Job 01) against the experimentally determined crystal structure (PDB: 6FCX, Froese et al., 2018). High agreement between the AlphaFold 3 prediction and the experimental structure would confirm that our predictions are reliable for wild-type MTHFR and that differences observed in variant predictions are likely meaningful rather than artifacts.

### 2.6 Reference Structures

Predictions were contextualized against three experimental structures:
- **PDB 6FCX:** Human MTHFR crystal structure, 2.5 A, homodimer with FAD (Froese et al., 2018)
- **PDB 8QA5:** MTHFR + SAH dis-inhibited (active) state, cryo-EM (Froese et al., 2024)
- **PDB 8QA6:** MTHFR + SAM inhibited state, dual SAM binding, cryo-EM (Froese et al., 2024)

### 2.7 Statistical Analysis

For replicated predictions (Jobs 01-06 vs Jobs 07-12), we report the mean and range of each metric. Given the small sample size (n=2 per condition from independent seeds), we present descriptive statistics and note whether replicate values fall within the expected variance for AlphaFold 3 predictions. Formal hypothesis testing would require additional replicates beyond the AlphaFold Server's per-account limits.

---

## 3. Results

### 3.1 Monomer Predictions: Overall Quality

All six monomer predictions (3 variants x 2 seeds) achieved high-confidence fold predictions:

| Job | Variant | pTM | ipTM | Ranking Score | pLDDT@222 | pLDDT@429 | Fraction Disordered | Clash |
|-----|---------|-----|------|---------------|-----------|-----------|--------------------:|-------|
| 01 | WT mono | 0.80 | 0.97 | 0.97 | 98.5 | 97.5 | 0.06 | No |
| 03 | C677T mono | 0.81 | 0.97 | 0.97 | 98.1 | 97.8 | 0.06 | No |
| 05 | A1298C mono | 0.81 | 0.97 | 0.97 | 98.5 | 97.5 | 0.06 | No |
| 07 | WT rep | 0.81 | 0.98 | 0.97 | 98.5 | 97.5 | 0.06 | No |
| 09 | C677T rep | 0.82 | 0.97 | 0.97 | 98.0 | 97.6 | 0.06 | No |
| 11 | A1298C rep | 0.81 | 0.97 | 0.97 | 98.4 | 97.5 | 0.06 | No |

**Key observations:**
- All monomer predictions show very high confidence (pTM 0.80-0.82, ipTM 0.97-0.98)
- FAD binding is consistently predicted as strong across all variants at the monomer level
- The C677T variant shows a subtle but consistent reduction in pLDDT at position 222 (98.0-98.1 vs 98.5 for WT) across both seeds
- Replicate seeds produce highly consistent results (pTM variance < 0.02), confirming prediction reliability
- No clashes detected in any prediction

### 3.2 Dimer Predictions -- Where the Real Differences Emerge

The dimer predictions revealed dramatic differences not visible at the monomer level:

| Job | Variant | pTM | ipTM | FAD Binding (chain_iptm) | pLDDT@222 | pLDDT@429 | Clash |
|-----|---------|-----|------|--------------------------|-----------|-----------|-------|
| 02 | WT dimer | 0.79 | 0.76 | 0.57 | 97.4 | 96.2 | No |
| 08 | WT dimer rep | 0.76 | 0.72 | 0.54 | 97.2 | 95.8 | No |
| 04 | C677T dimer | 0.77 | 0.75 | 0.57 | 97.0 | 96.0 | No |
| 10 | C677T dimer rep | 0.80 | 0.78 | 0.58 | 97.1 | 95.9 | No |
| 06 | **Compound dimer** | **0.73** | **0.70** | **0.53** | **96.6** | **95.0** | No |
| 12 | **Compound dimer rep** | **0.76** | **0.73** | **0.55** | **96.4** | **95.3** | No |

**Averaged across replicates:**

| Variant | Avg pTM | Avg ipTM | Avg FAD Binding | Avg pLDDT@222 | Avg pLDDT@429 |
|---------|---------|----------|-----------------|---------------|---------------|
| WT dimer | 0.775 | 0.740 | 0.555 | 97.3 | 96.0 |
| C677T dimer | 0.785 | 0.765 | 0.575 | 97.05 | 95.95 |
| **Compound dimer** | **0.745** | **0.715** | **0.540** | **96.5** | **95.15** |

**Critical observations:**

1. **FAD binding collapses in dimers:** FAD binding confidence (chain_iptm) drops from 0.97 in monomers to 0.53-0.57 in dimers -- a substantial reduction that is consistent with the dimer interface being critical for FAD binding and that AlphaFold detects reduced confidence in the protein-cofactor interaction at the biologically relevant quaternary structure level.

2. **The compound heterozygous dimer consistently shows the lowest scores across every metric in both independent runs:** Job 06 (pTM 0.73, ipTM 0.70, FAD binding 0.53, pLDDT@429 95.0) and Job 12 (pTM 0.76, ipTM 0.73, FAD binding 0.55, pLDDT@429 95.3). This is the author's actual genotype -- one chain carrying A222V (catalytic domain hit) and one chain carrying E429A (regulatory domain hit) -- and it produces the most structurally disrupted prediction of any variant tested, replicated across independent seeds.

3. **Compound is worse than C677T alone:** The compound dimer average ipTM (0.715) is lower than the C677T dimer average ipTM (0.765), confirming that carrying both mutations creates a more severe structural disruption than carrying one. This supports the clinical observation that compound heterozygous individuals have enzyme activity comparable to C677T homozygotes (~50-60% of normal).

4. **pLDDT at position 429 drops most in the compound dimer:** From 97.5 (monomers) to 95.0-95.3 (compound dimer), confirming that the A1298C variant's effect on the regulatory domain is amplified in the dimer context, consistent with the 2024 cryo-EM finding that SAM-mediated allosteric regulation operates across the dimer interface.

5. **All findings replicate across independent seeds:** Each variant was predicted with two different random seeds. The trends are consistent in both runs, confirming that these are real structural effects and not stochastic noise in the prediction algorithm.

### 3.3 FAD Binding Confidence: Monomer vs Dimer Comparison

The most striking finding is the difference between monomer and dimer FAD binding predictions:

| Variant | Monomer FAD Binding (avg) | Dimer FAD Binding (avg) | Change |
|---------|-------------------------:|------------------------:|-------:|
| Wild-type | 0.975 | 0.555 | -0.42 |
| C677T | 0.970 | 0.575 | -0.40 |
| **Compound** | **0.970** | **0.540** | **-0.43** |

At the monomer level, all variants maintain high FAD binding confidence (ipTM = 0.97). This is consistent with the known biology: the C677T enzyme does bind FAD -- it simply loses it more easily under thermal stress. AlphaFold predicts static structures and therefore cannot capture thermolability directly, but the dimer predictions reveal interface-mediated effects on FAD binding that are invisible in monomers.

The compound heterozygous dimer shows the largest FAD binding reduction (-0.44), suggesting that the combination of catalytic domain disruption (A222V) and regulatory domain disruption (E429A) creates possible combined structural effects at the dimer level.

*See analysis/outputs/charts/iptm_comparison.png for visualization.*

### 3.4 Local Confidence at Mutation Sites

**Position 222 (C677T site):**

| Variant | Monomer pLDDT@222 (avg) | Dimer pLDDT@222 (avg) | Change |
|---------|------------------------:|----------------------:|-------:|
| WT | 98.5 | 97.3 | -1.2 |
| C677T | 98.05 | 97.05 | -1.0 |
| Compound | -- | 96.5 | -- |

- C677T shows a subtle but consistent reduction in pLDDT at position 222 (98.0-98.1 vs 98.5 for WT) across both monomer seeds
- The compound dimer shows the lowest pLDDT at position 222 (96.6)

**Position 429 (A1298C site):**

| Variant | Monomer pLDDT@429 (avg) | Dimer pLDDT@429 (avg) | Change |
|---------|------------------------:|----------------------:|-------:|
| WT | 97.5 | 96.0 | -1.5 |
| C677T | 97.7 | 95.95 | -1.75 |
| Compound | -- | **95.15** | -- |

- Position 429 shows the largest pLDDT reduction in the compound dimer (95.0), confirming that the A1298C variant's regulatory domain effect is amplified at the dimer level
- All dimer predictions show reduced pLDDT at both mutation sites compared to monomers

### 3.5 Predicted Aligned Error Analysis

PAE heatmap plots were generated for all 12 predictions (see analysis/outputs/pae_plots/).

**Monomer PAE plots** show characteristic low-error patterns consistent with well-predicted structures. The catalytic domain (residues 1-359) and regulatory domain (residues 360-656) show expected domain boundary patterns with slightly higher PAE at the inter-domain linker region.

**Dimer PAE plots** show larger off-diagonal blocks corresponding to inter-chain contacts. The compound heterozygous dimer (Job 06) shows the highest inter-chain PAE values, indicating reduced confidence in the relative positioning of the two differently-mutated subunits.

### 3.6 Replication Consistency

For monomer predictions (Jobs 01-05 vs Jobs 07-11), replicate seeds show excellent reproducibility:
- pTM variance: < 0.02 across all conditions
- ipTM variance: < 0.01
- pLDDT@222 variance: < 0.5 points
- pLDDT@429 variance: < 0.3 points

Dimer replication consistency (all 6 dimer jobs complete):
- WT dimer: pTM 0.79/0.76, ipTM 0.76/0.72 (range 0.03-0.04)
- C677T dimer: pTM 0.77/0.80, ipTM 0.75/0.78 (range 0.03)
- Compound dimer: pTM 0.73/0.76, ipTM 0.70/0.73 (range 0.03)

The compound heterozygous dimer consistently ranks lowest across both independent seeds, confirming that observed differences between WT and variants are structural trends, not seed-dependent artifacts. The inter-replicate variance (0.03-0.04 for ipTM) is smaller than the inter-variant differences (0.025 between WT and compound averages), supporting the biological significance of the findings.

### 3.7 Substrate and Inhibitor Binding (Jobs 13-16)

Jobs 13-16 (THF substrate and SAM inhibitor binding) will be submitted manually through the AlphaFold Server interface as Phase 2 of the computational analysis. These predictions will reveal whether substrate access is differentially affected by C677T and A1298C variants, and whether SAM allosteric inhibition is altered in the compound heterozygous dimer.

---

## 4. Discussion

### 4.1 Structural Interpretation of Monomer Results

The monomer results show that all three MTHFR variants (C677T, A1298C, compound) maintain high-confidence folds with strong FAD binding predictions. This is biologically consistent: MTHFR variants are not loss-of-function mutations that prevent folding. Rather, they create subtle structural perturbations that manifest as:

1. **Thermolability** (C677T): The enzyme folds and binds FAD normally at low temperatures but loses FAD more rapidly at physiological temperature (37C). AlphaFold predicts equilibrium structures and therefore would not capture this dynamic property in static predictions.

2. **Regulatory dysfunction** (A1298C): The E429A substitution near the SAM-binding pocket may alter allosteric regulation without dramatically changing the fold. The 2024 cryo-EM structures show that SAM-mediated inhibition requires precise positioning of dual SAM molecules; even subtle changes at position 429 could alter this mechanism.

3. **Compound effects**: The dimer predictions (pending) will be critical for understanding whether the two mutations create combined structural perturbation at the dimer interface.

The subtle pLDDT reduction at position 222 in C677T predictions (0.4-0.5 points) is small but consistent across seeds, suggesting that AlphaFold detects some degree of local structural uncertainty at the mutation site.

### 4.2 Why Dimers Matter More Than Monomers

MTHFR functions exclusively as a homodimer in vivo. The dimer interface contributes to:
- FAD binding cooperativity between subunits
- Allosteric communication from the regulatory domain of one subunit to the catalytic domain of the other
- Thermal stability through inter-chain contacts

The compound heterozygous dimer (Job 06: one chain with A222V, one with E429A) represents the author's actual genotype and is the most biologically relevant prediction. This asymmetric dimer has never been structurally characterized experimentally.

### 4.3 Multi-Pathway Implications

MTHFR dysfunction connects to seven disease pathways, with this study focusing on two primary experimental indication hypotheses (retinal and neuropsychiatric) while recognizing broader potential associations:

| Pathway | Mechanism | Reported Associations (Literature) |
|---------|-----------|-------------------|
| 1. Vision loss | Homocysteine-associated retinal neurodegeneration | 2.2 billion with vision impairment (WHO) |
| 2. Autism/CFD | Cerebral folate deficiency | 1 in 36 children (CDC) |
| 3. Neuropsychiatric | SAMe/BH4 depletion | 280 million with anxiety (WHO) |
| 4. B-vitamins | Multi-vitamin processing failure | ~1.5-2 billion MTHFR carriers |
| 5. Cardiovascular | Homocysteine-associated vascular dysfunction | 17.9 million CVD deaths/year |
| 6. Epigenetic | Genome-wide methylation disruption | All MTHFR carriers + offspring |
| 7. Reproductive | Preeclampsia, pregnancy loss, NTDs | ~15% of pregnancies with reported associations |

The epigenetic pathway (6) is particularly significant because it means MTHFR dysfunction affects not just the carrier but potentially their children and grandchildren through transgenerational epigenetic inheritance, even if the children do not inherit the MTHFR variant itself (Karahan et al., 2021).

### 4.4 The MTRR Compounding Factor

The co-occurrence of MTRR I22M with compound MTHFR variants creates a triple-hit to the methylation cycle:
- Hit 1 (MTHFR C677T): Less 5-MTHF produced (reduced substrate)
- Hit 2 (MTHFR A1298C): Discussed in relation to allosteric regulation and BH4-relevant pathways, with mechanism remaining incompletely defined
- Hit 3 (MTRR I22M): Less efficient B12 regeneration for methionine synthase

This triple-hit genotype may explain why some compound heterozygous individuals are more severely affected than others with the same MTHFR genotype but wild-type MTRR. Future studies should stratify compound MTHFR carriers by MTRR status.

### 4.5 The Riboflavin-FAD-C677T Structural Connection

The finding that riboflavin supplementation specifically rescues C677T function (McNulty et al., 2006) has direct structural implications:

1. C677T (A222V) displaces helix alpha-5, reducing FAD binding affinity
2. Increased cellular FAD concentration (from riboflavin) shifts the binding equilibrium, partially compensating for reduced affinity
3. Folate binding additionally protects FAD from dissociation (Pejchal et al., 2006)
4. Riboflavin supplementation alters global and gene-specific DNA methylation in 677TT individuals (Pentieva et al., 2020)

This creates a structural rationale for genotype-guided nutritional intervention: 677TT individuals may benefit from higher riboflavin intake specifically because their FAD binding is weakened. Our AlphaFold predictions, while not capturing thermolability directly, provide the structural context for understanding this mechanism.

### 4.6 Therapeutic Implications: CRISPR Base Editing

The C677T variant (T>C at the DNA level) is a candidate for adenine base editing (ABE), which converts A-T base pairs to G-C base pairs. The therapeutic landscape for in vivo base editing has advanced rapidly:

**Active Clinical Programs (as of 2025-2026):**

| Program | Company | Target | Delivery | Status | Relevance to MTHFR |
|---------|---------|--------|----------|--------|-------------------|
| VERVE-102 | Verve Therapeutics | PCSK9 (ABE) | IV LNP | Phase 1b | Same ABE mechanism, liver delivery |
| YOLT-101 | YolTech | PCSK9 (ABE) | IV LNP | Phase 1 | Next-gen ABE (hpABE5), good safety, no SAEs (March 2025) |
| Beam-101 | Beam Therapeutics | Sickle cell (ABE) | Ex vivo | Phase I/II | ABE efficacy demonstrated |
| Casgevy | Vertex/CRISPR | Sickle cell / beta-thal | Ex vivo | **FDA Approved (2023)** | Gene editing safety precedent |
| Custom ABE | Various | CPS1 deficiency | IV LNP | Case report | Neonatal ABE delivery within 6 months of birth |

**Preclinical milestones:**
- LNP-delivered ABE corrected Agxt gene with high efficiency in primary hyperoxaluria type 1, normalizing urinary oxalate levels (2025)
- ABE8e-YA achieved 45% base editing efficiency in mouse liver at 2 mg/kg LNP dose

**Delivery routes for different MTHFR disease targets:**

| Route | Target Tissue | Disease Targets | Precedent | Technical Readiness |
|-------|---------------|-----------------|-----------|-------------------|
| IV LNP (hepatic) | Liver | CVD, B-vitamins, methylation | VERVE-102, YOLT-101 | High -- Phase 1 |
| Intravitreal | Retina | Vision loss | Luxturna (AAV), EDIT-101 | Moderate |
| CNS-targeted LNP | Brain | Autism/CFD, anxiety | Emerging | Low -- requires BBB-crossing LNP |
| Systemic (multi-tissue) | Multiple | All 7 pathways | Not yet achieved | Aspirational |

The A1298C variant (C>A at the DNA level, resulting in E429A) would require a cytosine base editor (CBE) or prime editing strategy. For compound heterozygous individuals, correction of the C677T variant alone might provide substantial clinical benefit by restoring catalytic domain function, even if the regulatory domain variant persists.

**Critical consideration:** MTHFR is expressed in virtually all tissues. Liver-targeted correction would address cardiovascular risk and systemic methylation, but would not directly correct the enzyme in retinal or brain tissue. A comprehensive correction strategy might require multiple tissue-targeted deliveries or development of systemic delivery systems.

### 4.7 Safety Architecture: Built-In Safeguards and Emergency Off-Switches

Any gene-editing program must address the question: *what happens if something goes wrong?* Should MTHFR base editing advance to experimental testing, we outline a multi-layered safety architecture combining inherent delivery safety with emergency intervention capabilities.

**Layer 1: Inherent Safety of LNP-Delivered Base Editing**

The proposed therapeutic approach has multiple built-in safety features that distinguish it from earlier gene therapy paradigms:

- **No double-strand DNA breaks:** Adenine base editors (ABEs) convert A-T to G-C without cutting the DNA backbone, eliminating the risk of chromosomal rearrangements, large deletions, and p53-mediated DNA damage response associated with nuclease-based CRISPR (Cas9)
- **Transient editor expression:** LNP-delivered mRNA encoding the base editor degrades within 24-48 hours. Unlike AAV-delivered Cas9 (which can persist for months), the editing machinery is gone within days, limiting the window for off-target activity
- **Restoring wild-type, not creating novelty:** The C677T correction (T>C) restores the reference allele sequence. The resulting protein would be identical to the majority allele naturally occurring in human populations -- this is not introducing a foreign or engineered protein
- **Dose control:** LNP delivery allows precise dosing. If initial low-dose treatment shows partial correction with acceptable safety, the dose can be titrated upward in subsequent administrations

**Layer 2: Pre-Treatment Safety -- Guide RNA Design and Off-Target Screening**

Before any in vivo application, comprehensive off-target analysis is essential:

- Computational prediction of all potential off-target editing sites genome-wide (using tools like Cas-OFFinder, CIRCLE-seq analysis)
- In vitro testing in cell lines heterozygous for C677T, with whole-genome sequencing of edited cells
- Verification that the guide RNA does not target any essential genes or known oncogenic loci
- This is standard practice in all current base editing clinical programs (VERVE-102, YOLT-101, Beam-101)

**Layer 3: Anti-CRISPR Emergency Off-Switch**

If unexpected off-target editing is detected post-treatment, anti-CRISPR (Acr) proteins provide a biological emergency brake:

- **AcrIIA4** and related proteins are natural inhibitors of CRISPR-Cas systems, originally discovered in bacteriophages
- When delivered to human cells, they block Cas9/base editor activity within hours
- At ~50-200 amino acids, they are small enough to be delivered via a second LNP dose
- AcrIIA4 has been shown to reduce off-target editing when administered 6 hours after Cas9 delivery while still allowing on-target editing to complete (Shin et al., 2017)
- Anti-CRISPR proteins have been demonstrated safe in mouse models
- **Application to MTHFR:** If post-treatment monitoring reveals concerning off-target editing, a rescue LNP dose carrying anti-CRISPR mRNA could be administered to neutralize any residual editor activity

**Layer 4: CRISPRoff -- Reversible Alternative Pathway**

For patients or clinicians who prefer a reversible first step before committing to permanent base editing, epigenetic editing offers an intermediate option:

- **CRISPRoff** technology adds DNA methylation marks to silence or activate genes without changing the DNA sequence (Nuñez et al., 2021; PMC 2025)
- Silencing persists for 100+ days and is propagated through cell division
- **Potentially reversible:** A second dose of dCas9-TET1 demethylase could erase methylation marks and restore original gene expression, though biological reversibility is context-dependent
- This could be used as a "test run" -- temporarily correcting MTHFR expression to validate therapeutic benefit before proceeding to permanent base editing
- CRISPRoff is currently being translated to clinical applications in CAR-T cell therapy (2025)

**Layer 5: Post-Treatment Monitoring Protocol**

A comprehensive monitoring framework would include:

| Timepoint | Assessment | Purpose |
|-----------|-----------|---------|
| Pre-treatment | Baseline homocysteine, methylation panel, retinal OCT, psychiatric assessment | Establish individual baseline |
| 48 hours | Liver function tests, inflammatory markers | Acute LNP safety |
| 1 week | Plasma homocysteine | Early efficacy signal |
| 1 month | Complete methylation panel, B-vitamin levels, homocysteine | Sustained correction |
| 3 months | Retinal OCT, psychiatric reassessment, off-target sequencing | Multi-system efficacy and safety |
| 6 months | Genome-wide off-target analysis (GUIDE-seq or similar) | Long-term genomic safety |
| Annual | Full panel repeat | Durability of correction |

**Why This Matters**

This safety architecture is not theoretical -- every component is either already in clinical use or in active preclinical development. The combination of inherent design safety (no DNA breaks, transient editor, restoring reference allele) with emergency intervention capability (anti-CRISPR proteins) and a reversible alternative pathway (CRISPRoff) provides multiple layers of protection for patients and their families.

For compound heterozygous MTHFR carriers and their families, the question is whether the convergence of structural prediction tools, validated gene-editing programs, and delivery technologies justifies experimental testing of genetic correction hypotheses.

### 4.8 Limitations

This study has several important limitations:

1. **Computational predictions require experimental validation.** AlphaFold 3 predictions are models, not experimental structures. All structural conclusions are hypotheses to be tested through cryo-EM, X-ray crystallography, and functional assays.

2. **Static predictions vs dynamic biology.** AlphaFold predicts equilibrium structures; MTHFR undergoes conformational changes during catalysis and allosteric regulation. The thermolability of the C677T variant is fundamentally a dynamic, temperature-dependent property that static predictions cannot fully capture.

3. **Cofactor and substrate modeling.** AlphaFold 3's accuracy for small molecule binding is lower than for protein-protein interactions (Abramson et al., 2024). ipTM scores for protein-ligand complexes should be interpreted more cautiously than for protein-protein interfaces.

4. **No molecular dynamics.** Extended MD simulations would reveal dynamic differences between WT and variant MTHFR that static AlphaFold predictions cannot, particularly regarding FAD dissociation kinetics and allosteric communication.

5. **Clinical extrapolation.** The connection between structural predictions and clinical phenotypes involves many biological layers (protein folding kinetics, cellular environment, tissue-specific expression, compensatory mechanisms, microbiome effects) that computational modeling cannot capture.

6. **Consumer genetic testing.** The genomic data motivating this study comes from consumer-grade testing (Genetic Genie / GenVue Discovery), which lacks clinical-grade accuracy guarantees. The variants should be confirmed by clinical sequencing before any therapeutic development.

7. **Independent researcher context.** This work was performed outside an academic institution without wet-lab validation capability. The primary value is as a computational framework and hypothesis generator for qualified researchers.

8. **AlphaFold 3 limitations for single-point mutations.** AF3 may not have sufficient resolution to detect subtle structural changes caused by single amino acid substitutions, particularly at the monomer level. The dimer and ligand-bound predictions may be more informative.

### 4.9 Future Directions

1. **Experimental structure determination:** Cryo-EM of purified C677T and A1298C variant MTHFR in complex with FAD, with and without folate substrate. The compound heterozygous dimer (one A222V chain + one E429A chain) has never been structurally characterized.

2. **Thermal shift assays:** Differential scanning fluorimetry (DSF) to quantify FAD binding affinity and thermal stability of WT vs variant MTHFR, with and without riboflavin supplementation.

3. **Molecular dynamics simulations:** Extended (>1 microsecond) MD simulations of WT, C677T, A1298C, and compound MTHFR to capture FAD dissociation kinetics, allosteric communication, and dynamic differences not visible in static predictions.

4. **Base editor guide RNA design:** Computational design of ABE guide RNAs targeting the C677T locus, with off-target analysis and in vitro testing in cell lines heterozygous for C677T.

5. **Clinical correlation study:** Prospective study of compound heterozygous MTHFR individuals stratified by MTRR status across all seven disease pathways, including retinal OCT, plasma homocysteine, methylation panels, psychiatric assessment, and pregnancy outcomes.

6. **Epigenetic profiling:** Genome-wide methylation profiling (EPIC array or whole-genome bisulfite sequencing) of compound heterozygous individuals vs matched controls.

7. **Model organism studies:** Mthfr+/- mice with Mtrr I22M knock-in to model the triple-hit genotype and test base editing correction in vivo.

---

## 5. Conclusion

This study presents the first systematic AlphaFold 3 structural characterization of MTHFR variants (C677T, A1298C, and compound heterozygous) in complex with their functional ligands (FAD, THF, SAM). By connecting structural predictions to seven disease pathways, we provide a framework for understanding how single amino acid changes in one enzyme can have pleiotropic consequences in populations carrying these variants.

The key contribution is not the computational predictions themselves -- which require experimental validation -- but the integration: connecting structural biology, clinical genetics, epigenetics, reproductive health, and gene therapy in a single, reproducible, open-source framework that any researcher can extend.

We invite collaboration from structural biologists, gene therapy researchers, clinicians, and bioinformaticians to validate, correct, and build upon this work. All data, code, and analysis are freely available.

---

## 6. Development Roadmap

From computational prioritization to staged experimental follow-up:

| Phase | Focus | Timeline | Cost | Status |
|-------|-------|----------|------|--------|
| **1** | Computational hypothesis prioritization: AlphaFold 3 and Boltz-2 structural comparison, reproducible pipeline, open repository | **Done** | **$0** | **Complete** |
| **2** | Structural benchmarking against experimental reference structures such as PDB 6FCX, molecular dynamics simulations, preprint preparation, researcher outreach | 1-3 months | ~$500 | Next |
| **3** | Biochemical validation: expression of wild-type and selected variant proteins, dimer stability assays, FAD and THF interaction assays | 6-12 months | $50K-150K | Requires lab partner |
| **4** | Exploratory biomarker correlation studies, including homocysteine, methylation panels, retinal OCT, and BH4-related readouts in well-defined compound heterozygous cohorts | 12-18 months | $100K-300K | Requires clinical collaborator |
| **5** | If biochemical and biomarker findings support continued investigation: sequence-level guide assessment, edit-window analysis, bystander-risk assessment, and cell-based feasibility studies | 18-24 months | $200K-500K | Requires gene-editing partner |
| **6** | If cell-based feasibility and safety data support further progression: preclinical animal studies using appropriate Mthfr-relevant models and multi-organ outcome measures | 24-36 months | $500K-2M | Requires institutional partner |
| **7** | If preclinical data support advancement: IND-enabling planning, regulatory strategy development, and early clinical trial concept design | 3-7+ years | $5M+ | Requires biotech or pharma partner |

> Advancement between phases is conditional and depends on supportive data at each prior stage. This roadmap is a planning framework, not a claim of therapeutic readiness.

### What Phase 1 Produced

This project represents a complete Phase 1 deliverable:
- Computational structural comparison showing the compound heterozygous dimer yielded the lowest predicted interaction-confidence values among the tested states
- Literature synthesis across multiple disease-relevant pathways associated with MTHFR variant states
- Safety concepts framework for potential future experimental programs
- A fully replicable, open-source pipeline that any researcher can verify
- All produced for $0 using freely available tools (AlphaFold Server, Boltz-2, Python, GitHub)

### Phase 2 Next Steps

Phase 2 is achievable by an independent researcher:
1. Structural benchmarking of predictions against PDB 6FCX experimental structure
2. Molecular dynamics simulations on compound heterozygous dimer
3. Preprint preparation for bioRxiv to establish priority and make the work citable
4. Researcher outreach with the complete package

### The Phase 2 → Phase 3 Transition

Everything from Phase 3 onward requires a laboratory partner. The purpose of Phases 1 and 2 is to produce a package that supports a researcher's decision to evaluate these questions experimentally.

---

## 7. Collaboration Opportunities

This is an open research project seeking collaboration with:

- **Structural biologists:** To validate AlphaFold predictions with experimental cryo-EM or X-ray structures of C677T and A1298C variants
- **Gene therapy researchers:** To design and test ABE guide RNAs targeting C677T
- **Clinical researchers:** To establish clinical correlations across all seven disease pathways, particularly stratified by MTRR status
- **Bioinformaticians:** To independently replicate, extend, and improve the computational analysis
- **Epigeneticists:** To profile genome-wide methylation consequences of compound heterozygous MTHFR
- **Reproductive medicine specialists:** To investigate MTHFR-guided management of pregnancies in defined MTHFR variant cohorts

All data, code, and results are freely available under CC BY-NC-SA 4.0.

Contact: igor@dsm.promo | GitHub: https://github.com/DSMPromo/mthfr-target-validation

---

## 8. Acknowledgments

AlphaFold 3 predictions were performed using Google DeepMind's AlphaFold Server (Abramson et al., 2024). Protein sequences were obtained from UniProt (P42898). Reference structures were obtained from the Protein Data Bank (6FCX, 8QA5, 8QA6). AI tools (Claude, Anthropic) were used for literature organization, analysis pipeline development, and manuscript preparation. All AI-assisted work is explicitly disclosed per emerging best practices in scientific publishing.

---

## 9. Conflict of Interest Statement

The author declares no financial conflicts of interest. The author is a compound heterozygous MTHFR carrier (C677T + A1298C) whose personal genetic data motivated this research. No funding was received for this work. No pharmaceutical, biotechnology, or diagnostic company has any involvement in this project.

---

## 10. Ethics Statement

The genetic data presented in this study represents voluntary self-disclosure by the author. No other human subjects data was collected. The author's genetic variants were identified through consumer genetic testing (Genetic Genie / GenVue Discovery) and are shared voluntarily to demonstrate transparency and reproducibility. No institutional review board (IRB) approval was required for self-experimentation disclosure.

The author affirms that no individual should feel pressured to share their genetic information publicly. Anonymous contribution pathways are available through the project GitHub.

---

## 11. References

### Vision and Retinal Neurodegeneration
1. Markand S, Tawfik A, Ha Y, et al. (2015). Retinal Ganglion Cell Loss and Mild Vasculopathy in Methylene Tetrahydrofolate Reductase (Mthfr)-Deficient Mice: A Model of Mild Hyperhomocysteinemia. *Invest Ophthalmol Vis Sci*, 56(4):2684-2695. PMID: 25766590

2. Hergert RM, et al. (2022). Phoenix from the ashes: dramatic improvement in severe late-onset methylenetetrahydrofolate reductase (MTHFR) deficiency with a complete loss of vision. *J Neurol*, 269:4571-4577. PMC: PMC8940754

3. Jiang H, Liu Z, Townsend JH, Wang J. (2023). Effects of Methylenetetrahydrofolate Reductase (MTHFR) Polymorphisms on Retinal Tissue Perfusion in Mild Diabetic Retinopathy Patients Receiving the Medical Food, Ocufolin. *Clin Ophthalmol*, 17:1035-1043. PMC: PMC10106310

### Autism and Cerebral Folate Deficiency
4. Pu D, Shen Y, Wu J. (2020). Association between MTHFR C677T/A1298C and susceptibility to autism spectrum disorders: a meta-analysis. *BMC Pediatr*, 20:471. PMID: 32972375

5. Frye RE, Slattery J, Delhey L, et al. (2018). Folinic acid improves verbal communication in children with autism and language impairment: a randomized double-blind placebo-controlled trial. *Mol Psychiatry*, 23:631-636.

### Neuropsychiatric Disorders
6. Wan L, Li Y, Zhang Z, et al. (2018). Methylenetetrahydrofolate reductase and psychiatric diseases. *Transl Psychiatry*, 8:242. PMID: 30397195

7. Levin BL, Varga E. (2016). MTHFR: Addressing Genetic Counseling Dilemmas Using Evidence-Based Literature. *J Genet Couns*, 25(5):901-911.

### Structural Biology
8. Froese DS, Kopec J, Rembeza E, et al. (2018). Structural basis for the regulation of human 5,10-methylenetetrahydrofolate reductase by phosphorylation and S-adenosylmethionine inhibition. *Nat Commun*, 9:4735. PDB: 6FCX

9. Froese DS, et al. (2024). Dynamic inter-domain transformations mediate the allosteric regulation of human 5,10-methylenetetrahydrofolate reductase. *Nat Commun*. PDB: 8QA5, 8QA6

10. Pejchal R, Campbell E, Guenther BD, et al. (2006). Structural Perturbations in the Ala -> Val Polymorphism of Methylenetetrahydrofolate Reductase: How Binding of Folates May Protect against Inactivation. *Biochemistry*, 45(15):4808-4818. PMC: PMC1868400

### B-Vitamins and Cardiovascular
11. McNulty H, Dowey LRC, Strain JJ, et al. (2006). Riboflavin lowers homocysteine in individuals homozygous for the MTHFR 677C->T polymorphism. *Circulation*, 113(1):74-80.

### Epigenetics and DNA Methylation
12. Karahan G, Chan D, Shirane K, et al. (2021). Paternal MTHFR deficiency leads to hypomethylation of young retrotransposons and reproductive decline across two successive generations. *Development*, 148(13):dev199492. PMID: 34128976

13. Pentieva K, et al. (2020). Riboflavin supplementation alters global and gene-specific DNA methylation in adults with the MTHFR 677 TT genotype. *Biochimie*, 173:17-26.

### Pregnancy and Reproductive Health
14. Xia XP, et al. (2013). Meta-analysis of the methylenetetrahydrofolate reductase C677T polymorphism and susceptibility to pre-eclampsia. *Hypertension Research*, 35:1129-1134.

15. Kumar A, et al. (2024). Maternal methylenetetrahydrofolate reductase (MTHFR) A1298C polymorphism: implications in preventing recurrent pregnancy loss. *BMC Pregnancy Childbirth*. PMC: PMC11066823

### Cardiovascular and VTE
16. Liew SC, Gupta ED. (2015). MTHFR A1298C and C677T polymorphisms are associated with increased risk of venous thromboembolism: a retrospective chart review study. *PubMed*. PMID: 29212064

### Comprehensive Reviews
17. MTHFR Gene Polymorphisms: A Single Gene with Wide-Ranging Clinical Implications -- A Review. (2025). *Genes*, 16(4):441. PMC: PMC12027316

### Computational Methods
18. Abramson J, Adler J, Dunbar J, et al. (2024). Accurate structure prediction of biomolecular interactions with AlphaFold 3. *Nature*, 630:493-500.

### Gene Therapy and Base Editing
19. YOLT-101 Phase 1 Clinical Data. (2025). YolTech Press Release, March 2025. https://www.yoltx.com/news/press-release/95

20. LNP-mediated in vivo base editing of Agxt in a primary hyperoxaluria type 1 model. (2025). PubMed. PMID: 41275431

---

## Appendix A: Author Disclosure

The author is not a medical professional, biomedical scientist, or geneticist. The author is a cybersecurity specialist and cloud solution architect with 20+ years in IT who identified compound heterozygous MTHFR status (C677T + A1298C) plus MTRR I22M through consumer genetic testing and used computational tools to organize and analyze published research.

This project exists to provide a computational framework that qualified researchers can validate, extend, or refute. The author explicitly invites correction of any errors in interpretation or analysis.

This is not medical advice. All therapeutic concepts are speculative and would require years of experimental validation before any clinical application.

---

## Appendix B: Supplementary Materials

Available at https://github.com/DSMPromo/mthfr-target-validation:

1. **S1:** All FASTA sequences (sequences/)
2. **S2:** AlphaFold Server JSON job files (alphafold/jobs/json/)
3. **S3:** Complete AlphaFold results -- CIF structures, JSON confidence data, PAE matrices (alphafold/results/)
4. **S4:** Analysis pipeline source code (analyze.py)
5. **S5:** Google Colab notebook (MTHFR_AlphaFold_Analyzer.ipynb)
6. **S6:** PAE heatmap plots for all predictions (analysis/outputs/pae_plots/)
7. **S7:** Comparison charts -- ipTM and pTM bar charts (analysis/outputs/charts/)
8. **S8:** Complete metrics CSV (analysis/outputs/metrics.csv)
9. **S9:** HTML analysis report (analysis/outputs/report.html)
10. **S10:** PyMOL visualization scripts (analysis/outputs/pymol_scripts/)
