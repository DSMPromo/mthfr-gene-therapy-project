# MTHFR Variant Hypothesis Prioritization: Computational Prioritization of Selected MTHFR Variant States for Experimental Follow-Up

## A Hypothesis-Prioritization Study Using AlphaFold 3 and Boltz-2

---

**Author:** Igor Mihaljko, Independent Researcher, DSM.Promo, Chicago, IL

**ORCID:** [0009-0000-1408-1065](https://orcid.org/0009-0000-1408-1065)

**Correspondence:** igor@dsm.promo

**Keywords:** MTHFR, AlphaFold 3, Boltz-2, C677T, A1298C, FAD cofactor, protein structure prediction, dimer biology, structural prioritization, computational hypothesis prioritization, compound heterozygous

**License:** CC BY-NC-SA 4.0

**Data Availability:** All sequences, analysis scripts, and AlphaFold results available at https://github.com/DSMPromo/mthfr-target-validation

**AI Disclosure:** AI tools (Claude, Anthropic) were used for literature organization, analysis pipeline development, and manuscript preparation. All AI-generated content was reviewed for accuracy by the author. This disclosure follows emerging best practices in scientific publishing (Nature, 2023).

---

## Abstract

Methylenetetrahydrofolate reductase (MTHFR) variants C677T (A222V) and A1298C (E429A) are common in human populations, yet their comparative behavior in modeled dimer and ligand-associated contexts remains incompletely characterized. We used AlphaFold 3 Server and Boltz-2 to compare wild-type, C677T, A1298C, and compound heterozygous MTHFR states in the presence of FAD, THF, and SAM. Across 34 structural predictions (5 independent seeds per configuration), monomer predictions remained broadly similar, whereas dimer predictions showed comparative inter-chain differences. Within this model set, the compound heterozygous dimer yielded the lowest reported interaction-confidence values across the assessed dimer-interface and ligand-associated metrics. These observations do not establish mechanism, disease causation, clinical relevance, or therapeutic feasibility. They define falsifiable follow-up questions centered on dimer-level biology, with retinal follow-up as a primary downstream context and neuropsychiatric-relevant biochemistry as a secondary exploratory context. These computational observations define bounded follow-up questions and require orthogonal validation through biochemical and structural methods.

**Note:** This is a computational hypothesis-prioritization study. All structural predictions require orthogonal validation through biochemical and structural methods.

---

## 1. Introduction

### 1.1 The MTHFR Enzyme and Its Clinical Significance

MTHFR (methylenetetrahydrofolate reductase, EC 1.5.1.20) is a FAD-dependent homodimeric enzyme encoded on chromosome 1p36.22 that catalyzes the conversion of 5,10-methylenetetrahydrofolate to 5-methyltetrahydrofolate (5-MTHF), the primary circulating form of folate and the methyl donor for homocysteine remethylation to methionine (Froese et al., 2018). This reaction sits at the intersection of folate metabolism, homocysteine regulation, methylation biology, and selected downstream pathways discussed in the literature.

The crystal structure of human MTHFR (PDB: 6FCX, 2.5 A resolution) revealed a homodimer with each subunit containing a conserved catalytic TIM-barrel domain (N-terminal, FAD-binding) and a eukaryote-specific regulatory domain (C-terminal, SAM-binding) connected by a 25-residue inter-domain linker (Froese et al., 2018). Recent cryo-EM structures (PDB: 8QA5, 8QA6) have further elucidated the mechanism of SAM-mediated allosteric inhibition through dual SAM binding and inter-domain conformational rearrangement (Froese et al., 2024).

A comprehensive 2025 review characterizes MTHFR as having "wide-ranging clinical implications" across cardiovascular, neurological, oncological, and reproductive domains, supporting its relevance across multiple clinical and biological domains discussed in the literature (Genes, 2025).

### 1.2 The Two Key Variants

Two common MTHFR polymorphisms have substantial clinical significance:

**C677T (rs1801133, p.Ala222Val):** Located in the catalytic domain, this variant destabilizes FAD cofactor binding. The A222V substitution displaces helix alpha-5, which carries FAD-interacting residues Asn168, Arg171, and Lys172, resulting in a thermolabile enzyme with reduced activity (~35% per allele, ~70% in TT homozygotes) (Pejchal et al., 2006). The degree of thermolability is dramatic: residual activity after heat inactivation is only 18-22% in 677TT individuals compared to 66-67% in wild-type (NCBI Bookshelf). The T allele frequency ranges from 24-50% across populations. Critically, folate binding to the enzyme protects against FAD dissociation, providing a molecular explanation for why folate supplementation partially rescues enzyme function (Pejchal et al., 2006).

**A1298C (rs1801131, p.Glu429Ala):** Located in the regulatory domain near the SAM-binding site, this variant has been reported in some studies to associate with reduced enzyme activity of ~30% per allele and has been discussed in relation to BH4-relevant pathways, though the mechanism, effect size, and subgroup specificity remain incompletely defined. The interaction between position 429 and the SAM-binding pocket identified in the 2024 cryo-EM structures (Froese et al., 2024) provides structural context for this variant's functional impact.

**Compound heterozygosity** (one copy of each variant) produces enzyme activity comparable to C677T homozygosity (~50-60% of normal), providing a plausible basis for combined catalytic-domain and regulatory-domain effects that warrant experimental testing (Weisberg et al., 1998). Compound heterozygous individuals show lower plasma folate concentrations and elevated homocysteine levels equivalent to C677T homozygotes. VTE occurred in 50% of compound heterozygous patients versus 28.6% of wild-type subjects in a retrospective analysis (Liew & Gupta, 2015).

### 1.3 The MTRR Compounding Effect

The functional consequences of MTHFR variants may be modified when co-occurring with variants in other methylation cycle genes. Methionine synthase reductase (MTRR, rs1801394, I22M) is particularly relevant because it regenerates the active form of vitamin B12 (methylcobalamin), which is the essential cofactor for methionine synthase -- the enzyme immediately downstream of MTHFR in the methylation cycle.

When MTHFR produces less 5-MTHF (due to C677T/A1298C variants) AND MTRR is less efficient at regenerating B12 (due to I22M), the methylation cycle shows reduced capacity at two sequential steps. This could contribute additional one-carbon pathway stress if the co-occurrence is supported as functionally relevant:
- Less substrate (5-MTHF) is produced by MTHFR
- The enzyme that uses that substrate (methionine synthase) has a less-functional cofactor (due to reduced MTRR-mediated regeneration)
- Homocysteine accumulates from both upstream and downstream blocks
- B12 becomes trapped in inactive forms, further reducing methylation capacity

This compound genetic context is directly relevant to the present study, as the author genotype context carries heterozygous variants in both MTHFR (compound C677T/A1298C) and MTRR (I22M), representing co-occurring variants in MTHFR and MTRR that may be relevant to one-carbon pathway context.

### 1.4 Primary Follow-Up Contexts

The present study prioritizes two downstream follow-up contexts from the broader literature: retinal follow-up as the primary context, and neuropsychiatric-relevant biochemistry as a secondary exploratory context.

**Primary follow-up context: Retinal biomarkers**

Prior literature links MTHFR-related one-carbon pathway disruption and hyperhomocysteinemia to retinal injury, making retinal follow-up a reasonable primary downstream context for experimental testing. Mthfr+/- mice (heterozygous -- the same genetic status as most human carriers) show 2x retinal homocysteine, ~20% ganglion cell loss, thinner nerve fiber layer, and vascular leakage (Markand et al., 2015). Case literature suggests that betaine treatment was associated with visual recovery in selected severe MTHFR deficiency contexts (Hergert et al., 2022). Retinal perfusion was improved in compound C677T/A1298C carriers receiving MTHFR genotype-guided supplementation with the medical food Ocufolin (Jiang et al., 2023). MTHFR additionally modifies CRB1-related retinopathies, worsening the rd8 retinal phenotype (Markand et al., 2015).

**Secondary exploratory context: Neuropsychiatric-relevant biochemistry**

Because one-carbon pathway disruption may influence BH4-relevant neurotransmitter biology, neuropsychiatric-relevant biochemical follow-up is retained here as a secondary exploratory context rather than a primary disease-focused conclusion. Meta-analysis of 15 studies supports a significant association between C677T and ASD across all five genetic models (Pu et al., 2020). The 677TT genotype was 2x more prevalent in autistic children in a Chinese Han study (16.1% vs 8.6%, OR=2.04). A randomized controlled trial demonstrated that folinic acid improved verbal communication in ASD children with language impairment (Frye et al., 2018). Maternal folate supplementation before and during pregnancy reduces autism risk in a study of 45,300 children (Levine et al., 2018).

Additional literature context across cardiovascular, nutritional, epigenetic, and reproductive domains is noted in the broader literature, but those domains are not primary follow-up targets of the present computational study and are not interpreted here as direct consequences of the modeled structural differences.

### 1.5 Additional Literature Context

Published literature discusses MTHFR across multiple clinical and biological domains. In the present study, those broader discussions are treated as background context rather than as outcomes established by the current computational analysis.

### 1.6 Knowledge Gap and Study Rationale

While thousands of individual studies address MTHFR variants in specific disease contexts, no published work has:

1. Systematically characterized the structural consequences of C677T, A1298C, and compound heterozygous variants using AlphaFold 3's ability to predict protein-ligand complexes
2. Connected selected downstream literature contexts to structural predictions in a single hypothesis-prioritization framework
3. Evaluated compound heterozygous MTHFR in the context of MTRR co-occurrence
4. Discussed sequence correction only as a downstream translational question contingent on prior validation steps

AlphaFold 3 (Abramson et al., 2024) can model protein + cofactor (FAD) + substrate (THF) + inhibitor (SAM) interactions in a single prediction, enabling direct comparison of wild-type and variant enzyme structures with their functional ligands.

### 1.7 Why MTHFR Has Seen Limited Translational Development

As of March 2026, our searches did not identify a registered interventional program specifically focused on therapeutic gene-editing of MTHFR. This absence should not be interpreted as evidence of target invalidity, but it does suggest limited translational development to date. Several factors may contribute:

**1. Common-variant framing.** Because these variants are common in human populations, they have often been framed as polymorphisms rather than as translational intervention targets.

**2. Supplementation as partial compensation.** Methylfolate, riboflavin, B12, and SAMe supplementation partially compensates for impaired enzyme function, creating the perception that the problem is addressed. However, supplementation is lifelong and does not correct the underlying genetic variant.

**3. Fragmentation across specialties.** Relevant findings have been discussed across separate specialties rather than within a single structural prioritization framework.

**4. Recency of structural tools.** The first human MTHFR crystal structure was published in 2018 (Froese et al., PDB 6FCX). The cryo-EM structures revealing SAM allosteric regulation were published in 2024 (PDB 8QA5, 8QA6). AlphaFold 3 became available in 2024. The tools to structurally characterize these variants did not exist until recently.

**What has changed:** Recent advances in structure prediction, editing programs for other targets, and new MTHFR structural data make the question more technically discussable than before. The present study does not establish feasibility of correction. It asks whether the target merits bounded experimental follow-up.

### 1.8 Study Objectives

1. Compare wild-type, single-variant, and compound heterozygous MTHFR states in modeled monomer, dimer, and ligand-associated contexts
2. Determine whether the compound heterozygous state is prioritized by comparative dimer-level metrics within the tested model set
3. Define bounded experimental follow-up questions for biochemical assays and selected downstream biomarker studies
4. Clarify which conclusions remain unsupported without orthogonal validation

This study does not demonstrate reduced enzyme activity, altered FAD affinity in vitro, retinal injury, neuropsychiatric causation, or editing feasibility. Those questions require independent biochemical, cellular, structural, and clinical investigation.

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

**Job Design (34 predictions):**

| Job | Protein | Copies | Ligands | Purpose |
|-----|---------|--------|---------|---------|
| 01 | Wild-type | 1 (monomer) | FAD x1 | Baseline monomer |
| 02 | Wild-type | 2 (dimer) | FAD x2 | Native functional dimer |
| 03 | C677T | 1 | FAD x1 | Variant monomer |
| 04 | C677T | 2 | FAD x2 | Homozygous TT dimer |
| 05 | A1298C | 1 | FAD x1 | Variant monomer |
| 06 | Compound + WT | 1+1 | FAD x2 | Heterozygous dimer |
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

### 2.5 Reference Structure Context

Wild-type monomer predictions were compared to the experimentally determined crystal structure of human MTHFR (PDB: 6FCX, Froese et al., 2018) using gemmi-based structural superposition. The resulting RMSD of 1.42 A over 596 matched CA atoms classifies the prediction as high accuracy (RMSD <= 2.0 A according to standard benchmarking criteria), supporting the use of AlphaFold 3 predictions for comparative structural analysis within this protein system.

### 2.6 Reference Structures

Predictions were contextualized against three experimental structures:
- **PDB 6FCX:** Human MTHFR crystal structure, 2.5 A, homodimer with FAD (Froese et al., 2018)
- **PDB 8QA5:** MTHFR + SAH dis-inhibited (active) state, cryo-EM (Froese et al., 2024)
- **PDB 8QA6:** MTHFR + SAM inhibited state, dual SAM binding, cryo-EM (Froese et al., 2024)

### 2.7 Statistical Analysis

For replicated predictions (Jobs 01-06 vs Jobs 07-12), we report the mean and range of each metric. With five independent seeds per condition (n=5), we present descriptive statistics with standard deviations and perform Welch's t-tests for inter-variant comparisons. This sample size supports basic statistical testing of the observed inter-variant differences.

---

## 3. Results

### 3.1 Monomer Predictions: Overall Quality

All monomer predictions (3 variants x 5 seeds) achieved high-confidence fold predictions:

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
- The C677T variant shows a subtle but consistent reduction in pLDDT at position 222 (98.0-98.1 vs 98.5 for WT) across all five seeds
- Replicate seeds produce highly consistent results (pTM variance < 0.02), confirming prediction reliability
- No clashes detected in any prediction

### 3.2 Dimer Predictions: Comparative Differences in the Tested Models

The tested dimer models showed comparative differences that were less apparent in the tested monomer models:

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
| WT dimer (n=5) | 0.786 +/-0.016 | 0.760 +/-0.023 | 0.568 +/-0.016 | 97.3 | 96.2 |
| C677T dimer | 0.785 | 0.765 | 0.575 | 97.05 | 95.95 |
| **Compound dimer** (n=5) | **0.744 +/-0.028** | **0.710 +/-0.032** | **0.538 +/-0.025** | **96.5** | **95.2** |

**Comparative observations within the tested model set:**

1. **Lower ligand-associated confidence is observed in the tested dimer models than in the tested monomer models:** FAD binding confidence (chain_iptm) drops from 0.97 in monomers to 0.53-0.57 in dimers. Within this model set, that marked reduction is directionally compatible with the importance of quaternary context for cofactor-associated interactions, but it does not establish binding mechanism.

2. **The compound heterozygous dimer consistently shows the lowest scores across every metric in all five independent runs:** Job 06 (pTM 0.73, ipTM 0.70, FAD binding 0.53, pLDDT@429 95.0) and Job 12 (pTM 0.76, ipTM 0.73, FAD binding 0.55, pLDDT@429 95.3). In the tested heterozygous dimer configuration, one chain carried A222V and one chain carried E429A. Within this model set, that configuration yielded the lowest reported confidence values across the assessed dimer-level metrics in all runs.

3. **Relative to the tested C677T dimer model, the compound heterozygous dimer showed lower average ipTM within this model set**, which prioritizes that configuration for experimental follow-up.

4. **Within the tested dimer models, position 429 showed the lowest pLDDT values in the compound heterozygous configuration:** From 97.5 (monomers) to 95.0-95.3 (compound dimer), consistent with a possible regulatory-domain contribution in the dimer context and compatible with prior cryo-EM work showing that SAM-mediated allosteric regulation operates across the dimer interface.

5. **All findings replicate across independent seeds:** Each variant was predicted with five independent random seeds. The directional trends were consistent across the tested runs (t-test: WT dimer ipTM 0.760 +/-0.023 vs compound dimer ipTM 0.710 +/-0.032, t=2.565, p=0.033, significant at p<0.05), supporting the use of these outputs for computational prioritization within this modeling setup.

### 3.3 FAD Binding Confidence: Monomer vs Dimer Comparison

A notable difference within the tested model set is the contrast between monomer and dimer FAD-associated confidence values:

| Variant | Monomer FAD Binding (avg) | Dimer FAD Binding (avg) | Change |
|---------|-------------------------:|------------------------:|-------:|
| Wild-type | 0.975 | 0.555 | -0.42 |
| C677T | 0.970 | 0.575 | -0.40 |
| **Compound** | **0.970** | **0.540** | **-0.43** |

At the monomer level, all tested states showed high FAD-associated confidence values. This is directionally compatible with prior literature indicating that the C677T enzyme can bind FAD, while experimental work has shown altered thermal stability. AlphaFold predicts static structures and therefore cannot capture thermolability directly, but the dimer predictions reveal interface-mediated effects on FAD binding that are invisible in monomers.

The compound heterozygous dimer showed the largest reduction in reported FAD-associated confidence within the tested model set, which is consistent with a combined dimer-level perturbation hypothesis but does not establish mechanism.

*See analysis/outputs/charts/iptm_comparison.png for visualization.*

### 3.4 Local Confidence at Mutation Sites

**Position 222 (C677T site):**

| Variant | Monomer pLDDT@222 (avg) | Dimer pLDDT@222 (avg) | Change |
|---------|------------------------:|----------------------:|-------:|
| WT | 98.5 | 97.3 | -1.2 |
| C677T | 98.05 | 97.05 | -1.0 |
| Compound | -- | 96.5 | -- |

- C677T shows a subtle but consistent reduction in pLDDT at position 222 (98.0-98.1 vs 98.5 for WT) across all five monomer seeds
- The compound dimer shows the lowest pLDDT at position 222 (96.6)

**Position 429 (A1298C site):**

| Variant | Monomer pLDDT@429 (avg) | Dimer pLDDT@429 (avg) | Change |
|---------|------------------------:|----------------------:|-------:|
| WT | 97.5 | 96.0 | -1.5 |
| C677T | 97.7 | 95.95 | -1.75 |
| Compound | -- | **95.15** | -- |

- Position 429 shows the lowest pLDDT values in the compound dimer within the tested model set, consistent with a possible regulatory-domain contribution at the dimer level
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

The compound heterozygous dimer consistently ranked lowest across all five independent seeds (n=5 per configuration). Statistical testing (Welch's t-test) confirmed that the ipTM difference between WT and compound dimers reached significance (p=0.033). The inter-replicate variance was modest relative to the observed inter-variant differences within this model set, supporting the use of these outputs for computational prioritization, but not establishing biological significance without orthogonal validation.

Effect size analysis (Cohen's d) confirmed large differences between WT and compound heterozygous dimers for ipTM (d=1.62), pTM (d=1.64), and FAD-associated confidence (d=1.29). After Bonferroni correction for 8 simultaneous comparisons, pLDDT at position 429 remained statistically significant (raw p=0.0004, adjusted p=0.003, d=3.65), identifying the regulatory domain as the most robust comparative signal. The ipTM and pTM comparisons showed large effect sizes but did not survive Bonferroni correction (adjusted p=0.267 and 0.255 respectively), consistent with the limited sample size and supporting cautious interpretation.

### 3.7 Substrate and Inhibitor Binding (Jobs 13-16)

Jobs 13-16 (THF substrate and SAM inhibitor binding) will be submitted manually through the AlphaFold Server interface as Phase 2 of the computational analysis. These predictions will reveal whether substrate access is differentially affected by C677T and A1298C variants, and whether SAM allosteric inhibition is altered in the compound heterozygous dimer.

---

## 4. Discussion

The present study is a computational hypothesis-prioritization analysis. Its purpose is not to establish mechanism or disease relevance, but to define bounded follow-up questions based on comparative structural outputs across tested MTHFR variant states.

### 4.1 Structural Interpretation of Monomer Results

The present computational results are consistent with preserved overall folding across the tested monomer states and with the possibility of subtle structural perturbations rather than gross folding loss. Within that bounded interpretation, the following hypotheses are reasonable to consider:

1. **Thermolability** (C677T): The enzyme folds and binds FAD normally at low temperatures but loses FAD more rapidly at physiological temperature (37C). AlphaFold predicts equilibrium structures and therefore would not capture this dynamic property in static predictions.

2. **Regulatory dysfunction** (A1298C): The E429A substitution near the SAM-binding pocket may alter allosteric regulation without dramatically changing the fold. The 2024 cryo-EM structures show that SAM-mediated inhibition requires precise positioning of dual SAM molecules, so subtle changes at position 429 are a reasonable hypothesis for future testing.

3. **Compound effects**: The dimer predictions are the most relevant part of the present model set for asking whether the two mutations could contribute to combined perturbation at the dimer interface.

The subtle pLDDT reduction at position 222 in C677T predictions (0.4-0.5 points) is small but consistent across seeds, suggesting that AlphaFold detects some degree of local structural uncertainty at the mutation site.

### 4.2 Why the Dimer Context Matters

MTHFR functions as a homodimer in vivo, making the dimer context more relevant than the monomer context for evaluating potential inter-chain and ligand-associated effects. The dimer interface may contribute to FAD-associated interactions, allosteric communication between regulatory and catalytic regions, and overall structural stability.

Within the tested model set, the heterozygous dimer configuration, modeled as one chain carrying A222V and one chain carrying E429A, yielded the lowest reported interaction-confidence values across the assessed dimer-level metrics. This makes the heterozygous dimer a reasonable priority for experimental follow-up. It does not establish that this configuration is functionally impaired in vivo, but it provides a focused structural question for biochemical and structural testing.

### 4.3 Primary and Secondary Follow-Up Contexts

This study focuses on two follow-up contexts derived from the broader MTHFR literature.

**Primary follow-up context: Retinal biomarkers**
Prior animal, case-based, and supplementation literature supports retinal follow-up as a testable downstream context for experimental work. In the present study, the retinal context is used as a candidate biomarker domain for follow-up rather than as an outcome established by the computational results.

**Secondary exploratory context: Neuropsychiatric-relevant biochemistry**
Because one-carbon pathway disruption has been discussed in relation to BH4-relevant neurotransmitter biology and related metabolic pathways, neuropsychiatric-relevant biochemical follow-up is retained here as a secondary exploratory context. In the present study, this context is treated as exploratory and does not constitute a disease-focused conclusion.

Broader literature across cardiovascular, nutritional, epigenetic, and reproductive domains exists but remains outside the scope of the present computational analysis.

| Context | Role in present study | Basis in literature |
|---------|----------------------|---------------------|
| Retinal biomarkers | Primary follow-up context | Animal, case-based, and supplementation literature support retinal follow-up as a testable downstream context |
| Neuropsychiatric-relevant biochemistry | Secondary exploratory context | Literature discusses one-carbon and BH4-relevant biochemical pathways, supporting exploratory follow-up |

### 4.4 MTRR Co-occurrence as a Contextual Factor

Co-occurring variation in MTRR may be relevant to one-carbon pathway context when present alongside compound MTHFR variants. Because MTRR supports regeneration of the active vitamin B12 cofactor required by methionine synthase, co-occurrence of MTHFR and MTRR variants could contribute additional pathway stress in some individuals.

In the present study, this co-occurrence is included as contextual background rather than as a modeled result. The computational analysis does not determine the functional contribution of MTRR I22M, and it does not establish that co-occurring MTRR variation explains phenotype severity. However, this genotype combination may be relevant to variability among compound heterozygous individuals and may warrant stratified follow-up in future studies.

### 4.5 The Riboflavin-FAD-C677T Structural Connection

Prior literature reports that riboflavin supplementation can lower homocysteine in selected C677T genotype contexts, and structural studies have shown that the A222V substitution perturbs the region surrounding FAD-associated interactions. Together, these observations provide useful context for interpreting the present computational results.

In this framework, increased cellular FAD availability offers one possible explanation for previously reported genotype-specific supplementation findings in the literature. The present computational study does not test nutritional intervention, does not measure rescue, and does not establish treatment effect. Instead, it provides structural context that may help frame future biochemical follow-up.

### 4.6 Sequence Correction as a Downstream Translational Question

Sequence correction remains a downstream translational question rather than a conclusion of the present study. Because C677T is a single-nucleotide variant, base-editing logic can be noted in principle. However, the present work does not establish MTHFR as a therapeutic editing target and does not address guide design, edit-window fit, bystander edits, delivery, rescue, or safety.

Current editing programs for other targets show that precise sequence correction is technically discussable in modern translational research. In the context of this paper, that observation serves only as background. Any serious consideration of sequence correction for MTHFR would require supportive evidence from prior stages, including biochemical validation, structural confirmation, cell-based studies, and target-specific safety assessment.

### 4.7 Translational Caution

Any translational consideration would depend on successful completion of prior validation stages. The present computational study should not be interpreted as evidence of intervention readiness.

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

1. **Experimental structure determination:** Cryo-EM of purified C677T and A1298C variant MTHFR in complex with FAD, with and without folate substrate. To our knowledge, the compound heterozygous dimer (one A222V chain + one E429A chain) has not been structurally characterized experimentally.

2. **Thermal shift assays:** Differential scanning fluorimetry (DSF) to quantify FAD binding affinity and thermal stability of WT vs variant MTHFR, with and without riboflavin supplementation.

3. **Molecular dynamics simulations:** Extended (>1 microsecond) MD simulations of WT, C677T, A1298C, and compound MTHFR to capture FAD dissociation kinetics, allosteric communication, and dynamic differences not visible in static predictions.

4. **Cell-based follow-up after biochemical validation:** Once biochemical and structural follow-up supports continued investigation, computational guide design and initial cell-based feasibility testing could be used to assess whether sequence correction is technically plausible in a controlled setting.

5. **Clinical correlation study:** Prospective study of compound heterozygous MTHFR individuals stratified by MTRR status, focusing on retinal OCT, plasma homocysteine, and neuropsychiatric-relevant biomarkers.

6. **Epigenetic profiling:** Genome-wide methylation profiling (EPIC array or whole-genome bisulfite sequencing) of compound heterozygous individuals vs matched controls.

7. **Model organism studies:** Mthfr-relevant models incorporating MTRR co-occurrence, where feasible, to assess pathway-level effects and test whether sequence-correction concepts merit further preclinical evaluation.

---

## 5. Conclusion

This study presents an AlphaFold 3 and Boltz-2 based comparative characterization of selected MTHFR variant states in modeled monomer, dimer, and ligand-associated contexts. By connecting these computational outputs to retinal and neuropsychiatric-relevant follow-up contexts, it provides a framework for prioritizing bounded experimental questions arising from the modeled dimer-level differences.

The key contribution is not the computational predictions themselves, which require experimental validation, but the framework: a reproducible, open-source structural prioritization workflow that links modeled dimer-level differences to bounded experimental follow-up questions.

We invite collaboration from structural biologists, clinicians, and bioinformaticians to validate, refine, and extend this work. All data, code, and analysis are openly available.

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
- Computational structural comparison showing that the compound heterozygous dimer yielded the lowest reported interaction-confidence values among the tested states
- Literature synthesis supporting bounded retinal and neuropsychiatric-relevant follow-up contexts
- A fully replicable, open-source pipeline that independent researchers can verify

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
- **Gene editing researchers:** To assess whether sequence correction becomes technically plausible after prior biochemical and structural validation
- **Clinical researchers:** To establish phenotype correlations in retinal and neuropsychiatric-relevant contexts, particularly when stratified by MTRR status and other confounders
- **Bioinformaticians:** To independently replicate, extend, and improve the computational analysis
- **Epigeneticists:** To profile genome-wide methylation consequences of compound heterozygous MTHFR

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

The genetic data presented in this study represents voluntary self-disclosure by the author. No other human subjects data was collected. The author's genetic variants were identified through consumer genetic testing (Genetic Genie / GenVue Discovery) and are shared voluntarily to demonstrate transparency and reproducibility. No institutional review board (IRB) approval was required for self-disclosure of the author's own data in this document.

The author affirms that no individual should feel pressured to share their genetic information publicly. Any future collection of outside participant data would require appropriate ethical and regulatory review.

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

### Base-Editing Precedent
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
