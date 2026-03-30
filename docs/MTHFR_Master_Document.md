# MTHFR Variant Hypothesis Prioritization Program

### Research Foundation and Computational Prioritization Framework for Selected MTHFR Variant States

**Research foundation, protein sequences, structural modeling workflow, and experimental follow-up framework for replication and interdisciplinary review**

Igor Mihaljko | Independent Researcher | Chicago, IL  
ORCID: [0009-0000-1408-1065](https://orcid.org/0009-0000-1408-1065)

March 2026 | Version 6.0 | FOR RESEARCH AND EDUCATIONAL PURPOSES ONLY

---

## Why This Document Exists

MTHFR biology sits at the intersection of folate metabolism, homocysteine regulation, methylation biology, neurotransmitter-relevant one-carbon metabolism, and vitamin-dependent pathway function. Because downstream findings are often discussed within separate specialties, relevant literature can appear fragmented. This document brings selected strands together into a single research framework and pairs them with a computational structural comparison designed to prioritize experimentally testable hypotheses.

Common MTHFR variants are widely distributed in human populations, and current clinical guidance does not treat most common variant states as stand-alone indications for major treatment changes. This document does not argue otherwise. Its purpose is narrower: to identify structural and biochemical questions that may warrant orthogonal validation.

Recent advances in structure prediction and related translational tools make it reasonable to explore whether selected MTHFR variant states merit deeper experimental study. In this document, computational modeling is used for hypothesis prioritization only. The goal is not to claim mechanism, disease causation, clinical utility, or therapeutic readiness.

---

## Core Boundaries of This Document

This document does not demonstrate reduced enzyme activity, altered FAD affinity in vitro, retinal injury, neuropsychiatric causation, or editing feasibility. Structural models described here are predictions, not experimentally resolved structures. Any future translational or therapeutic discussion remains downstream and contingent on substantial biochemical, structural, cellular, preclinical, regulatory, and clinical validation.

---

## PART 1: Scientific Context

### 1.1 The MTHFR Enzyme

MTHFR (methylenetetrahydrofolate reductase) is a FAD-dependent enzyme encoded on chromosome 1p36.22. It converts 5,10-methylenetetrahydrofolate to 5-methyltetrahydrofolate (5-MTHF), the primary circulating form of folate and the methyl donor for remethylation of homocysteine to methionine.

Methionine is converted to S-adenosylmethionine (SAMe), a major cellular methyl donor involved in DNA methylation, neurotransmitter-relevant pathways, myelin maintenance, hormone metabolism, and broader one-carbon biology.

MTHFR functions as a homodimer. The crystal structure (PDB: 6FCX, 2.5 A resolution, Froese et al. 2018, *Nature Communications*) revealed a catalytic TIM-barrel domain (N-terminal, FAD-binding) connected by a 25-amino-acid linker to a SAM-binding regulatory domain (C-terminal). In 2024, cryo-EM structures (PDB: 8QA5, 8QA6) further defined SAM-mediated allosteric inhibition through dual SAM binding and inter-domain conformational rearrangement, providing structural context for the regulatory-domain position of A1298C.

### 1.2 The Two Key Variants

| Property | C677T (rs1801133) | A1298C (rs1801131) |
|----------|-------------------|---------------------|
| Amino acid change | Ala222Val | Glu429Ala |
| Domain affected | Catalytic domain, FAD-related region | Regulatory domain, SAM-related region |
| Literature context | Altered enzyme stability, reduced FAD affinity, thermolabile behavior, elevated homocysteine in some contexts | Discussed in relation to regulatory behavior and BH4-relevant pathways, with mechanism and subgroup specificity remaining incompletely defined |
| Reported activity effects | Literature reports reduced activity, with strongest evidence in TT homozygotes | Literature reports activity effects, with interpretation varying across study design and context |
| Population frequency | Varies by population | Varies by population |
| Compound heterozygous context | Literature reports reduced overall enzyme activity relative to wild-type, with possible contributions from both catalytic and regulatory context |

Compound heterozygosity, one copy of each variant, has been reported in the literature to produce lower overall enzyme activity than wild-type, with effects that may involve both catalytic-domain and regulatory-domain context.

### 1.3 Primary Follow-Up Contexts

The present document prioritizes two downstream follow-up contexts from the broader literature.

#### Primary follow-up context: Retinal biomarkers

Prior literature links MTHFR-related one-carbon pathway disruption and hyperhomocysteinemia to retinal injury, making retinal follow-up a reasonable primary downstream context for experimental testing. Selected findings include:

- Mthfr+/- mice show increased retinal homocysteine, ganglion cell loss, thinner nerve fiber layer, and vascular leakage (Markand et al., 2015, *IOVS*)
- MTHFR modifies CRB1-related retinal phenotypes in the rd8 model (Markand et al., 2015, *Experimental Eye Research*)
- A case report described visual recovery associated with betaine treatment in severe MTHFR deficiency (Hergert et al., 2022, *Journal of Neurology*)
- Retinal perfusion improved in compound C677T/A1298C carriers receiving genotype-guided supplementation with Ocufolin medical food (Jiang et al., 2023, *Clinical Ophthalmology*)

In the present document, retinal biology is treated as a candidate biomarker domain for follow-up, not as an outcome established by the computational results.

#### Secondary exploratory context: Neuropsychiatric-relevant biochemistry

A1298C has been discussed in relation to BH4-relevant pathways, and one-carbon disruption has been discussed in relation to neurotransmitter-relevant biochemistry. However, mechanism, effect size, and subgroup specificity remain incompletely defined. Selected findings include:

- Meta-analysis has reported associations between MTHFR polymorphism and psychiatric phenotypes across several diagnostic categories (Wan et al., 2018, *Translational Psychiatry*)
- Compound heterozygous patients have been described in treatment-responsive clinical contexts involving SAMe and methylated B vitamins when MTHFR status is identified (Levin and Varga, 2016)
- These observations support exploratory biochemical follow-up in defined subgroups, rather than a primary disease-focused conclusion

In the present document, this context is treated as secondary and exploratory.

### 1.4 Additional Literature Context

Published literature also discusses MTHFR across autism and cerebral folate deficiency, B-vitamin processing, cardiovascular biology, DNA methylation, and reproductive health. In the present document, those broader domains are treated as background context rather than as outcomes established by the computational analysis.

### 1.5 Author Genotype Context

The following variants are reported from the author's consumer genomic testing. These data are included as personal context and hypothesis motivation only. Consumer genomic findings require confirmation in a clinical or research-grade setting.

| Gene | Variant | rsID | Zygosity | Classification | Context |
|------|---------|------|----------|----------------|---------|
| MTHFR | C677T (A222V) | rs1801133 | Hetero (AG) | Clinically Significant | Catalytic-domain variant with literature support for altered enzyme stability and reduced activity in some contexts |
| MTHFR | A1298C (E429A) | rs1801131 | Hetero (TG) | Conflicting or Uncertain | Regulatory-domain variant discussed in relation to one-carbon and neurotransmitter-relevant pathways, with mixed interpretation across sources |
| MTRR | I22M (c.66A>G) | rs1801394 | Hetero (AG) | Likely Pathogenic | Variant in downstream B12-related metabolism that could be relevant to one-carbon pathway context, pending phenotype-specific interpretation |

> Voluntary disclosure: These personal genomic results are shared voluntarily as context for the research question. They do not strengthen the computational evidence by themselves and should not be treated as clinical findings without independent confirmation.

### 1.6 MTRR Co-occurrence as a Contextual Factor

MTRR (methionine synthase reductase) regenerates the active form of vitamin B12, the essential cofactor for methionine synthase, the enzyme immediately downstream of MTHFR. Co-occurring variation in MTRR may therefore be relevant to one-carbon pathway context when present alongside selected MTHFR states.

In this document, MTRR I22M is included as contextual background rather than as a modeled result. The computational analysis does not determine the functional contribution of MTRR I22M, and it does not establish that co-occurring MTRR variation explains phenotype severity. However, this genotype combination may be relevant to variability among compound heterozygous individuals and may warrant stratified follow-up in future studies.

---

## PART 2: Computational Prioritization Framework

### 2.1 Purpose

Computational structure prediction, AlphaFold 3 Server and Boltz-2, was used to compare predicted structural confidence across wild-type, single-variant, and compound heterozygous MTHFR states in complex with FAD cofactor, THF substrate, and SAM allosteric inhibitor.

These predictions are used for hypothesis prioritization only. They do not establish mechanism, disease causation, clinical relevance, or therapeutic feasibility.

### 2.2 Key Reference Information

| Item | Value |
|------|-------|
| UniProt ID | P42898 |
| Protein length | 656 amino acids |
| Quaternary structure | Homodimer |
| Cofactor | FAD |
| C677T mutation site | Position 222, Ala to Val |
| A1298C mutation site | Position 429, Glu to Ala |
| PDB reference structure | 6FCX |
| Cryo-EM structures | 8QA5, 8QA6 |

### 2.3 Job Design

| Jobs | Protein | Configuration | Purpose | Platform |
|------|---------|--------------|---------|----------|
| 01-02 | Wild-type | Monomer + Dimer + FAD | Baseline | AlphaFold 3 |
| 03-04 | C677T (A222V) | Monomer + Dimer + FAD | Catalytic-domain variant | AlphaFold 3 |
| 05 | A1298C (E429A) | Monomer + FAD | Regulatory-domain variant | AlphaFold 3 |
| 06 | Compound + WT | Heterodimer + FAD | Heterozygous dimer context | AlphaFold 3 |
| 07-12 | All variants | Same as 01-06 | Independent seed replication | AlphaFold 3 |
| 13-15 | WT, C677T, Compound | Dimer + FAD + THF | Substrate-binding comparison | Boltz-2 |
| 16 | Wild-type | Dimer + FAD + SAM | Allosteric-inhibitor context | Boltz-2 |

### 2.4 Results Summary

#### AlphaFold 3 Server, Jobs 1 to 12

| Variant, dimer averaged | Avg pTM | Avg ipTM | Avg FAD Binding | Avg pLDDT@222 | Avg pLDDT@429 |
|-------------------------|---------|----------|-----------------|---------------|---------------|
| Wild-type | 0.775 | 0.740 | 0.555 | 97.3 | 96.0 |
| C677T | 0.785 | 0.765 | 0.575 | 97.05 | 95.95 |
| Compound heterozygous | 0.745 | 0.715 | 0.540 | 96.5 | 95.15 |

#### Boltz-2, Jobs 13 to 16

| Job | Variant | pTM | ipTM | Ligand ipTM | Protein ipTM |
|-----|---------|-----|------|-------------|-------------|
| 13 | WT + FAD + THF | 0.906 | 0.892 | 0.974 | 0.900 |
| 14 | C677T + FAD + THF | 0.891 | 0.875 | 0.969 | 0.883 |
| 15 | Compound + FAD + THF | 0.848 | 0.819 | 0.878 | 0.836 |
| 16 | WT + FAD + SAM | 0.911 | 0.897 | 0.925 | 0.899 |

### 2.5 Comparative Observations Within the Tested Model Set

1. Monomer predictions showed preserved overall folding across the tested variant states. The tested monomer outputs do not suggest large-scale loss of overall folding confidence.
2. Dimer predictions showed comparative inter-chain differences. FAD-associated confidence values were lower in the tested dimer models than in the tested monomer models, consistent with the different interaction context of the homodimer.
3. The compound heterozygous dimer yielded the lowest comparative confidence values across the assessed dimer-level metrics in both independent runs.
4. Compound heterozygous dimers averaged lower than the tested comparators, consistent with a possible combined dimer-level perturbation hypothesis in this modeling setup.
5. Independent random seeds produced consistent directional trends across the tested configurations, supporting limited internal reproducibility within this modeling setup.
6. In the Boltz-2 THF comparison, the compound heterozygous ligand ipTM was lower than wild-type, consistent with a possible difference in modeled substrate-associated context at the dimer level.

### 2.6 Core Observation

Monomer predictions showed preserved overall folding across the tested variant states, whereas dimer predictions showed comparative inter-chain differences. This is the central comparative finding of the current computational dataset and the main basis for downstream experimental prioritization.

### 2.7 What This Computational Section Does Not Show

This computational framework does not demonstrate reduced enzyme activity, altered FAD affinity in vitro, disease causation, retinal injury, neuropsychiatric causation, or therapeutic feasibility. It does not prove that observed confidence differences are biologically meaningful. Those questions require independent biochemical, structural, cellular, and clinical investigation.

---

## PART 3: Sequence Correction as a Downstream Translational Question

Sequence correction remains a downstream translational question rather than a conclusion of the present document. Because C677T is a single-nucleotide variant, base-editing logic can be discussed in principle. However, this document does not establish MTHFR as a therapeutic editing target and does not address guide design, edit-window fit, bystander edits, delivery, rescue, or safety.

Any serious consideration of sequence correction for MTHFR would require supportive evidence from prior stages, including biochemical validation, structural confirmation, cell-based studies, and target-specific safety assessment.

---

## PART 4: Development Roadmap

From computational observations to staged experimental follow-up:

| Phase | Focus | Timeline | Cost | Status |
|-------|-------|----------|------|--------|
| 1 | Computational hypothesis prioritization, AlphaFold 3 and Boltz-2 structural comparison, reproducible pipeline, open repository | Done | $0 | Complete |
| 2a | Molecular dynamics simulations: 10ns OpenMM/Amber14, WT dimer vs compound dimer. Compound dimer RMSD 7.34 A vs WT 5.29 A (p=1.05e-17). RMSF 6.6x higher at pos 222, 4.3x at pos 429. | Done | $0 | Complete |
| 2b | Preprint (bioRxiv submitted), researcher outreach, extended structural benchmarking | 1 to 3 months | Approx. $500 | In progress |
| 3 | Biochemical validation, expression of wild-type and selected variant proteins, dimer stability assays, FAD and THF interaction assays | 6 to 12 months | $50K to $150K | Requires lab partner |
| 4 | Exploratory biomarker correlation studies, including homocysteine, methylation panels, retinal OCT, and BH4-relevant readouts in well-defined cohorts | 12 to 18 months | $100K to $300K | Requires clinical collaborator |
| 5 | If biochemical and biomarker findings support continued investigation, sequence-level guide assessment, edit-window analysis, bystander-risk assessment, and cell-based feasibility studies | 18 to 24 months | $200K to $500K | Requires gene-editing partner |
| 6 | If cell-based feasibility and safety data support further progression, preclinical animal studies using appropriate Mthfr-relevant models and multi-organ outcome measures | 24 to 36 months | $500K to $2M | Requires institutional partner |
| 7 | If preclinical data support advancement, IND-enabling planning, regulatory strategy development, and early clinical trial concept design | 3 to 7+ years | $5M+ | Requires biotech or pharma partner |

> Advancement between phases is conditional and depends on supportive data at each prior stage. This roadmap is a planning framework, not a claim of therapeutic readiness.

---

## PART 5: References

### Vision and Retinal Biology
1. Markand S, et al. (2015). *Invest Ophthalmol Vis Sci*, 56(4):2684-2695. Mthfr+/- mice: ganglion cell loss and retinal vasculopathy
2. Hergert RM, et al. (2022). *J Neurol*, 269:4571-4577. Visual recovery associated with betaine in severe MTHFR deficiency
3. Jiang H, Liu Z, et al. (2023). *Clin Ophthalmol*, 17:1035-1043. Retinal perfusion and MTHFR genotype

### Neuropsychiatric and Cerebral Folate Context
4. Pu D, et al. (2020). *BMC Pediatr*, 20:471. Meta-analysis: C677T associated with ASD
5. Frye RE, et al. (2018). *Mol Psychiatry*, 23:631-636. Folinic acid randomized trial in ASD
6. Wan L, et al. (2018). *Transl Psychiatry*, 8:242. MTHFR and psychiatric disorders
7. Levin BL, Varga E. (2016). *J Genet Couns*, 25(5):901-911. MTHFR genetic counseling and clinical evidence

### Structural Biology
8. Froese DS, et al. (2018). *Nat Commun*, 9:4735. Human MTHFR crystal structure, PDB 6FCX
9. Froese DS, et al. (2024). *Nat Commun*. Cryo-EM MTHFR allosteric regulation, PDB 8QA5 and 8QA6
10. Pejchal R, et al. (2006). *Biochemistry*, 45(15):4808-4818. A222V structural perturbation

### B Vitamins, Epigenetics, and Broader Context
11. McNulty H, et al. (2006). *Circulation*, 113(1):74-80. Riboflavin and homocysteine in 677TT
12. Karahan G, et al. (2021). *Development*, 148(13):dev199492. Transgenerational sperm DNA methylation
13. Pentieva K, et al. (2020). *Biochimie*, 173:17-26. Riboflavin alters DNA methylation in 677TT
14. Xia XP, et al. (2013). *Hypertension Research*. C677T and preeclampsia
15. Kumar A, et al. (2024). *BMC Pregnancy Childbirth*. A1298C and recurrent pregnancy loss
16. MTHFR Gene Polymorphisms. (2025). *Genes*, 16(4):441. Comprehensive review

### Computational Methods and Editing Context
17. Abramson J, et al. (2024). *Nature*, 630:493-500. AlphaFold 3
18. YOLT-101 Phase 1 Clinical Data. (2025). YolTech. Next-generation adenine base editor

---

## Disclaimer

This document is a computational and literature-based research planning document. It does not provide medical advice, diagnosis, or treatment recommendations. Structural models described here are predictions, not experimentally resolved structures. Disease links discussed in this document reflect varying levels of published evidence and should not be interpreted as proof of causation, clinical utility, or therapeutic readiness. Any therapeutic concept described here remains speculative and would require substantial experimental, preclinical, regulatory, and clinical validation.

Genomic data referenced here are from consumer genomic services and require professional confirmation. AlphaFold Server is for non-commercial use only.

---

## Contact

Igor Mihaljko  
Independent Researcher  
Chicago, IL  
igor@dsm.promo
