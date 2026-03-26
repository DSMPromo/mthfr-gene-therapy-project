# MTHFR Variant Target Validation Program

### Integrated Research Foundation and Computational Prioritization Framework for Selected MTHFR Variant States

**Research foundation, protein sequences, structural modeling workflow, and collaboration framework for replication and experimental follow-up**

Igor Mihaljko | CEO & Founder, DSM.Promo | Chicago, IL
ORCID: [0009-0000-1408-1065](https://orcid.org/0009-0000-1408-1065)

March 2026 | Version 5.0 | FOR RESEARCH AND EDUCATIONAL PURPOSES ONLY

---

## Why These Findings Merit Cross-Disciplinary Review

MTHFR biology sits at the intersection of folate metabolism, homocysteine regulation, methylation, neurotransmitter synthesis, and vitamin-dependent one-carbon metabolism. Because downstream effects are often studied within separate clinical specialties, findings relevant to MTHFR-associated biology are frequently discussed in fragmented ways. This document brings those strands together into a single research framework and pairs them with a computational structural comparison designed to prioritize experimentally testable hypotheses.

Common MTHFR variants are widely distributed in human populations, and current clinical guidance does not treat most common variant states as stand-alone indications for major treatment changes. For example, the CDC states that common MTHFR variants such as C677T are not a reason to avoid folic acid, and notes limited evidence that A1298C alone significantly affects folate processing.

Recent advances in structure prediction and gene-editing make it reasonable to explore whether selected MTHFR variant states warrant deeper experimental study. In this document, computational modeling is used for hypothesis prioritization only. The goal is not to claim therapeutic readiness, but to identify structural and biochemical questions suitable for orthogonal validation.

---

## PART 1: The Science

### 1.1 The MTHFR Enzyme

MTHFR (methylenetetrahydrofolate reductase) is a FAD-dependent enzyme encoded on chromosome 1p36.22. It converts 5,10-methylenetetrahydrofolate to 5-methyltetrahydrofolate (5-MTHF), the primary circulating form of folate and the methyl donor for remethylation of homocysteine to methionine.

Methionine is converted to S-adenosylmethionine (SAMe), the universal methyl donor for over 200 reactions: DNA methylation, neurotransmitter synthesis (serotonin, dopamine, norepinephrine, GABA), myelin maintenance, hormone metabolism, and detoxification.

The enzyme functions as a homodimer. The crystal structure (PDB: 6FCX, 2.5 A resolution, Froese et al. 2018, *Nature Communications*) reveals a catalytic TIM-barrel domain (N-terminal, FAD binding) connected by a 25-amino-acid linker to a SAM-binding regulatory domain (C-terminal) that mediates dimerization and allosteric inhibition. In 2024, cryo-EM structures (PDB: 8QA5, 8QA6) revealed the mechanism of SAM-mediated allosteric inhibition through dual SAM binding and inter-domain conformational rearrangement, providing structural context for the A1298C variant's position in the regulatory domain.

### 1.2 The Two Key Variants

| Property | C677T (rs1801133) | A1298C (rs1801131) |
|----------|-------------------|---------------------|
| Amino acid change | Ala222Val (A→V at position 222) | Glu429Ala (E→A at position 429) |
| Domain affected | Catalytic domain (FAD binding site) | Regulatory domain (SAM binding / BH4-relevant pathways) |
| Enzyme activity | Literature reports ~35% reduction per allele; ~70% in TT homozygotes | Literature reports activity reduction, with interpretation varying across study design and context |
| Primary association | Altered enzyme stability, reduced FAD affinity, thermolabile enzyme, elevated homocysteine | Discussed in relation to BH4-relevant pathways and one-carbon metabolism, with mechanism and subgroup specificity remaining incompletely defined |
| T allele frequency | ~50% Latino, 24-40% European, 26-37% East Asian, ~6.6% African | ~25% European, varies by population |
| Compound heterozygous | Literature reports reduced overall enzyme activity relative to wild-type, with effects that may involve both catalytic and regulatory context |

### 1.3 Associated Disease-Relevant Pathways

Selected MTHFR variant states have been associated in the literature with multiple downstream biological domains. These associations vary in strength, and causation is not established uniformly across all conditions.

**Pathway 1: Retinal and Vascular Biology (Primary Experimental Focus)**

Prior literature has linked MTHFR-related one-carbon pathway disruption and hyperhomocysteinemia to retinal injury. Key published findings:

- Mthfr+/- mice show 2x retinal homocysteine, ~20% ganglion cell loss, thinner nerve fiber layer, and vascular leakage, without elevated IOP (Markand et al., 2015, *IOVS*)
- MTHFR modifies CRB1-related retinopathies, worsening rd8 retinal phenotype (Markand et al., 2015, *Exp Eye Res*)
- Case report: betaine treatment was associated with visual recovery in a patient with severe MTHFR deficiency (Hergert et al., 2022, *J Neurology*)
- Retinal perfusion was improved in compound C677T/A1298C carriers receiving MTHFR genotype-guided supplementation with Ocufolin medical food (Jiang et al., 2023, *Clin Ophthalmol*)

**Pathway 2: Neuropsychiatric Biology (Primary Experimental Focus)**

A1298C has been discussed in relation to BH4-relevant pathways. BH4 is the cofactor for tryptophan hydroxylase (serotonin synthesis) and tyrosine hydroxylase (dopamine synthesis). However, the mechanism, effect size, and subgroup specificity remain incompletely defined. Key published findings:

- MTHFR polymorphism is associated with schizophrenia, major depression, autism, bipolar disorder, and ADHD (Wan et al., 2018, *Transl Psychiatry*)
- Compound heterozygous patients show treatment-responsive symptoms with SAMe and methylated B vitamins when MTHFR status is identified (Levin & Varga, 2016)
- This may warrant exploratory study in defined subgroups where one-carbon metabolism and BH4-relevant biology are hypothesized to be relevant

**Pathway 3: Autism Spectrum Disorder and Cerebral Folate Deficiency**

Reduced 5-MTHF transport across the blood-brain barrier via folate receptor alpha is associated with cerebral folate deficiency. Meta-analysis of 15 studies supports a significant association between C677T and ASD across all five genetic models (Pu et al., 2020). A randomized controlled trial showed that folinic acid improved verbal communication in ASD children with language impairment (Frye et al., 2018).

**Pathway 4: B-Vitamin Processing**

The C677T variant has reduced FAD affinity. B2 supplementation partially compensates by increasing cellular FAD availability. Riboflavin supplementation lowered homocysteine specifically in 677TT homozygotes (McNulty et al., 2006, *Circulation*).

**Pathway 5: Cardiovascular Biology**

Elevated homocysteine is associated with increased cardiovascular risk. Homocysteine-related vascular mechanisms observed in retinal literature may also be relevant to systemic vascular biology, though this project does not test that directly.

**Pathway 6: Genome-Wide DNA Methylation**

In mouse models, Mthfr deficiency has been associated with transgenerational sperm DNA methylation effects persisting across F1 and F2 generations (Karahan et al., 2021, *Development*). Riboflavin supplementation alters global and gene-specific DNA methylation in adults with the MTHFR 677TT genotype (Pentieva et al., 2020). These findings support further investigation of epigenetic consequences.

**Pathway 7: Pregnancy and Reproductive Health**

MTHFR variants have been associated with preeclampsia (meta-analysis of 51 studies, Xia et al., 2013), recurrent pregnancy loss (Kumar et al., 2024), and neural tube defects. The reproductive pathway is notable because a carrier parent can pass both the genetic variant and potential epigenetic consequences to offspring, which warrants investigation as a potential transgenerational mechanism.

---

## PART 2: Author Genotype Context

The following variants are reported from the author's consumer genomic testing (Genetic Genie / GenVue Discovery). These data are included as personal context and hypothesis motivation only. Consumer genomic findings require confirmation in a clinical or research-grade setting.

| Gene | Variant | rsID | Zygosity | Classification | Context |
|------|---------|------|----------|----------------|---------|
| MTHFR | C677T (A222V) | rs1801133 | Hetero (AG) | Clinically Significant | Catalytic-domain variant with literature support for altered enzyme stability and reduced activity in some contexts |
| MTHFR | A1298C (E429A) | rs1801131 | Hetero (TG) | Conflicting/Uncertain | Regulatory-domain variant discussed in relation to one-carbon and neurotransmitter-relevant pathways, with mixed clinical interpretation across sources |
| MTRR | I22M (c.66A>G) | rs1801394 | Hetero (AG) | Likely Pathogenic | Variant in downstream B12-related metabolism that could be relevant to one-carbon pathway context, pending phenotype-specific interpretation |

> **Voluntary Disclosure:** I am sharing my genetic variant data voluntarily because I believe transparency strengthens this research. This is my personal decision. You should never feel pressured to share your genetic information publicly.

### The MTRR Context

MTRR (methionine synthase reductase) regenerates the active form of vitamin B12, the essential cofactor for methionine synthase, the enzyme immediately downstream of MTHFR. This combination could place additive pressure on one-carbon metabolism by affecting upstream folate handling and downstream B12-dependent recycling. Whether this co-occurrence contributes to phenotypic variability among compound heterozygous individuals warrants investigation.

---

## PART 3: AlphaFold 3 and Boltz-2 Structural Modeling

### 3.1 Purpose

Computational structure prediction (AlphaFold 3 Server and Boltz-2 via Tamarind Bio) was used to compare predicted structural confidence across wild-type, single-variant, and compound heterozygous MTHFR states in complex with FAD cofactor, THF substrate, and SAM allosteric inhibitor.

These predictions are used for hypothesis prioritization only. They do not establish mechanism, disease causation, or therapeutic feasibility.

### 3.2 Key Reference Information

| Item | Value |
|------|-------|
| UniProt ID | P42898 |
| Protein length | 656 amino acids |
| Quaternary structure | Homodimer |
| Cofactor | FAD (flavin adenine dinucleotide) |
| C677T mutation site | Position 222: Ala → Val |
| A1298C mutation site | Position 429: Glu → Ala |
| PDB reference structure | 6FCX (2018, X-ray, 2.5 A) |
| Cryo-EM structures | 8QA5 (SAH dis-inhibited), 8QA6 (SAM inhibited) |

### 3.3 Job Design

| Jobs | Protein | Configuration | Purpose | Platform |
|------|---------|--------------|---------|----------|
| 01-02 | Wild-type | Monomer + Dimer + FAD | Baseline | AlphaFold 3 |
| 03-04 | C677T (A222V) | Monomer + Dimer + FAD | Catalytic domain variant | AlphaFold 3 |
| 05 | A1298C (E429A) | Monomer + FAD | Regulatory domain variant | AlphaFold 3 |
| 06 | Compound + WT | Heterodimer + FAD | Compound heterozygous state | AlphaFold 3 |
| 07-12 | All variants | Same as 01-06 | Replication with independent seeds | AlphaFold 3 |
| 13-15 | WT, C677T, Compound | Dimer + FAD + THF | Substrate binding comparison | Boltz-2 |
| 16 | Wild-type | Dimer + FAD + SAM | Allosteric inhibitor binding | Boltz-2 |

### 3.4 Results Summary

**AlphaFold 3 Server (Jobs 1-12):**

| Variant (Dimer, averaged) | Avg pTM | Avg ipTM | Avg FAD Binding | Avg pLDDT@222 | Avg pLDDT@429 |
|---------------------------|---------|----------|-----------------|---------------|---------------|
| Wild-type | 0.775 | 0.740 | 0.555 | 97.3 | 96.0 |
| C677T | 0.785 | 0.765 | 0.575 | 97.05 | 95.95 |
| Compound heterozygous | 0.745 | 0.715 | 0.540 | 96.5 | 95.15 |

**Boltz-2 (Jobs 13-16):**

| Job | Variant | pTM | ipTM | Ligand ipTM | Protein ipTM |
|-----|---------|-----|------|-------------|-------------|
| 13 | WT + FAD + THF | 0.906 | 0.892 | 0.974 | 0.900 |
| 14 | C677T + FAD + THF | 0.891 | 0.875 | 0.969 | 0.883 |
| 15 | Compound + FAD + THF | 0.848 | 0.819 | 0.878 | 0.836 |
| 16 | WT + FAD + SAM | 0.911 | 0.897 | 0.925 | 0.899 |

**Key observations:**

1. Monomer predictions showed preserved overall folding across the tested variant states (ipTM 0.97-0.98); the tested monomer predictions do not suggest large-scale loss of overall folding confidence
2. Dimer predictions show inter-chain differences: FAD binding confidence drops from 0.97 (monomer) to 0.53-0.58 (dimer), consistent with known cooperativity in the homodimer
3. The compound heterozygous dimer yielded the lowest comparative confidence values across the assessed metrics in both independent runs
4. Compound heterozygous dimers averaged lower than the tested comparators (ipTM 0.715 vs WT 0.740 and C677T 0.765), consistent with a possible combined dimer-level perturbation in this modeling setup
5. Independent random seeds produced consistent directional trends across the tested configurations, supporting limited internal reproducibility within this modeling setup
6. THF substrate binding: compound het ligand ipTM 0.878 vs 0.974 WT, consistent with reduced substrate access at the dimer level

**Core observation:** Monomer predictions showed preserved overall folding across the tested variant states, whereas dimer predictions show inter-chain differences. This is the central comparative finding of the current computational dataset and the main basis for downstream experimental prioritization.

---

## PART 4: Safety Concepts for Future Experimental Programs

Any future gene-editing program would need to address safety. The following are examples of safety layers used or discussed in broader editing research that could be relevant to MTHFR correction if it advances to experimental testing. Each would require target-specific validation.

| Layer | Concept | Precedent |
|-------|---------|-----------|
| 1. Inherent design | Base editing avoids DNA double-strand breaks; LNP-delivered mRNA degrades within 48 hours | Used in VERVE-102, BEAM-101 |
| 2. Off-target screening | Cas-OFFinder, GUIDE-seq, whole-genome sequencing before any intervention | Standard in current editing programs |
| 3. Anti-CRISPR proteins | Emergency neutralization of residual editor activity (e.g., AcrIIA4) | Demonstrated in preclinical models |
| 4. Reversible epigenetic editing | CRISPRoff as a non-permanent first step; potentially reversible under laboratory conditions | In early translational development |
| 5. Post-intervention monitoring | Biomarker panels, off-target sequencing at defined timepoints | Standard clinical protocol design |

At the target base, a candidate correction strategy would aim to restore the reference allele sequence rather than introduce a novel coding change.

At the sequence level, C677T appears theoretically compatible with adenine base editing (ABE) logic, which converts A-T base pairs to G-C base pairs. However, edit window, guide design, bystander edits, tissue targeting, and functional rescue would all need empirical validation.

---

## PART 5: Development Roadmap

From computational observations to experimental validation:

| Phase | Focus | Timeline | Cost | Status |
|-------|-------|----------|------|--------|
| 1 | Computational hypothesis prioritization, AlphaFold 3 and Boltz-2 structural comparison, reproducible pipeline, open repository | Done | $0 | Complete |
| 2 | Structural benchmarking against experimental reference structures such as PDB 6FCX, molecular dynamics simulations, preprint preparation, researcher outreach | 1-3 months | ~$500 | Next |
| 3 | Biochemical validation: expression of wild-type and selected variant proteins, dimer stability assays, FAD and THF interaction assays | 6-12 months | $50K-150K | Requires lab partner |
| 4 | Exploratory biomarker correlation studies, including homocysteine, methylation panels, retinal OCT, and BH4-related readouts in well-defined compound heterozygous cohorts | 12-18 months | $100K-300K | Requires clinical collaborator |
| 5 | If biochemical and biomarker findings support continued investigation: sequence-level guide assessment, edit-window analysis, bystander-risk assessment, and cell-based feasibility studies | 18-24 months | $200K-500K | Requires gene-editing partner |
| 6 | If cell-based feasibility and safety data support further progression: preclinical animal studies using appropriate Mthfr-relevant models and multi-organ outcome measures | 24-36 months | $500K-2M | Requires institutional partner |
| 7 | If preclinical data support advancement: IND-enabling planning, regulatory strategy development, and early clinical trial concept design | 3-7+ years | $5M+ | Requires biotech or pharma partner |

> Advancement between phases is conditional and depends on supportive data at each prior stage. This roadmap is a planning framework, not a claim of therapeutic readiness.

---

## PART 6: Researcher Outreach

### Target Researchers

| Researcher | Institution | Expertise | Rationale |
|------------|------------|-----------|-----------|
| Dr. Sylvia Smith | Augusta University | MTHFR + retinal degeneration | Published the Mthfr+/- mouse retina papers |
| Dr. Richard Frye | Phoenix Children's Hospital | Folate metabolism + ASD | Leading autism/CFD researcher, folinic acid trials |
| Dr. D. Sean Froese | U. of Zurich | MTHFR crystal structure | Solved PDB 6FCX, can validate predictions against experimental data |
| Dr. Kiran Musunuru | U. of Pennsylvania | CRISPR base editing in vivo | Pioneer of LNP-delivered base editing |
| Dr. David Liu | Broad Institute / Harvard | Invented base editors | Created ABE/CBE technology |

---

## References

### Vision and Retinal Biology
1. Markand S, et al. (2015). *Invest Ophthalmol Vis Sci*, 56(4):2684-2695. Mthfr+/- mice: ganglion cell loss, retinal vasculopathy
2. Hergert RM, et al. (2022). *J Neurol*, 269:4571-4577. Case report: visual recovery associated with betaine in severe MTHFR deficiency
3. Jiang H, Liu Z, et al. (2023). *Clin Ophthalmol*, 17:1035-1043. Retinal perfusion and MTHFR genotype

### Autism and Cerebral Folate Deficiency
4. Pu D, et al. (2020). *BMC Pediatr*, 20:471. Meta-analysis: C677T associated with ASD
5. Frye RE, et al. (2018). *Mol Psychiatry*, 23:631-636. Folinic acid RCT in ASD

### Neuropsychiatric
6. Wan L, et al. (2018). *Transl Psychiatry*, 8:242. MTHFR and psychiatric disorders
7. Levin BL, Varga E. (2016). *J Genet Couns*, 25(5):901-911. MTHFR genetic counseling

### Structural Biology
8. Froese DS, et al. (2018). *Nat Commun*, 9:4735. Human MTHFR crystal structure (PDB: 6FCX)
9. Froese DS, et al. (2024). *Nat Commun*. Cryo-EM: MTHFR + SAM allosteric regulation (PDB: 8QA5, 8QA6)
10. Pejchal R, et al. (2006). *Biochemistry*, 45(15):4808-4818. A222V structural perturbation

### B-Vitamins and Cardiovascular
11. McNulty H, et al. (2006). *Circulation*, 113(1):74-80. Riboflavin and homocysteine in 677TT

### Epigenetics
12. Karahan G, et al. (2021). *Development*, 148(13):dev199492. Transgenerational sperm DNA methylation
13. Pentieva K, et al. (2020). *Biochimie*, 173:17-26. Riboflavin alters DNA methylation in 677TT

### Pregnancy
14. Xia XP, et al. (2013). *Hypertension Research*. C677T and preeclampsia
15. Kumar A, et al. (2024). *BMC Pregnancy Childbirth*. A1298C and recurrent pregnancy loss

### Reviews
16. MTHFR Gene Polymorphisms (2025). *Genes*, 16(4):441. Comprehensive review

### Computational Methods
17. Abramson J, et al. (2024). *Nature*, 630:493-500. AlphaFold 3

### Gene Editing Precedent
18. YOLT-101 Phase 1 Clinical Data (2025). YolTech. Next-gen adenine base editor

---

## Disclaimer

This document is a computational and literature-based research planning document. It does not provide medical advice, diagnosis, or treatment recommendations. Structural models described here are predictions, not experimentally resolved structures. Disease links discussed in this document reflect varying levels of published evidence and should not be interpreted as proof of causation, clinical utility, or therapeutic readiness. Any therapeutic concept described here remains speculative and would require substantial experimental, preclinical, regulatory, and clinical validation.

Genomic data are from consumer genomic services and require professional confirmation. AlphaFold Server is for non-commercial use only.
