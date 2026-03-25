# MTHFR Gene Therapy Platform: One Gene, Five Disease Pathways, Billions of Lives

> **An open research project calling on real doctors and scientists to help turn a cybersecurity guy's genetic discovery into something that could change lives.**

[![License: CC BY-NC-SA 4.0](https://img.shields.io/badge/License-CC%20BY--NC--SA%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by-nc-sa/4.0/)
[![Status: Active Research](https://img.shields.io/badge/Status-Active%20Research-green.svg)]()
[![Looking for: PhD Collaborators](https://img.shields.io/badge/Looking%20for-PhD%20Collaborators-orange.svg)]()

---

## Who I Am (And Who I Am Not)

**I'm not a biohacker. I'm not trying to be something I'm not.**

My name is Igor Mihaljko. I'm a cybersecurity specialist, cloud solution architect, and infrastructure architect with 20+ years in IT and 14 Microsoft certifications. I run an AI marketing agency called DSM.Promo in Chicago.

I'm just a guy who thinks outside the box. I looked at my genetic test results, started asking questions that crossed the boundaries between medical specialties, and used AI tools to organize what I found. This project exists because I believe the next generation deserves better -- and sometimes it takes someone from outside a field to see what the people inside it can't.

**This project started with my own DNA.**

I took a consumer genetic test and discovered I carry two MTHFR mutations: C677T (heterozygous) and A1298C (heterozygous). This is called compound heterozygous status. When I started researching what this means, I found something that shocked me: the same broken enzyme that was potentially affecting my health is connected to five major disease categories -- and nobody had put the full picture together in one place.

I'm sharing everything I've found -- the science, the protein sequences, the structural modeling protocol, the analysis workflow -- so that real doctors and real scientists can take this further. I've done the computational groundwork using AI tools. I need people with labs, with expertise, with the ability to validate and build on this.

**This is my small contribution. Maybe together we can change lives for people who need it right now.**

---

## Important: Read This First

**This project is for research and educational purposes only.** Nothing here is medical advice. I am not qualified to give medical advice. All therapeutic concepts described are speculative and would need years of validation by qualified researchers before any human application.

The genomic data used in this project comes from consumer genetic testing (Genetic Genie / GenVue Discovery), which has no guarantees of clinical accuracy. Anyone interested in their MTHFR status should consult a qualified genetic counselor and get clinical-grade testing.

**But here's what I believe:** everyone should know their MTHFR status. A simple genetic test could explain years of unexplained symptoms -- anxiety that won't respond to medication, fatigue that doesn't make sense, vision problems, B-vitamin "deficiencies" that supplements don't seem to fix. Up to 40% of the world's population carries at least one MTHFR variant. Most don't know it.

---

## The Problem Nobody Connected

Modern medicine is organized by organ system. Ophthalmologists treat eyes. Psychiatrists treat anxiety. Neurologists study autism. Cardiologists manage heart disease. Nutritionists advise on B-vitamins.

**Nobody steps back to see that the same broken enzyme -- MTHFR -- feeds all five problems.**

I'm not the first person to study MTHFR. Thousands of papers have been published on individual aspects of this enzyme. But what I haven't found is a single document that connects all five disease pathways to one genetic root cause and proposes a unified correction strategy. That's what this project does.

## My DNA: The Key to This Research

Everything in this project traces back to my actual genetic variants, confirmed by consumer genetic testing:

| Gene | Variant | rsID | Zygosity | Classification | Relevance |
|------|---------|------|----------|----------------|-----------|
| **MTHFR** | **C677T (A222V)** | rs1801133 | **Hetero (AG)** | Clinically Significant | Catalytic domain -- FAD cofactor binding destabilized. All 5 disease targets. |
| **MTHFR** | **A1298C (E429A)** | rs1801131 | **Hetero (TG)** | Conflicting/Uncertain | Regulatory domain -- BH4/neurotransmitter pathway impaired. Anxiety, B-vitamins. |
| **MTRR** | **I22M (c.66A>G)** | rs1801394 | **Hetero (AG)** | Likely Pathogenic | B12 metabolism cofactor -- compounds the MTHFR methylation defect. |

> **Voluntary Disclosure:** I am sharing my genetic variant data voluntarily because I believe transparency strengthens this research. This is my personal decision. You should never feel pressured to share your genetic information publicly. If you want to contribute your own data to this project, you can do so anonymously.

**The compound heterozygous MTHFR status is the key.** C677T hits the catalytic domain (less folate processed). A1298C hits the regulatory domain (less neurotransmitter cofactor produced). Together, they create a dual failure that affects methylation, homocysteine, neurotransmitters, and B-vitamin processing simultaneously.

**This is not just my problem.** Up to 40% of the global population -- potentially 1.5 to 2 billion people -- carry at least one of these variants. Many don't know. Many are suffering from symptoms their doctors can't fully explain.

**Get tested.** It's a simple blood test or saliva test. Ask your doctor about MTHFR C677T and A1298C genotyping. Know your status.

---

## The Five Disease Pathways -- One Broken Enzyme

| # | Disease Target | How MTHFR Causes It | Who's Affected |
|---|---------------|---------------------|----------------|
| 1 | **Vision Loss** | Elevated homocysteine damages retinal cells and blood vessels in the eye | 2.2 billion with vision impairment (WHO) |
| 2 | **Autism / Cerebral Folate Deficiency** | Low methylfolate can't cross blood-brain barrier, developing brain starved of folate | 1 in 36 children have ASD (CDC) |
| 3 | **Anxiety and Depression** | Low SAMe leads to insufficient serotonin, dopamine, GABA; too much glutamate causes chronic anxiety | 280 million with anxiety disorders (WHO) |
| 4 | **B-Vitamin Processing Failure** | Can't convert B2, B6, B9, B12 to active forms even with normal diet | ~1.5-2 billion MTHFR carriers |
| 5 | **Heart Disease and Stroke** | Elevated homocysteine damages blood vessels throughout the body | 17.9 million CVD deaths/year (WHO) |

**The thesis:** One CRISPR base editing correction of MTHFR could address all five pathways through a single genetic intervention.

---

## What I've Built (And What I Need Help With)

### What's Done
- Comprehensive research document connecting all five disease pathways with published peer-reviewed evidence
- All four protein sequences (wild-type, C677T, A1298C, compound) ready for AlphaFold structural modeling
- Complete AlphaFold 3 submission plan (16 jobs across 3 days) -- step by step, anyone can follow
- Analysis workflow with exact PyMOL commands and Python scripts
- Metrics template for recording and comparing results
- Researcher outreach plan with specific collaboration targets

### What I Need From Real Scientists
- **Validation:** Are my interpretations of the published research correct?
- **Structural biology:** Do the AlphaFold predictions make biological sense?
- **Wet-lab work:** Can someone test base editor guide RNAs targeting C677T and A1298C?
- **Clinical perspective:** Is this therapeutically viable? What would a preclinical path look like?
- **Replication:** Can independent researchers reproduce the AlphaFold results?

**I have the AI tools and the systems thinking. I need the biology expertise. Let's work together.**

---

## Project Structure

```
mthfr-gene-therapy-project/
|-- README.md                    <-- You are here
|-- DISCLAIMER.md                <-- Important legal and medical disclaimers
|-- CONTRIBUTING.md              <-- How you can help
|-- LICENSE                      <-- CC BY-NC-SA 4.0
|
|-- docs/
|   |-- MTHFR_Definitive_Master_Document.docx
|       Complete research document with ALL science, ALL sequences,
|       ALL protocols, ALL analysis steps -- everything in one file
|
|-- sequences/                   <-- Ready to paste into AlphaFold Server
|   |-- MTHFR_wildtype.fasta     <-- Normal enzyme (656 amino acids)
|   |-- MTHFR_C677T.fasta        <-- A222V variant (catalytic domain)
|   |-- MTHFR_A1298C.fasta       <-- E429A variant (regulatory domain)
|   |-- MTHFR_compound.fasta     <-- Both mutations (my genotype)
|
|-- alphafold/
|   |-- jobs/
|   |   |-- submission_plan.md   <-- 16 jobs, step-by-step instructions
|   |-- results/                 <-- Put your downloaded results here
|
|-- analysis/
|   |-- analysis_workflow.md     <-- PyMOL commands, Python scripts
|   |-- metrics_template.csv     <-- Record your confidence scores here
|   |-- figures/                 <-- Save your comparison images here
|
|-- outreach/
    |-- target_researchers.md    <-- Who to contact for collaboration
    |-- email_template.md        <-- Ready-to-customize outreach email
```

## Quick Start: Replicate This Research

**Anyone with a Google account can do this. It's free. It takes 3 days.**

### Day 1: Run AlphaFold Jobs 1-6
1. Go to [alphafoldserver.com](https://alphafoldserver.com) and sign in with Google
2. Open `sequences/MTHFR_wildtype.fasta` -- copy the sequence (everything after the > header line)
3. In AlphaFold Server, click "Add entity" then select "Protein" and paste the sequence
4. Click "Add entity" again, select "Ligand", and choose "FAD" from the dropdown
5. Click "Continue and preview job", name it "Job01_WT_mono", and submit
6. Follow `alphafold/jobs/submission_plan.md` for all 6 Day 1 jobs

### Day 2: Replicate with New Seeds
Clone all 6 jobs using the three-dot menu and "Clone and reuse" option, then submit

### Day 3: Add Substrates
Run 4 more jobs adding folate (THF) and SAM ligands (see submission plan for details)

### After All Jobs Complete
Follow `analysis/analysis_workflow.md` to compare results in PyMOL

---

## Key References

1. Markand S, et al. (2015). *IOVS*, 56(4):2684-2695 -- Mthfr+/- mice: 20% ganglion cell loss, retinal vasculopathy
2. Hergert RM, et al. (2022). *J Neurology*, 269:4571-4577 -- Complete blindness reversed with betaine in MTHFR deficiency
3. Liu Z, et al. (2023). *Clin Ophthalmol*, 17:1035-1043 -- Retinal perfusion improved by MTHFR genotype
4. Pu D, et al. (2020). *BMC Pediatrics*, 20:471 -- Meta-analysis: MTHFR C677T significantly associated with ASD
5. Wan L, et al. (2018). *Transl Psychiatry*, 8:242 -- MTHFR and psychiatric diseases
6. Levin BL and Varga E. (2016). *Integr Med*, 15(2):48-52 -- Compound het anxiety treated with SAMe
7. Froese DS, et al. (2018). *Nat Commun*, 9:4735 -- Human MTHFR crystal structure, SAM regulation (PDB: 6FCX)
8. Pejchal R, et al. (2006). *Biochemistry*, 45(15):4808-4818 -- Structural perturbation of A222V: helix α5 displacement disrupts FAD binding
9. Froese DS, et al. (2024). *Nat Commun* -- Cryo-EM structures of MTHFR + SAM inhibited/active states (PDB: 8QA5, 8QA6)
10. Frye RE, et al. (2018). *Mol Psychiatry*, 23:631-636 -- Folinic acid improved verbal communication in ASD (randomized controlled trial)
11. Abramson J, et al. (2024). *Nature*, 630:493-500 -- AlphaFold 3: accurate structure prediction of biomolecular interactions

**Database Links:**
[UniProt P42898](https://www.uniprot.org/uniprotkb/P42898/entry) | [PDB 6FCX](https://www.rcsb.org/structure/6FCX) | [ClinVar rs1801133](https://www.ncbi.nlm.nih.gov/clinvar/?term=rs1801133) | [ClinVar rs1801131](https://www.ncbi.nlm.nih.gov/clinvar/?term=rs1801131)

---

## How You Can Help

See [CONTRIBUTING.md](CONTRIBUTING.md) for detailed instructions.

**If you're a scientist:** We need molecular biologists, ophthalmologists, neuroscientists, psychiatrists, structural biologists, and bioinformaticians. Open an issue or email me.

**If you're not a scientist:** Share this project. Get tested for MTHFR. Tell your doctor about it. The more people who know, the more pressure there is to take this seriously.

---

## License

[Creative Commons Attribution-NonCommercial-ShareAlike 4.0](https://creativecommons.org/licenses/by-nc-sa/4.0/)

---

## Contact

**Igor Mihaljko** | CEO and Founder, DSM.Promo | Chicago, IL | igor@dsm.promo

---

*I'm not a doctor. I'm not a biohacker. I'm a cybersecurity guy and cloud architect who looked at his DNA results, fell down a research rabbit hole, and realized that the dots nobody was connecting could potentially help billions of people. I don't have a lab. I don't have a PhD. I built this project to hand it to the real scientists. Everything is here -- the sequences, the protocols, the analysis pipeline. Take it. Build on it. Prove me wrong or prove me right.*

*If you're a researcher who can take this further, please reach out. If you're a parent of a child with autism, or someone living with unexplained anxiety, or someone losing their vision -- get tested for MTHFR. Know your status. It might explain everything.*

*I want to help future generations not have to go through what so many people are going through right now.*
