# AlphaFold Server Submission Plan

> **16 jobs across 3 days.** Follow exactly for reproducible results.

## Before You Start

1. Go to [alphafoldserver.com](https://alphafoldserver.com)
2. Sign in with your Google account
3. You get **20 jobs per day** (resets at midnight BST)
4. Each job takes **3–8 minutes** and produces **5 ranked predictions**

## Important Notes

- **Verify mutation positions** against UniProt P42898 canonical sequence before submitting
- The C677T variant (p.Ala222Val) position should be verified by searching for "A" at the expected location
- Set seed to **"Auto"** for all jobs unless doing exact replication
- The FASTA files in `sequences/` folder contain the ready-to-paste sequences

---

## DAY 1: Core Comparisons (6 jobs)

### Job 1: Wild-Type Monomer + FAD
| Entity | Type | Input | Copies |
|--------|------|-------|--------|
| 1 | Protein | Paste `MTHFR_wildtype.fasta` sequence | 1 |
| 2 | Ligand | Select **FAD** from dropdown | 1 |

**Purpose:** Baseline single-chain structure with FAD cofactor bound.

---

### Job 2: Wild-Type Homodimer + FAD
| Entity | Type | Input | Copies |
|--------|------|-------|--------|
| 1 | Protein | Paste `MTHFR_wildtype.fasta` sequence | **2** |
| 2 | Ligand | Select **FAD** | **2** |

**Purpose:** Native functional dimer — this is how the enzyme works in your body.

---

### Job 3: C677T Variant Monomer + FAD
| Entity | Type | Input | Copies |
|--------|------|-------|--------|
| 1 | Protein | Paste `MTHFR_C677T.fasta` sequence | 1 |
| 2 | Ligand | Select **FAD** | 1 |

**Purpose:** See how the A222V mutation destabilizes FAD binding.

---

### Job 4: C677T Homodimer + FAD (TT genotype)
| Entity | Type | Input | Copies |
|--------|------|-------|--------|
| 1 | Protein | Paste `MTHFR_C677T.fasta` sequence | **2** |
| 2 | Ligand | Select **FAD** | **2** |

**Purpose:** Models the most severe homozygous TT individuals.

---

### Job 5: A1298C Variant Monomer + FAD
| Entity | Type | Input | Copies |
|--------|------|-------|--------|
| 1 | Protein | Paste `MTHFR_A1298C.fasta` sequence | 1 |
| 2 | Ligand | Select **FAD** | 1 |

**Purpose:** See how E429A affects the regulatory domain.

---

### Job 6: Compound Heterozygous Dimer — YOUR GENOTYPE ⭐
| Entity | Type | Input | Copies |
|--------|------|-------|--------|
| 1 | Protein | Paste `MTHFR_compound.fasta` sequence | 1 |
| 2 | Protein | Paste `MTHFR_wildtype.fasta` sequence | 1 |
| 3 | Ligand | Select **FAD** | **2** |

**Purpose:** THIS MODELS YOUR ACTUAL ENZYME. One mutant chain + one wild-type chain as a heterodimer. This is what's happening in the cells of every compound heterozygous carrier.

---

## DAY 2: Replication Seeds (6 jobs)

Clone Jobs 1–6 using the **"Clone and reuse"** button (three-dot menu next to each job in History).

1. Click the three-dot menu → "Clone and reuse"
2. In the "Confirm and submit" dialog, keep seed on **"Auto"**
3. Submit without changes
4. Repeat for all 6 jobs

**Purpose:** Second independent prediction for statistical comparison. Compare ranking_score between original and clone. Use the higher-ranked model.

---

## DAY 3: Substrate & Inhibitor Binding (4 jobs)

### Job 13: WT Dimer + FAD + Folate Substrate
| Entity | Type | Input | Copies |
|--------|------|-------|--------|
| 1 | Protein | Paste `MTHFR_wildtype.fasta` | **2** |
| 2 | Ligand | **FAD** | **2** |
| 3 | Ligand | **THF** (tetrahydrofolate) | **2** |

---

### Job 14: C677T Dimer + FAD + Folate
| Entity | Type | Input | Copies |
|--------|------|-------|--------|
| 1 | Protein | Paste `MTHFR_C677T.fasta` | **2** |
| 2 | Ligand | **FAD** | **2** |
| 3 | Ligand | **THF** | **2** |

---

### Job 15: Compound Dimer + FAD + Folate ⭐
| Entity | Type | Input | Copies |
|--------|------|-------|--------|
| 1 | Protein | Paste `MTHFR_compound.fasta` | 1 |
| 2 | Protein | Paste `MTHFR_wildtype.fasta` | 1 |
| 3 | Ligand | **FAD** | **2** |
| 4 | Ligand | **THF** | **2** |

---

### Job 16: WT Dimer + FAD + SAM (Allosteric Inhibitor)
| Entity | Type | Input | Copies |
|--------|------|-------|--------|
| 1 | Protein | Paste `MTHFR_wildtype.fasta` | **2** |
| 2 | Ligand | **FAD** | **2** |
| 3 | Ligand | **SAM** (S-adenosylmethionine) | **2** |

**Purpose:** Model allosteric inhibition. The SAM-binding regulatory domain is where A1298C strikes.

---

## After All Jobs Complete

1. Download all ZIP files
2. Organize into folders (see project structure in README)
3. Follow `analysis/analysis_workflow.md` for next steps
4. Record all metrics in `analysis/metrics_template.csv`
