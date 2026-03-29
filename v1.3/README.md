# MTHFR Gene Therapy Project - v1.3

## Purpose

v1.3 extends the AlphaFold Server dimer modeling campaign to 10 independent seed runs for WT and compound heterozygous dimers, and 7 runs for C677T dimers. Runs 1-5 were generated in earlier versions; this batch adds runs 6-10 (WT/compound) and runs 6-7 (C677T).

## Jobs (12 total)

| # | File | Variant | Type | Chains | FAD |
|---|------|---------|------|--------|-----|
| 1 | wt_dimer_run6.json | Wild-type | Homodimer | WT x2 | 2 |
| 2 | wt_dimer_run7.json | Wild-type | Homodimer | WT x2 | 2 |
| 3 | wt_dimer_run8.json | Wild-type | Homodimer | WT x2 | 2 |
| 4 | wt_dimer_run9.json | Wild-type | Homodimer | WT x2 | 2 |
| 5 | wt_dimer_run10.json | Wild-type | Homodimer | WT x2 | 2 |
| 6 | compound_dimer_run6.json | Compound (C677T+A1298C) | Heterodimer | COMPOUND x1 + WT x1 | 2 |
| 7 | compound_dimer_run7.json | Compound (C677T+A1298C) | Heterodimer | COMPOUND x1 + WT x1 | 2 |
| 8 | compound_dimer_run8.json | Compound (C677T+A1298C) | Heterodimer | COMPOUND x1 + WT x1 | 2 |
| 9 | compound_dimer_run9.json | Compound (C677T+A1298C) | Heterodimer | COMPOUND x1 + WT x1 | 2 |
| 10 | compound_dimer_run10.json | Compound (C677T+A1298C) | Heterodimer | COMPOUND x1 + WT x1 | 2 |
| 11 | c677t_dimer_run6.json | C677T (homozygous) | Homodimer | C677T x2 | 2 |
| 12 | c677t_dimer_run7.json | C677T (homozygous) | Homodimer | C677T x2 | 2 |

## Mutations

- **C677T**: A222V substitution in the catalytic domain (position 222: A -> V)
- **A1298C**: E429A substitution in the regulatory domain (position 429: E -> A)
- **Compound**: Both A222V and E429A on the same chain, paired with one WT chain (heterodimer)

## Folder Structure

```
v1.3/
  jobs/
    json/           # 12 individual job files + ALL_12_v1.3_JOBS.json
  results/          # AlphaFold Server output (to be populated)
  README.md
```

## Combined File

`ALL_12_v1.3_JOBS.json` contains all 12 jobs in a single array for batch reference.

## Cumulative Run Counts After v1.3

| Variant | Total Runs |
|---------|-----------|
| WT dimer | 10 (runs 1-10) |
| Compound heterodimer | 10 (runs 1-10) |
| C677T dimer | 7 (runs 1-7) |
