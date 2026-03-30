# Changelog

## v1.3 (March 29, 2026) -- COMPLETE
- Extended to 10 independent seeds per configuration (64 total AlphaFold predictions)
- 3 of 4 key metrics now survive Bonferroni correction (was 1 of 4 with 5 seeds)
- ipTM: p=0.003, Bonferroni p=0.035 (significant)
- pTM: p=0.003, Bonferroni p=0.031 (significant)
- pLDDT@429: p<0.000001, Bonferroni p=0.000005 (highly significant)
- 10ns molecular dynamics simulations completed (OpenMM, Amber14, Google Colab A100)
- MD results: compound dimer RMSD 7.34 vs WT 5.29 A (p=1.05e-17)
- RMSF at pos 222: 9.76 A vs 1.47 A (6.6x more flexible in compound)
- RMSF at pos 429: 6.96 A vs 1.61 A (4.3x more flexible in compound)
- Published: Zenodo DOI 10.5281/zenodo.19318627
- Published: Preprints.org ID 205673
- ORCID linked: 0009-0000-1408-1065
- Researcher outreach emails sent (Dr. Froese, Dr. Smith)
- Consolidated all results into single alphafold/results_all/ folder (64 results + 4 Boltz-2)
- Added run_md.py, validate_structure.py, MTHFR_MD_Simulation.ipynb
- Split all figures into individual small files (no oversized images)

## v1.2 (March 29, 2026)
- Extended to 5 independent seeds per configuration (34 total predictions)
- Added statistical testing: Welch's t-test, Mann-Whitney U, Cohen's d, Bonferroni correction
- Key result: pLDDT@429 WT vs compound dimer p=0.0004 (Bonferroni adjusted p=0.003)
- RMSD validation against PDB 6FCX for all variants (all < 2.0 A)
- Added confidence interval and RMSD comparison charts
- Full per-residue pLDDT comparison plot
- Consolidated results into single alphafold/results_all/ folder
- Clean naming convention: {variant}_{type}_run{N}
- Added validate_structure.py script

## v1.0 (March 25, 2026)
- Initial 12-job AlphaFold 3 analysis (2 seeds per configuration)
- 4 Boltz-2 substrate/inhibitor binding predictions (THF, SAM)
- Cross-platform validation (AlphaFold 3 + Boltz-2)
- Core finding: compound heterozygous dimer shows lowest confidence values
- Automated analysis pipeline (analyze.py)
- Publication figure generator (generate_figures.py)
- HTML report generator
- Google Colab notebook
- Language compliance scanner (21 rules)
- 20 verified references
- Research paper draft (submission-ready)
- Master document v6.0
- ORCID: 0009-0000-1408-1065
- bioRxiv submission: BIORXIV/2026/715059

## v0.1 (March 25, 2026)
- Initial project setup
- Protein sequences verified against UniProt P42898
- AlphaFold Server JSON files created
- README and documentation framework
