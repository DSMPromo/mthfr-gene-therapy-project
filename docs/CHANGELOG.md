# Changelog

## v1.6 (April 3, 2026) -- COMPLETE
- Retina-only focus: removed neuropsychiatric secondary exploratory context from all documents
- Retinal biomarkers now sole primary validation pathway across all documents
- Added "Why Retinal Focus?" section to README with 4 key rationale points
- Retinal pathway prioritized for: quantifiable OCT biomarkers, established Mthfr+/- mouse models, genotype-response correlation in human studies, case evidence of visual recovery
- Roadmap Phase 4 updated: retinal OCT (nerve fiber layer, macular perfusion) as primary endpoints
- Collaboration section focused on ophthalmologists/retinal researchers
- Neuropsychiatric, cardiovascular, reproductive, epigenetic domains retained as background literature context only
- All 4 documents + PDF updated consistently
- Research Paper: Discussion Section 4.3 rewritten as "Retinal Follow-Up as the Primary Validation Context"

## v1.5 (April 3, 2026) -- COMPLETE
- Full document rebuild: all 4 documents + Colab notebook updated and cross-verified
- Critical data correction: all averaged values recomputed as true n=10 from metrics.csv
- Fixed 8 orphaned references in Research Paper (Wan, McNulty, Karahan, Pentieva, Xia, Kumar, YOLT-101, LNP base editing) -- all now cited in body text
- Fixed Appendix B paths: json/ → json_all/, results/ → results_all/
- Master Document updated to v7.0 with Boltz-2 results table and full 100ns MD results section
- Added 3 missing references to Master Document (Weisberg 1998, Liew & Gupta 2015, Levine 2018)
- README: updated "Where we are now" from Phase 1 to Phase 2a complete
- README: added verify_md.py and analysis/md_results/ to project structure tree
- README: fixed Master Document version reference (v6.0 → v7.0)
- Supplementary statistics regenerated from true n=10 seed data
- Boltz-2 results table added to Research Paper Section 3.7 (was empty placeholder)
- All documents verified for cross-consistency of averaged values, phase status, and version numbers

## v1.4 (April 3, 2026) -- COMPLETE
- Extended molecular dynamics to 100ns (was 10ns) on RTX 4090
- Fixed PBC artifact: per-chain RMSD analysis instead of whole dimer
- 100ns results: compound dimer more compact than WT (RMSD 6.22 vs 7.16 A)
- Equilibrium RMSD (>50ns): compound 6.88 vs WT 8.17 A (Cohen's d = 4.31)
- Per-chain RMSD: WT A=6.89/B=7.44, Compound A=5.67/B=6.76
- t-test (equilibrium): p < 1e-323 (highly significant)
- Both systems reach equilibrium at ~62-66ns
- Added verify_md.py: 34 independent verification checks (validated 2x)
- Added strided trajectory loading for memory-safe analysis of 50+ GB DCD files
- Added RMSD CSV exports for reproducibility
- Updated all documents with corrected 100ns findings
- Updated run_md.py with PBC correction and per-chain analysis

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
