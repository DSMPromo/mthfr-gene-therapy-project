# Contributing to the MTHFR Gene Therapy Platform Project

Thank you for your interest in contributing! This is an open research project and we welcome collaborators from all backgrounds.

## How to Contribute

### 1. Replicate AlphaFold Results
The most immediately valuable contribution is independent replication of our structural predictions.

**Steps:**
1. Go to [alphafoldserver.com](https://alphafoldserver.com) (free, requires Google account)
2. Use the FASTA sequences in `sequences/` folder
3. Follow the submission plan in `alphafold/jobs/submission_plan.md`
4. Record your confidence metrics in the template at `analysis/metrics_template.csv`
5. Submit a pull request with your results

### 2. Validate with Experimental Data
If you have access to wet-lab facilities:
- Express wild-type and variant MTHFR and measure FAD binding affinity
- Compare thermal stability (the C677T variant is known to be thermolabile)
- Test base editor guide RNAs targeting the C677T and A1298C sites
- Retinal cell culture models with MTHFR knockdown/correction

### 3. Extend the Analysis
- Run molecular dynamics simulations on AlphaFold-predicted structures
- Perform binding energy calculations for FAD in WT vs variant contexts
- Analyze evolutionary conservation of positions 222 and 429
- Compare AlphaFold predictions against the experimental PDB 6FCX structure

### 4. Expand the Literature Review
- Add new publications connecting MTHFR to any of the five disease targets
- Identify additional MTHFR variants worth modeling
- Find clinical trial data on methylfolate/folinic acid interventions

### 5. Improve Documentation
- Fix errors, improve clarity, add diagrams
- Translate key documents into other languages
- Create video tutorials for the AlphaFold submission process

## What We Especially Need

| Role | Why |
|------|-----|
| Molecular biologist (CRISPR/base editing) | Design and test guide RNAs for C677T/A1298C correction |
| Ophthalmology researcher | Validate retinal degeneration connection, access to Mthfr mouse models |
| Autism/neurodevelopment researcher | Validate cerebral folate deficiency mechanism |
| Structural biologist | Compare AlphaFold predictions to experimental structures |
| Bioinformatician | Large-scale analysis of MTHFR variants across populations |
| Patent attorney | IP strategy for MTHFR base editing therapeutic claims |
| Grant writer | NIH/NSF/foundation funding applications |

## Pull Request Process

1. Fork the repository
2. Create a branch (`feature/your-contribution`)
3. Add your work with clear documentation
4. Submit a pull request describing what you did and what it means
5. We'll review and merge

## Code of Conduct

- All contributions must be for non-commercial research/educational purposes
- Be respectful and constructive
- Cite sources for all scientific claims
- Never overstate the clinical implications of computational predictions
- Remember: this is for the people — the kids with autism, the adults with unexplained anxiety, the millions losing their vision. Keep that mission front and center.

## Questions?

Open an issue or contact igor@dsm.promo
