# Selection Analysis

This directory contains scripts for analyzing selection in ampliconic gene families using PAML and visualises the selection patterns.

## Scripts Overview

### 1. `PAML_site_tests.ipynb`
Performs site-specific positive selection tests using PAML's codeml program.

**Input:**
- Multiple sequence alignments from `scripts/01_ampliconic_detection/clustering_ampliconic_X.ipynb` and `clustering_ampliconic_Y.ipynb`

**Output:**
- PAML site test results for each ampliconic cluster 
- Summary table indicating which gene families have positively selected sites and their omega (dN/dS) values across the phylogeny

**Note on execution:**
This notebook is computationally intensive. We recommend executing it via command line rather than running it interactively:
```bash
jupyter nbconvert --to notebook --execute PAML_site_tests.ipynb --output executed_PAML_site_tests.ipynb
```


### 2. `site_test_domain_mapping.ipynb`
Maps positively selected sites onto protein domain structures.

**Input:**
- Summary table of site test results from `PAML_site_tests.ipynb`
- Protein domain annotation files (TSV format) from `data/reference/protein_domains/`
- Specific sites of interest (identified from results from `PAML_site_tests.ipynb`)

**Output:**
- Visualization of positively selected sites mapped to human protein domain structures
- Generates Figure 5 and Supplementary Figure 7 of the manuscript

**Usage:**
Specify the sites of interest from the PAML summary table and the corresponding protein domain file. The script will create a domain diagram with positively selected sites highlighted.

### 3. `selection_heatmap_phylogenies.R`
Visualizes selection patterns across species and phylogenies.

**Input:**
- Bootstrapped pairwise dN/dS table from `scripts/01_ampliconic_detection/bootstrap_dnds.py`. Files are in `data/intermediate/dnds/`
- Branch test phylogenies (embedded in the script)

**Output:**
- Heatmaps showing pairwise dN/dS ratios across all species comparisons
- Phylogenetic trees with branch-specific dN/dS values from PAML branch tests
- Figures saved to `results/figures/`

**Note:**
The phylogenies from branch tests are hard-coded within the script with their dN/dS values already annotated.

## Workflow Summary
```
Alignments (from 01_ampliconic_detection)
          ↓
    PAML_site_tests.ipynb
          ↓
  Site test summary table
          ↓
site_test_domain_mapping.R  +  Protein domain files
          ↓
    Domain plots (Fig 5, Supp Fig 7)

Bootstrap dN/dS table (from 01_ampliconic_detection)
          ↓
selection_heatmap_phylogenies.R
          ↓
   Heatmaps + phylogeny plots (Fig 4)



