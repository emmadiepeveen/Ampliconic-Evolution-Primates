# Gene Conversion Analysis

This directory contains the script for detecting and analyzing gene conversion events in ampliconic gene families.

## Scripts Overview

### `gene_conversion_tests.ipynb`
Analyzes gene conversion patterns by examining sequence similarity and GC content in ampliconic gene families.

**Input:**
- `gene_details` files from `scripts/01_ampliconic_detection/ampliconic_clustering_X/Y.ipynb` (in `data/intermediate/clusters/`)
- Reference genome assemblies (same as used in multicopy gene detection)
- Codon alignments from `scripts/01_ampliconic_detection/clustering_X.ipynb` and `clustering_Y.ipynb`

**Analysis steps:**
1. **Full coding sequence extraction:** Extracts complete coding sequences for all genes within ampliconic clusters using the gene_details files and reference genomes
2. **Pairwise sequence alignment:** Aligns all full-length genes within each cluster and calculates pairwise sequence differences
3. **GC3 content calculation:** Computes GC content at third codon positions (GC3) for each gene using the codon alignments
4. **Statistical testing:** Evaluates the significance of relationships between sequence similarity, copy number and distance using permutation tests
5. **Visualization:** Generates gene conversion plots for Figure 3 of the manuscript

**Output:**
- Gene conversion analysis results
- GC3 content per gene
- Pairwise sequence difference within ampliconic cluster per species 
- Figure 3 - gene conversion plots  
- Supplementary Table 6 (gene conversion summary statistics) 


## Workflow Summary
```
gene_details + Reference genomes
          ↓
Extract full coding sequences
          ↓
     Align genes
          ↓
Calculate pairwise differences

Codon alignments (from clustering)
          ↓
   Calculate GC3 content
          ↓
Permutation tests for significance
          ↓
  Gene conversion plots (Figure 3)
```

