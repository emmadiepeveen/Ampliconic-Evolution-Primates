**Goal**
This script is to find gene copies of ampliconic gene families in long-read de novo assemblies. 


**Input:**
- Reference genome assemblies (same as used in multicopy gene detection)
- Long read de novo assemblies from Chimpanzees, Humans and gorillas (Long reads used are from Porsborg, Peter Soerud, Anders Poulsen Charmouh, Vinod Kumar Singh, et al. 2025. ‘Long-Read Sequencing of Primate Testis and Human Sperm Allows Identification of Recombination Events in Individuals’. Nature Communications 16 (1): 10337. https://doi.org/10.1038/s41467-025-65248-3.)

**Analysis steps:**
1. **Identify all X-linked contigs:** Use only assembled contigs that have X chromosome content 
2. **Pull out reference queries of the reference T2T genome:** get all coding sequences of all loci in an ampliconic cluster from the reference genome
3. **BLAT reference mRNAs against each individual's X-contigs:** Find all possible hits of the ampliconic copies
4. **Filter for chimeric alignments:** exclude all alignments that span over multiple copies 
5. **Cluster hits into discrete loci:** Overlapping hits from different loci are merged into single copies
6. **Pull each possible copy interval:** Make a separate FASTA file for every possible copy
7. **Align the reference protein to each copy to get its CDS:** Align protein of reference to the possible copies' sequence and find its CDS
8. **Copy number estimation for copies with ORF and no frameshifts:** Only keep copies with ORF and no frameshits 
9. **Calculate pairwise pN/pS:** From each ampliconic cluster per individual
10. **Build gene trees per family per species:** inferred a maximum-likelihood tree for each family and species from the codon alignment






**Output:**
- Ampliconic copy number estimation (Figure 6)
- pN/pS for each ampliconic family per individual (Supplementary Table 12)
- gene trees per family per species (Supplementary Figure 8-11)


## Workflow Summary
```
Identify all X-linked contigs (run wfmash.sh)
          ↓
Pull out reference queries of the reference T2T genome
          ↓
BLAT reference mRNAs against each individual's X-contigs
          ↓
Filter for chimeric alignments

Cluster hits into discrete loci
          ↓
Pull each possible copy interval 
          ↓
Align the reference protein to each copy to get its CDS
          ↓
Copy number estimation for copies with ORF and no frameshifts
          ↓
Calculate pairwise pN/pS from each ampliconic cluster per individual
          ↓
Build gene trees per family per species
```
