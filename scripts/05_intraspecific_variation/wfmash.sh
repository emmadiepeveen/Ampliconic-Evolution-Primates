#!/bin/bash
#SBATCH --job-name=wfmash
#SBATCH --account 
#SBATCH --output=logs/wfmash_%j.out
#SBATCH --error=logs/wfmash_%j.err
#SBATCH --time=03:00:00
#SBATCH --nodes=4
#SBATCH --ntasks=4
#SBATCH --array=0-5
#SBATCH --cpus-per-task=4
#SBATCH --mem=40G

set -euo pipefail
date
echo "job $SLURM_JOB_ID  task $SLURM_ARRAY_TASK_ID  on $(hostname)"

# ===== species config — edit these 4 lines to switch species =================
INDIVIDUALS=(indi1 indi2 indi3)
SPECIES_DIR=HomSap # change to species 
SPECIES_TAG=human
REF_BUILD=GCF_009914755.1
# ===============================================================================

REF_X=/path/to/genome/X_chromosome/directory/${SPECIES_DIR}_X.fasta
REF_Y=/path/to/genome/Y_chromosome/directory/${SPECIES_DIR}_Y.fasta
GENOME_TMPL=path/to/haploid/denovo/assembly/IND_haploid.fa
OUTDIR=/pop_data/$SPECIES_TAG

BLOCK_MIN=20000      # wfmash -l : drop mapping blocks shorter than this
MIN_X_ALN=20000      # keep a contig only if its TOTAL X alignment >= this
MIN_X_FRAC=0.30      #  >= this fraction of the contig maps to X
MAP_PID=95           # wfmash -p : minimum mapping identity (%)
THREADS=${SLURM_CPUS_PER_TASK:-10}

# ----- pick this task's individual --------------------------------------------
IND=${INDIVIDUALS[$SLURM_ARRAY_TASK_ID]}
GENOME=${GENOME_TMPL//IND/$IND}
ODIR=$OUTDIR/$IND
mkdir -p "$ODIR"
echo "individual: $IND"
echo "genome:     $GENOME"

XPAF=$ODIR/${IND}_vs_X.paf
YPAF=$ODIR/${IND}_vs_Y.paf

# ----- map assembly (query) to X and Y (target), approximate mapping only -----
/wfmash -m -p "$MAP_PID" -l "$BLOCK_MIN" -t "$THREADS" "$REF_X" "$GENOME" > "$XPAF"
wfmash -m -p "$MAP_PID" -l "$BLOCK_MIN" -t "$THREADS" "$REF_Y" "$GENOME" > "$YPAF"

# ----- keep X-dominant contigs -------------------------------------------------
# PAF col 1 = query contig, col 2 = contig length, col 11 = block length.
# Sum block length per contig for X and Y; keep a contig if its X total clears
# MIN_X_ALN, is >= its Y total, and covers >= MIN_X_FRAC of the contig.
KEEP=$ODIR/${IND}_Xcontigs.txt
awk -v MINX="$MIN_X_ALN" -v FRAC="$MIN_X_FRAC" '
  NR==FNR { x[$1] += $11; qlen[$1] = $2; next }   # first file = X paf
  { y[$1] += $11 }                                # second file = Y paf
  END { for (c in x)
          if (x[c] >= MINX && x[c] >= y[c] && x[c] >= FRAC * qlen[c]) print c }
' "$XPAF" "$YPAF" > "$KEEP"
echo "X-linked contigs kept: $(wc -l < "$KEEP")"

# ----- subset the assembly to those contigs (pure awk, no extra tools) --------
XFA=$ODIR/${IND}_Xcontigs.fa
awk 'NR==FNR { keep[$1]=1; next }
     /^>/    { name=substr($1,2); p=(name in keep) }
     p' "$KEEP" "$GENOME" > "$XFA"

# ----- index for later stages --------------------------------------------------
if command -v samtools >/dev/null 2>&1; then
    samtools faidx "$XFA"
    echo "kept total: $(awk '{s+=$2} END{printf "%.0f Mb\n", s/1e6}' "$XFA.fai")"
fi

echo "done -> $XFA"
date
