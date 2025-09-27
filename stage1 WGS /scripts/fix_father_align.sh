#!/bin/bash

# Script: fix_father_alignment.sh
# Purpose: Repair and realign father's WGS reads due to errors in previous alignment script.
# Notes:
#   - Previous alignment produced incomplete or corrupted father.bam
#   - Possible causes: interrupted BWA alignment, permission issues, or file corruption
#   - This script regenerates father.bam from father.fastq, sorts, and indexes it
# Tools: BWA, Samtools

# --- Paths ---
REF="/home/yossra_bioinformatics/Mariam/project1/Homo_sapiens_assembly38.fasta"  # your local copy
READS_DIR="/home/yossra_bioinformatics/Mariam/project1/qc_report/fastp_reports/repaired_reads"  # trimmed & repaired reads
OUT_DIR="/home/yossra_bioinformatics/Mariam/project1/alignment_map"

# --- Input files ---
R1="$READS_DIR/father_R1.repaired.fastq.gz"
R2="$READS_DIR/father_R2.repaired.fastq.gz"

# --- Output files ---
SAM="$OUT_DIR/father.sam"
BAM="$OUT_DIR/father.bam"
SORTED_BAM="$OUT_DIR/father_sorted.bam"

echo ">>> Starting father alignment repair..."

# --- Step 1: Align with BWA-MEM ---
echo ">>> Aligning father reads..."
bwa mem -t 4 $REF $R1 $R2 > $SAM
if [ $? -ne 0 ]; then
    echo "!!! BWA alignment failed. Check input reads and reference."
    exit 1
fi

# --- Step 2: Convert SAM to BAM ---
echo ">>> Converting SAM to BAM..."
samtools view -b $SAM > $BAM
if [ $? -ne 0 ]; then
    echo "!!! SAM to BAM conversion failed."
    exit 1
fi

# --- Step 3: Sort BAM ---
echo ">>> Sorting BAM..."
samtools sort -@4 -o $SORTED_BAM $BAM
if [ $? -ne 0 ]; then
    echo "!!! BAM sorting failed."
    exit 1
fi

# --- Step 4: Index BAM ---
echo ">>> Indexing sorted BAM..."
samtools index $SORTED_BAM
if [ $? -ne 0 ]; then
    echo "!!! BAM indexing failed."
    exit 1
fi

echo ">>> Father alignment repair complete!"
echo "Sorted and indexed BAM: $SORTED_BAM"

