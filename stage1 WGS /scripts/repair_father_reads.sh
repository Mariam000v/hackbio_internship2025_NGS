#!/bin/bash

# Script: father_repair.sh
# Purpose: Repair father's paired-end FASTQ reads using BBTools repair.sh
# Notes:
#   - This ensures read pairs are synchronized (no missing mates).
#   - Outputs repaired paired-end files and singletons.

# --- Paths ---
READS_DIR="/home/yossra_bioinformatics/Mariam/project1/qc_report/fastp_reports"
REPAIR_DIR="$READS_DIR/repaired_reads"

# --- Input files ---
R1="$READS_DIR/father_R1.trim.fastq.gz"
R2="$READS_DIR/father_R2.trim.fastq.gz"

# --- Output files ---
R1_REP="$REPAIR_DIR/father_R1.repaired.fastq.gz"
R2_REP="$REPAIR_DIR/father_R2.repaired.fastq.gz"
SINGLE="$REPAIR_DIR/father_singletons.fastq.gz"

# --- Make repair directory if not exists ---
mkdir -p $REPAIR_DIR

echo ">>> Repairing father FASTQ reads..."
repair.sh in1=$R1 in2=$R2 out1=$R1_REP out2=$R2_REP outs=$SINGLE
echo ">>> Repair complete: $R1_REP, $R2_REP, $SINGLE"

