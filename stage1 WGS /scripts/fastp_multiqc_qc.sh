#!/bin/bash
# fastp_multiqc_qc.sh
# Purpose: Run fastp QC + trimming for child and father, then summarize with MultiQC
#          Includes both raw and trimmed FASTQ in MultiQC
# Author: Mariam
# Date: 2025-09-17

# ------------------------------
# Settings
# ------------------------------
THREADS=4
OUTDIR=$(pwd)
FASTQ_DIR=/data/human_stage_1

# Input files
CHILD_R1=$FASTQ_DIR/child_1.fastq.gz
CHILD_R2=$FASTQ_DIR/child_2.fastq.gz
FATHER_R1=$FASTQ_DIR/father_1.fastq.gz
FATHER_R2=$FASTQ_DIR/father_2.fastq.gz

# Create output folders
mkdir -p $OUTDIR/fastp_reports
mkdir -p $OUTDIR/multiqc_report

# ------------------------------
# 1. Run fastp on child
# ------------------------------
echo "Running fastp on child..."
fastp \
  -i $CHILD_R1 \
  -I $CHILD_R2 \
  -o $OUTDIR/fastp_reports/child_R1.trim.fastq.gz \
  -O $OUTDIR/fastp_reports/child_R2.trim.fastq.gz \
  -h $OUTDIR/fastp_reports/child_fastp.html \
  -j $OUTDIR/fastp_reports/child_fastp.json \
  -w $THREADS

# ------------------------------
# 2. Run fastp on father
# ------------------------------
echo "Running fastp on father..."
fastp \
  -i $FATHER_R1 \
  -I $FATHER_R2 \
  -o $OUTDIR/fastp_reports/father_R1.trim.fastq.gz \
  -O $OUTDIR/fastp_reports/father_R2.trim.fastq.gz \
  -h $OUTDIR/fastp_reports/father_fastp.html \
  -j $OUTDIR/fastp_reports/father_fastp.json \
  -w $THREADS

# ------------------------------
# 3. Summarize with MultiQC including raw + trimmed FASTQ
# ------------------------------
echo "Running MultiQC to summarize all reports..."
multiqc $FASTQ_DIR $OUTDIR/fastp_reports -o $OUTDIR/multiqc_report

echo "QC completed."
echo "Fastp outputs: $OUTDIR/fastp_reports"
echo "MultiQC summary: $OUTDIR/multiqc_report/multiqc_report.html"

