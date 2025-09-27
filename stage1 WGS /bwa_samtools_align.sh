#!/bin/bash

# Script: bwa_samtools_align.sh
# Purpose: Align trimmed reads to reference genome, then sort & index BAM
# Tools: BBTools (repair.sh), BWA, Samtools

# --- Paths ---
REF="/home/yossra_bioinformatics/Mariam/project1/Homo_sapiens_assembly38.fasta"
READS_DIR="/home/yossra_bioinformatics/Mariam/project1/qc_report/fastp_reports"
OUT_DIR="/home/yossra_bioinformatics/Mariam/project1/alignment_map"
REPAIR_DIR="$READS_DIR/repaired_reads"

# --- Create output directories ---
mkdir -p $OUT_DIR
mkdir -p $REPAIR_DIR

# --- Step 1: Index reference genome ---
if [ ! -f "${REF}.bwt" ]; then
    echo ">>> Indexing reference genome with BWA ..."
    bwa index $REF
    echo ">>> Reference indexing complete."
else
    echo ">>> Reference already indexed. Skipping."
fi

# --- Samples to process ---
samples=("child" "father")

# --- Step 2: Repair with BBTools ---
for sample in "${samples[@]}"; do
    R1="$READS_DIR/${sample}_R1.trim.fastq.gz"
    R2="$READS_DIR/${sample}_R2.trim.fastq.gz"
    R1_REP="$REPAIR_DIR/${sample}_R1.repaired.fastq.gz"
    R2_REP="$REPAIR_DIR/${sample}_R2.repaired.fastq.gz"

    if [ ! -f "$R1_REP" ] || [ ! -f "$R2_REP" ]; then
        echo ">>> Repairing reads for $sample ..."
        repair.sh in1=$R1 in2=$R2 out1=$R1_REP out2=$R2_REP outs=$REPAIR_DIR/${sample}_singletons.fastq.gz
        echo ">>> Repair complete: $R1_REP , $R2_REP"
    else
        echo ">>> Repaired files for $sample already exist. Skipping."
    fi
done

# --- Step 3: Alignment loop ---
for sample in "${samples[@]}"; do
    echo ">>> Aligning $sample ..."

    # Input repaired files
    R1="$REPAIR_DIR/${sample}_R1.repaired.fastq.gz"
    R2="$REPAIR_DIR/${sample}_R2.repaired.fastq.gz"

    # Output files
    SAM="$OUT_DIR/${sample}.sam"
    BAM="$OUT_DIR/${sample}.bam"
    SORTED_BAM="$OUT_DIR/${sample}_sorted.bam"

    # Step 1: Align with BWA-MEM
    bwa mem -t 4 $REF $R1 $R2 > $SAM

    # Step 2: Convert SAM to BAM
    samtools view -b $SAM > $BAM

    # Step 3: Sort BAM
    samtools sort -o $SORTED_BAM $BAM

    # Step 4: Index BAM
    samtools index $SORTED_BAM

    echo ">>> $sample alignment complete: $SORTED_BAM"
done

echo "All alignments finished."

