#!/bin/bash

# Directories
TRIM_DIR="trim/trimmed"
STAR_INDEX="ref/star_index"
ALIGN_DIR="alignment_star"
LOG_DIR="alignment_logs"

# Make output directories
mkdir -p $ALIGN_DIR $LOG_DIR

# Maximum number of parallel STAR processes
MAX_JOBS=4
job_count=0

# Loop through all _1.trimmed.fastq.gz files
for f1 in $TRIM_DIR/*_1.trimmed.fastq.gz
do
    sample=$(basename "$f1" _1.trimmed.fastq.gz)
    f2="$TRIM_DIR/${sample}_2.trimmed.fastq.gz"

    echo "🔹 Aligning $sample ..."

    STAR \
        --runThreadN 8 \
        --genomeDir $STAR_INDEX \
        --readFilesIn $f1 $f2 \
        --readFilesCommand zcat \
        --outFileNamePrefix $ALIGN_DIR/${sample}_ \
        --outSAMtype BAM SortedByCoordinate \
        > $LOG_DIR/${sample}.log 2>&1 &

    job_count=$((job_count+1))

    # Limit number of parallel jobs
    if [ "$job_count" -ge "$MAX_JOBS" ]; then
        wait
        job_count=0
    fi
done

# Wait for any remaining background jobs
wait

echo "🎉 All alignments finished. BAM files are in $ALIGN_DIR"

