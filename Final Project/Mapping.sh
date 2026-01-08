#!/bin/bash

# Directories
TRIM_DIR=~/proj_lncRNA_immuno/data/trimmed_reads
MAP_DIR=~/proj_lncRNA_immuno/data/mapping
LOG_DIR=~/proj_lncRNA_immuno/results/star_logs
INDEX=~/proj_lncRNA_immuno/refdata/hg38/STAR_index

# Number of parallel jobs
PARALLEL_JOBS=2
THREADS=8

# Counter for parallel jobs
count=0

# Loop through trimmed FASTQ files
for fq1 in $TRIM_DIR/*_1_trimmed.fastq.gz; do
    fq2=${fq1/_1_trimmed.fastq.gz/_2_trimmed.fastq.gz}
    sample=$(basename $fq1 | sed 's/_1_trimmed.fastq.gz//')

    echo "Mapping sample: $sample"

    # Run STAR in background and write all logs to LOG_DIR
    STAR --runThreadN $THREADS \
         --genomeDir $INDEX \
         --readFilesIn $fq1 $fq2 \
         --readFilesCommand zcat \
         --outFileNamePrefix $MAP_DIR/${sample}_ \
         --outSAMtype BAM SortedByCoordinate \
         --outSAMattributes All \
         --outStd Log > $LOG_DIR/${sample}_STAR.log 2>&1 &

    # Increment counter
    count=$((count+1))

    # If we reach PARALLEL_JOBS, wait for them to finish
    if [ $count -ge $PARALLEL_JOBS ]; then
        wait
        count=0
    fi

done

# Wait for any remaining jobs
wait

echo "All samples mapped successfully."
