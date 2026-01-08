#This is a script to redo mapping for some failed samples.

#!/bin/bash
# Remap failed RNA-Seq samples using STAR with trimmed reads
# Run up to 2 samples in parallel; stop immediately if any fails

set -euo pipefail

# Number of threads
THREADS=4

# Paths
GENOME_DIR=~/proj_lncRNA_immuno/refdata/hg38/STAR_index
TRIMMED_DIR=~/proj_lncRNA_immuno/data/trimmed_reads
OUTPUT_DIR=~/proj_lncRNA_immuno/data/mapping
LOG_DIR=~/proj_lncRNA_immuno/results/star_logs

# Function to run STAR for one sample
run_star() {
    SAMPLE=$1
    echo "=============================="
    echo "Processing sample: $SAMPLE"

    READ1="$TRIMMED_DIR/${SAMPLE}_1_trimmed.fastq.gz"
    READ2="$TRIMMED_DIR/${SAMPLE}_2_trimmed.fastq.gz"
    BAM_FILE="$OUTPUT_DIR/${SAMPLE}_Aligned.sortedByCoord.out.bam"
    LOG_FILE="$LOG_DIR/${SAMPLE}_STAR.log"
    TMP_DIR="$LOG_DIR/${SAMPLE}_tmp"

    # Skip if BAM already exists
    if [[ -f "$BAM_FILE" ]]; then
        echo "⚠️ BAM already exists for $SAMPLE. Skipping..."
        return
    fi

    # Check if trimmed reads exist
    if [[ ! -f "$READ1" || ! -f "$READ2" ]]; then
        echo "❌ Warning: trimmed reads missing for $SAMPLE. Skipping..."
        return
    fi

    # Create temporary directory for STAR
    mkdir -p "$TMP_DIR"

    # Run STAR
    if STAR --runThreadN "$THREADS" \
            --genomeDir "$GENOME_DIR" \
            --readFilesIn "$READ1" "$READ2" \
            --readFilesCommand zcat \
            --outFileNamePrefix "$OUTPUT_DIR/${SAMPLE}_" \
            --outSAMtype BAM SortedByCoordinate \
            --outSAMattributes All \
            --outTmpDir "$TMP_DIR" > "$LOG_FILE" 2>&1; then
        echo "✅ Done mapping: $SAMPLE"
        rm -rf "$TMP_DIR"  # cleanup temp dir after success
    else
        echo "❌ STAR failed for $SAMPLE. Removing partial files and logs..."
        rm -f "$OUTPUT_DIR/${SAMPLE}"_*
        rm -rf "$TMP_DIR"
        rm -f "$LOG_DIR/${SAMPLE}"_*
        echo "🛑 STAR failed. Exiting."
        exit 1
    fi
}

# Run samples (max 2 in parallel)
N=0
while read -r SAMPLE; do
    run_star "$SAMPLE" &
    ((N++))
    if [[ $N -ge 2 ]]; then
        wait
        N=0
    fi
done < ~/proj_lncRNA_immuno/scripts/failed_samples.txt
wait

echo "✅ All failed samples processed successfully."
