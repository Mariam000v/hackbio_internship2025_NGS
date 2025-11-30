#!/bin/bash

# make output directories if not exist
mkdir -p trimmed fastp_reports

# loop through all paired-end files
for file in *_1.fastq.gz
do
    # get the sample name (remove _1.fastq.gz)
    sample=$(basename "$file" _1.fastq.gz)

    echo "🔹 Processing $sample ..."

    fastp \
        -i ${sample}_1.fastq.gz \
        -I ${sample}_2.fastq.gz \
        -o trimmed/${sample}_1.trimmed.fastq.gz \
        -O trimmed/${sample}_2.trimmed.fastq.gz \
        -q 20 -u 30 -n 10 -l 50 \
        -h fastp_reports/${sample}_fastp.html \
        -j fastp_reports/${sample}_fastp.json \
        -w 8
done

echo "✅ All trimming finished. Results are in 'trimmed/' and reports in 'fastp_reports/'"

