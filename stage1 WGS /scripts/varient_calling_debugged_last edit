#!/bin/bash

# Set project base dir
proj_dir=~/Mariam/project1

# -------------------------
# Add Read Groups
# -------------------------
if [ ! -f $proj_dir/alignment_map/child_rg.bam ]; then
    gatk AddOrReplaceReadGroups \
        -I $proj_dir/alignment_map/child.bam \
        -O $proj_dir/alignment_map/child_rg.bam \
        -ID child -LB lib1 -PL ILLUMINA -PU unit1 -SM child
fi

if [ ! -f $proj_dir/alignment_map/father_rg.bam ]; then
    gatk AddOrReplaceReadGroups \
        -I $proj_dir/alignment_map/father.bam \
        -O $proj_dir/alignment_map/father_rg.bam \
        -ID father -LB lib1 -PL ILLUMINA -PU unit1 -SM father
fi

# -------------------------
# Mark Duplicates
# -------------------------
if [ ! -f $proj_dir/alignment_map/child_dedup.bam ]; then
    gatk MarkDuplicates \
        -I $proj_dir/alignment_map/child_rg.bam \
        -O $proj_dir/alignment_map/child_dedup.bam \
        -M $proj_dir/alignment_map/child_dedup.metrics
fi

if [ ! -f $proj_dir/alignment_map/father_dedup.bam ]; then
    gatk MarkDuplicates \
        -I $proj_dir/alignment_map/father_rg.bam \
        -O $proj_dir/alignment_map/father_dedup.bam \
        -M $proj_dir/alignment_map/father_dedup.metrics
fi

# -------------------------
# Build BAM Index
# -------------------------
gatk BuildBamIndex -I $proj_dir/alignment_map/child_dedup.bam
gatk BuildBamIndex -I $proj_dir/alignment_map/father_dedup.bam

# -------------------------
# Base Quality Score Recalibration (BQSR)
# -------------------------
mkdir -p $proj_dir/BQSR

# Child
gatk BaseRecalibrator \
  -R /data/Homo_sapiens_assembly38.fasta \
  -I $proj_dir/alignment_map/child_dedup.bam \
  --known-sites /data/Homo_sapiens_assembly38.known_indels.vcf.gz \
  -O $proj_dir/BQSR/child_recal_data.table

gatk ApplyBQSR \
  -R /data/Homo_sapiens_assembly38.fasta \
  -I $proj_dir/alignment_map/child_dedup.bam \
  --bqsr-recal-file $proj_dir/BQSR/child_recal_data.table \
  -O $proj_dir/BQSR/child_recal.bam

# Father
gatk BaseRecalibrator \
  -R /data/Homo_sapiens_assembly38.fasta \
  -I $proj_dir/alignment_map/father_dedup.bam \
  --known-sites /data/Homo_sapiens_assembly38.known_indels.vcf.gz \
  -O $proj_dir/BQSR/father_recal_data.table

gatk ApplyBQSR \
  -R /data/Homo_sapiens_assembly38.fasta \
  -I $proj_dir/alignment_map/father_dedup.bam \
  --bqsr-recal-file $proj_dir/BQSR/father_recal_data.table \
  -O $proj_dir/BQSR/father_recal.bam

# -------------------------
# HaplotypeCaller
# -------------------------
mkdir -p $proj_dir/VCF

# Child
gatk HaplotypeCaller \
  -R /data/Homo_sapiens_assembly38.fasta \
  -I $proj_dir/BQSR/child_recal.bam \
  -O $proj_dir/VCF/child.g.vcf.gz \
  -ERC GVCF

# Father
gatk HaplotypeCaller \
  -R /data/Homo_sapiens_assembly38.fasta \
  -I $proj_dir/BQSR/father_recal.bam \
  -O $proj_dir/VCF/father.g.vcf.gz \
  -ERC GVCF

# -------------------------
# Combine GVCFs
# -------------------------
gatk CombineGVCFs \
  -R /data/Homo_sapiens_assembly38.fasta \
  --variant $proj_dir/VCF/child.g.vcf.gz \
  --variant $proj_dir/VCF/father.g.vcf.gz \
  -O $proj_dir/VCF/family_combined.g.vcf.gz

# -------------------------
# Joint Genotyping
# -------------------------
gatk GenotypeGVCFs \
  -R /data/Homo_sapiens_assembly38.fasta \
  -V $proj_dir/VCF/family_combined.g.vcf.gz \
  -O $proj_dir/VCF/family_variants.vcf.gz
