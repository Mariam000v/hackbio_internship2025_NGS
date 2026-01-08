
# RNA-seq Mapping Workflow (hg38, STAR)

This document describes the RNA-seq alignment workflow used to map trimmed paired-end reads to the human reference genome (hg38) using STAR and GENCODE v49 annotation.  

---

## 1. Reference preparation

```bash
# Move to the reference data folder
cd ~/proj_lncRNA_immuno/refdata/hg38

# Download the GTF annotation
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz

# Unzip the annotation
gunzip gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz

# Check files
ls -lh

##2. STAR genome index generation
# Make a folder for the STAR index
mkdir -p ~/proj_lncRNA_immuno/refdata/hg38/STAR_index

# Run STAR genome index generation
STAR --runThreadN 8 \
     --runMode genomeGenerate \
     --genomeDir ~/proj_lncRNA_immuno/refdata/hg38/STAR_index \
     --genomeFastaFiles ~/proj_lncRNA_immuno/refdata/hg38/hg38.fasta \
     --sjdbGTFfile ~/proj_lncRNA_immuno/refdata/hg38/gencode.v49.chr_patch_hapl_scaff.annotation.gtf


##3. Read alignment
# Run STAR mapping for all samples
./scripts/mapping.sh

##4. Identification of failed alignments
# Find BAM files with 0 size and create a list of failed samples
find . -name "*.bam" -size 0 | sed 's|^\./||; s/_Aligned.*//' > scripts/failed_samples.txt

##5. Inspection of failed samples
cd ~/proj_lncRNA_immuno

# Check associated files for failed samples
while read sample; do
    find data/mapping results/star_logs -name "${sample}*"
done < scripts/failed_samples.txt

##6. Cleanup of failed outputs (with logging)
cd ~/proj_lncRNA_immuno

# Create a log file for deletions
logfile="scripts/deleted_failed_samples.log"
echo "Deletion log - $(date)" > $logfile

# Delete files for failed samples while keeping a log
while read sample; do
    echo "Checking files for: $sample"

    # Find matching files
    files=$(find data/mapping results/star_logs -name "${sample}*")
    
    if [ -n "$files" ]; then
        echo "Deleting the following files for $sample:" | tee -a $logfile
        echo "$files" | tee -a $logfile
        
        # Actually delete the files
        echo "$files" | xargs rm -rf
    else
        echo "No files found for $sample" | tee -a $logfile
    fi
done < scripts/failed_samples.txt

##7. Verification and file counts

# Verify that deleted files no longer exist
cd ~/proj_lncRNA_immuno

while read sample; do
    echo "Checking for: $sample"
    find data/mapping results/star_logs -name "${sample}*"
done < scripts/failed_samples.txt

# Count remaining files
echo "Files in data/mapping:"
find data/mapping -type f | wc -l

echo "Files in results/star_logs:"
find results/star_logs -type f | wc -l

echo "Number of log files:"
find results/star_logs -name "*.log" | wc -l

echo "Number of failed samples:"
wc -l scripts/failed_samples.txt


## 8. Re-alignment of failed samples
# Make the re-mapping script executable
chmod +x scripts/redo_mapping.sh

# Run re-mapping for failed samples
./scripts/redo_mapping.sh


