#!/bin/bash
# This script copies the human genome reference files from project1 to project2

# Create ref directory if not exist
mkdir -p ref

# Copy all fasta and index files
cp ../project1/Homo_sapiens_assembly38.fasta* ref/

echo "✅ Reference files copied to ref/"
