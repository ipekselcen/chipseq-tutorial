#!/bin/bash
# 02_download.sh - Download ChIP-seq data from SRA

set -euo pipefail

echo "Downloading data from GEO (GSE263808)..."

mkdir -p data/raw
cd data/raw

# Sample mapping (replace with actual SRR numbers from GEO)
# Visit: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE263808

samples=(
    "nOPC_rep1_ChIP:SRR_PLACEHOLDER"
    "nOPC_rep1_Input:SRR_PLACEHOLDER"
    "nOPC_rep2_ChIP:SRR_PLACEHOLDER"
    "nOPC_rep2_Input:SRR_PLACEHOLDER"
    "nOPC_rep3_ChIP:SRR_PLACEHOLDER"
    "nOPC_rep3_Input:SRR_PLACEHOLDER"
    "aOPC_rep1_ChIP:SRR_PLACEHOLDER"
    "aOPC_rep1_Input:SRR_PLACEHOLDER"
    "aOPC_rep2_ChIP:SRR_PLACEHOLDER"
    "aOPC_rep2_Input:SRR_PLACEHOLDER"
    "aOPC_rep3_ChIP:SRR_PLACEHOLDER"
    "aOPC_rep3_Input:SRR_PLACEHOLDER"
)

for entry in "${samples[@]}"; do
    name="${entry%%:*}"
    srr="${entry##*:}"
    
    if [[ $srr == "SRR_PLACEHOLDER" ]]; then
        echo "Skipping $name - add real SRR number"
        continue
    fi
    
    echo "Downloading $name ($srr)..."
    
    # Download with prefetch and convert to FASTQ
    prefetch $srr
    fasterq-dump $srr --split-files --threads 8
    
    # Rename and compress
    mv ${srr}_1.fastq ${name}_1.fastq
    mv ${srr}_2.fastq ${name}_2.fastq
    pigz -p 4 ${name}_1.fastq &
    pigz -p 4 ${name}_2.fastq &
    wait
done

cd ../..
echo "Download complete"
