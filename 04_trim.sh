#!/bin/bash
# 04_trim.sh - Adapter trimming

set -euo pipefail

echo "Trimming adapters..."

mkdir -p results/trimmed

for r1 in data/raw/*_1.fastq.gz; do
    sample=$(basename $r1 _1.fastq.gz)
    r2="data/raw/${sample}_2.fastq.gz"
    
    echo "Processing $sample..."
    
    trim_galore \
        --paired \
        --quality 20 \
        --length 36 \
        --cores 4 \
        --output_dir results/trimmed \
        $r1 $r2
done

echo "Trimming complete"
