#!/bin/bash
# 06_filter.sh - Filter and deduplicate BAMs

set -euo pipefail

echo "Filtering BAM files..."

BLACKLIST="data/reference/mm10-blacklist.v2.bed"

for bam in results/aligned/*.bam; do
    sample=$(basename $bam .bam)
    
    # Skip if already filtered
    if [[ $sample == *"_filtered" ]]; then
        continue
    fi
    
    echo "Processing $sample..."
    
    # Filter by quality and remove blacklist
    samtools view -@ 4 -b -q 20 -F 1804 $bam \
        | bedtools intersect -v -abam - -b $BLACKLIST \
        > results/aligned/${sample}_temp.bam
    
    # Remove duplicates
    picard MarkDuplicates \
        I=results/aligned/${sample}_temp.bam \
        O=results/aligned/${sample}_filtered.bam \
        M=results/aligned/${sample}_dup_metrics.txt \
        REMOVE_DUPLICATES=true \
        ASSUME_SORTED=true \
        VALIDATION_STRINGENCY=LENIENT
    
    samtools index results/aligned/${sample}_filtered.bam
    rm results/aligned/${sample}_temp.bam
done

echo "Filtering complete"
