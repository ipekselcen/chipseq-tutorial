#!/bin/bash
# 05_align.sh - Align reads with Bowtie2

set -euo pipefail

echo "Aligning reads..."

mkdir -p results/aligned
GENOME="data/reference/mm10"

for r1 in results/trimmed/*_1_val_1.fq.gz; do
    sample=$(basename $r1 _1_val_1.fq.gz)
    r2="results/trimmed/${sample}_2_val_2.fq.gz"
    out="results/aligned/${sample}"
    
    echo "Aligning $sample..."
    
    bowtie2 \
        -x $GENOME \
        -1 $r1 \
        -2 $r2 \
        --threads 8 \
        --very-sensitive \
        --no-discordant \
        --no-mixed \
        2> ${out}.log \
        | samtools view -@ 4 -bS - \
        | samtools sort -@ 4 -o ${out}.bam -
    
    samtools index ${out}.bam
    samtools flagstat ${out}.bam > ${out}_stats.txt
done

echo "Alignment complete"
