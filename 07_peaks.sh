#!/bin/bash
# 07_peaks.sh - Peak calling with MACS2

set -euo pipefail

echo "Calling peaks with MACS2..."

mkdir -p results/peaks

# Define ChIP-Input pairs
samples=(
    "nOPC_rep1:nOPC_rep1_ChIP:nOPC_rep1_Input"
    "nOPC_rep2:nOPC_rep2_ChIP:nOPC_rep2_Input"
    "nOPC_rep3:nOPC_rep3_ChIP:nOPC_rep3_Input"
    "aOPC_rep1:aOPC_rep1_ChIP:aOPC_rep1_Input"
    "aOPC_rep2:aOPC_rep2_ChIP:aOPC_rep2_Input"
    "aOPC_rep3:aOPC_rep3_ChIP:aOPC_rep3_Input"
)

for entry in "${samples[@]}"; do
    IFS=':' read -r name chip input <<< "$entry"
    
    chip_bam="results/aligned/${chip}_filtered.bam"
    input_bam="results/aligned/${input}_filtered.bam"
    
    if [ ! -f "$chip_bam" ]; then
        echo "Skipping $name - BAM not found"
        continue
    fi
    
    echo "Calling peaks for $name..."
    
    # Call broad peaks (for histone marks)
    macs2 callpeak \
        -t $chip_bam \
        -c $input_bam \
        -f BAMPE \
        -g mm \
        --broad \
        --broad-cutoff 0.01 \
        -n $name \
        --outdir results/peaks
done

# Convert bedGraph to bigWig
CHROMSIZES="data/reference/mm10.chrom.sizes"

for bdg in results/peaks/*_treat_pileup.bdg; do
    name=$(basename $bdg _treat_pileup.bdg)
    sort -k1,1 -k2,2n $bdg > results/peaks/${name}_sorted.bdg
    bedGraphToBigWig results/peaks/${name}_sorted.bdg $CHROMSIZES results/peaks/${name}.bw
    rm results/peaks/${name}_sorted.bdg
done

echo "Peak calling complete"
