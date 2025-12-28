#!/bin/bash
# 03_qc.sh - Quality control with FastQC

set -euo pipefail

echo "Running quality control..."

mkdir -p results/qc/raw

# Run FastQC
fastqc data/raw/*.fastq.gz \
    --outdir results/qc/raw \
    --threads 8

# Summarize with MultiQC
multiqc results/qc/raw \
    --outdir results/qc \
    --filename multiqc_raw \
    --force

echo "QC complete. Check: results/qc/multiqc_raw.html"
