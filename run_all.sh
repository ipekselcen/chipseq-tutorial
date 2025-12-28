#!/bin/bash
# run_all.sh - Execute complete ChIP-seq pipeline

set -euo pipefail

echo "Starting ChIP-seq analysis pipeline"
echo "===================================="

# Setup
echo "[1/10] Setting up reference genome..."
bash scripts/01_setup.sh

# Download
echo "[2/10] Downloading data..."
bash scripts/02_download.sh

# QC
echo "[3/10] Quality control..."
bash scripts/03_qc.sh

# Trim
echo "[4/10] Trimming adapters..."
bash scripts/04_trim.sh

# Align
echo "[5/10] Aligning reads..."
bash scripts/05_align.sh

# Filter
echo "[6/10] Filtering BAMs..."
bash scripts/06_filter.sh

# Peaks
echo "[7/10] Calling peaks..."
bash scripts/07_peaks.sh

# Differential
echo "[8/10] Differential binding..."
bash scripts/08_differential.sh

# Annotate
echo "[9/10] Annotating peaks..."
Rscript scripts/09_annotate.R

# Visualize
echo "[10/10] Creating visualizations..."
Rscript scripts/10_visualize.R

echo ""
echo "===================================="
echo "Pipeline complete!"
echo "===================================="
echo "Check results in results/ directory"
