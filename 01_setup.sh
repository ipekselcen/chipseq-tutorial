#!/bin/bash
# 01_setup.sh - Download and index mm10 genome

set -euo pipefail

echo "Setting up reference genome..."

REFDIR="data/reference"
mkdir -p $REFDIR
cd $REFDIR

# Download mm10 genome
if [ ! -f "mm10.fa" ]; then
    echo "Downloading mm10 genome..."
    wget -q http://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz
    gunzip mm10.fa.gz
fi

# Download chromosome sizes
if [ ! -f "mm10.chrom.sizes" ]; then
    echo "Downloading chromosome sizes..."
    wget -q http://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes
fi

# Download blacklist
if [ ! -f "mm10-blacklist.v2.bed" ]; then
    echo "Downloading blacklist regions..."
    wget -q https://github.com/Boyle-Lab/Blacklist/raw/master/lists/mm10-blacklist.v2.bed.gz
    gunzip mm10-blacklist.v2.bed.gz
fi

# Build Bowtie2 index
if [ ! -f "mm10.1.bt2" ]; then
    echo "Building Bowtie2 index (30-60 min)..."
    bowtie2-build --threads 8 mm10.fa mm10
fi

cd ../..
echo "Reference setup complete"
