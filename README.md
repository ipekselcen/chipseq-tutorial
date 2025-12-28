# ChIP-seq Analysis: H4K8ac in Adult vs Neonatal OPCs

Complete ChIP-seq analysis pipeline from raw reads to biological insights, analyzing histone H4 lysine 8 acetylation patterns in oligodendrocyte progenitor cells.

## Dataset

**GEO Accession:** GSE263808  
**Publication:** Dansu & Selcen et al. (2024), *Journal of Cell Biology*  
**DOI:** 10.1083/jcb.202308064  
**Organism:** Mouse (*Mus musculus*)  
**Cell Type:** Oligodendrocyte progenitor cells (PDGFRα+)  
**Comparison:** Adult OPCs (P60) vs Neonatal OPCs (P5)  
**Target:** Histone H4 Lysine 8 Acetylation (H4K8ac)  
**Sequencing:** Illumina HiSeq 4000 (paired-end)  
**Samples:** 3 biological replicates per condition + Input controls

## Quick Start

```bash
# Clone repository
git clone https://github.com/ipekselcen/chipseq-tutorial.git
cd chipseq-tutorial

# Create environment
conda env create -f environment.yml
conda activate chipseq

# Run complete pipeline
bash scripts/run_pipeline.sh
```

**Runtime:** 6-8 hours (depending on system)

## Workflow

```
1. Download data          → 01_download_data.sh
2. Quality control        → 02_quality_control.sh
3. Trim adapters          → 03_trim_reads.sh
4. Align reads            → 04_align_reads.sh
5. Filter & deduplicate   → 05_filter_bams.sh
6. Peak calling           → 06_call_peaks.sh
7. Differential binding   → 07_differential_binding.R
8. Annotation & analysis  → 08_annotation.R
9. Visualization          → 09_visualization.R
```

## Key Biological Findings

From Dansu & Selcen et al. (2024):

**35,820 differentially bound regions** (FDR < 0.01)

* Greater H4K8ac occupancy in adult OPCs at:
  * **Progenitor genes** (*Hes5*, *Gpr17*)
  * **Metabolic genes** (*Txnip*, *Ptgds*)
  * **Myelin genes** (*Cnp*, *Mog*)

* Peak distribution:
  * 43.2% at promoters
  * 34.2% at introns
  * 18.3% at distal intergenic regions

* >60% overlap between H4K8ac occupancy and upregulated transcripts

## Tools

**Preprocessing:** FastQC, Trim Galore  
**Alignment:** Bowtie2  
**Post-processing:** SAMtools, Picard  
**Peak Calling:** MACS2  
**Quality Control:** deepTools, phantompeakqualtools  
**Differential Binding:** DiffBind  
**Annotation:** ChIPseeker, ChIPpeakAnno  
**Functional Analysis:** GREAT, clusterProfiler  
**Visualization:** ggplot2, Gviz, IGV

## Output Files

```
results/
├── qc/                    # Quality control reports
│   ├── fastqc/           # Raw read QC
│   ├── alignment_stats/  # Mapping statistics
│   └── fingerprint/      # ChIP quality metrics
├── aligned/               # BAM files (filtered & indexed)
├── peaks/                 # MACS2 peak calls
│   ├── narrowPeak/       # Peak coordinates
│   ├── summits/          # Peak summits
│   └── bedgraph/         # Coverage tracks
├── differential/          # DiffBind results
│   ├── consensus_peaks/  # Merged peak sets
│   └── diff_binding/     # Adult vs Neonatal
├── annotation/            # Genomic feature annotation
│   ├── peak_annotation/  # ChIPseeker results
│   └── functional/       # GO/pathway enrichment
└── figures/               # All plots (300 DPI PNG/PDF)
    ├── qc_plots/
    ├── genome_browser/
    ├── heatmaps/
    └── enrichment/
```

## Requirements

**Software:**

* R ≥ 4.0
* Bowtie2 ≥ 2.4.0
* SAMtools ≥ 1.15
* FastQC ≥ 0.11.9
* MACS2 ≥ 2.2.7
* deepTools ≥ 3.5.0
* bedtools ≥ 2.30.0

**R Packages:**

```r
# Install Bioconductor packages
BiocManager::install(c("DiffBind", "ChIPseeker", "ChIPpeakAnno",
                       "clusterProfiler", "org.Mm.eg.db", 
                       "TxDb.Mmusculus.UCSC.mm10.knownGene",
                       "GenomicRanges", "rtracklayer", "Gviz"))

# Install CRAN packages
install.packages(c("ggplot2", "pheatmap", "dplyr", "tidyr", 
                   "viridis", "RColorBrewer"))
```

**Reference Genome:**

* Mouse genome (mm10/GRCm38)
* Download automatically via script or manually from UCSC/Ensembl

## Project Structure

```
chipseq-tutorial/
├── README.md
├── GETTING_STARTED.md
├── environment.yml
├── data/
│   ├── raw/              # FASTQ files (from SRA)
│   └── reference/        # Genome and blacklist
├── scripts/
│   ├── 01_download_data.sh
│   ├── 02_quality_control.sh
│   ├── 03_trim_reads.sh
│   ├── 04_align_reads.sh
│   ├── 05_filter_bams.sh
│   ├── 06_call_peaks.sh
│   ├── 07_differential_binding.R
│   ├── 08_annotation.R
│   ├── 09_visualization.R
│   ├── utils/
│   │   ├── download_genome.sh
│   │   └── sample_info.txt
│   └── run_pipeline.sh
└── results/
    ├── qc/
    ├── aligned/
    ├── peaks/
    ├── differential/
    ├── annotation/
    └── figures/
```

## Analysis Highlights

### 1. Peak Calling Strategy
- Used MACS2 with Input control normalization
- Broad peak calling for histone modifications
- FDR < 0.01 threshold

### 2. Quality Control Metrics
- NSC (Normalized Strand Cross-correlation)
- RSC (Relative Strand Cross-correlation)  
- FRiP (Fraction of Reads in Peaks)
- Library complexity assessment

### 3. Differential Binding Analysis
- DiffBind consensus peak set generation
- DESeq2-based differential analysis
- Multiple testing correction (FDR)

### 4. Functional Annotation
- ChIPseeker for genomic feature distribution
- GREAT for functional enrichment
- Gene ontology and pathway analysis

## Biological Context

Adult oligodendrocyte progenitors (aOPCs) differ from neonatal progenitors (nOPCs) in their:

* **Proliferation rate** - Lower in adults
* **Transcriptional profile** - More similar to mature oligodendrocytes
* **Epigenetic landscape** - Enriched H4K8ac marks

This tutorial demonstrates how **H4K8ac**, an activating histone mark, occupies genomic regions corresponding to genes that regulate:

1. Progenitor identity maintenance
2. Metabolic processes
3. Myelin gene expression

Understanding these epigenetic differences is crucial for developing targeted myelin repair therapies for neurodegenerative diseases.

## Citation

**Data from:**  
Dansu DK*, Selcen I*, Sauma S, Prentice E, Huang D, Li M, Moyon S, Casaccia P. (2024). "Histone H4 acetylation differentially modulates proliferation in adult oligodendrocyte progenitors." *Journal of Cell Biology*, 223(11):e202308064. doi: 10.1083/jcb.202308064

*Co-first authors

**Dataset:**  
Available at GEO: [GSE263808](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE263808)

## About This Tutorial

Created as part of a comprehensive bioinformatics portfolio demonstrating expertise in:
* ChIP-seq experimental design and analysis
* Epigenetics and chromatin biology
* R/Bioconductor workflows
* Reproducible research practices

**Author:** Ipek Selcen, PhD  
**Website:** [ipekselcen.github.io](https://ipekselcen.github.io)  
**Contact:** [GitHub](https://github.com/ipekselcen)

## License

This tutorial is licensed under MIT License. The original data is from Dansu & Selcen et al. (2024) and subject to GEO data access policies.
