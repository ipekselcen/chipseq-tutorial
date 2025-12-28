# ChIP-seq Analysis: H4K8ac in Oligodendrocyte Progenitors

Complete ChIP-seq workflow from raw reads to differential binding analysis.

## Dataset

**Publication:** Dansu & Selcen et al. (2024) *J Cell Biol* 223(11):e202308064  
**GEO:** GSE263808  
**Comparison:** Adult OPCs (P60) vs Neonatal OPCs (P5)  
**Target:** H4K8ac histone modification  
**Samples:** 3 biological replicates + Input controls

## Quick Start

```bash
git clone https://github.com/ipekselcen/chipseq-tutorial.git
cd chipseq-tutorial

conda env create -f environment.yml
conda activate chipseq

bash scripts/run_all.sh
```

## Workflow

```
Raw FASTQ → QC → Trim → Align → Filter → Call Peaks → Differential → Annotate
```

| Step | Script | Tool |
|------|--------|------|
| 1 | `01_setup.sh` | Download reference genome |
| 2 | `02_download.sh` | Get data from SRA |
| 3 | `03_qc.sh` | FastQC + MultiQC |
| 4 | `04_trim.sh` | Trim Galore |
| 5 | `05_align.sh` | Bowtie2 |
| 6 | `06_filter.sh` | SAMtools + Picard |
| 7 | `07_peaks.sh` | MACS2 |
| 8 | `08_differential.sh` | MACS2 + diffReps |
| 9 | `09_annotate.R` | ChIPseeker |
| 10 | `10_visualize.R` | Plots |

## Key Results

- **35,820** differential H4K8ac sites (FDR < 0.01)
- **95%** enriched in adult OPCs
- Peaks at progenitor genes (*Hes5*, *Gpr17*), metabolic genes (*Txnip*), myelin genes (*Cnp*, *Mog*)

## Output

```
results/
├── qc/                  # Quality reports
├── aligned/             # BAM files
├── peaks/               # MACS2 peaks + bigWig
├── differential/        # Adult vs neonatal sites
├── annotation/          # Gene lists + GO terms
└── figures/             # All plots
```

## Citation

```
Dansu DK*, Selcen I*, et al. (2024). Histone H4 acetylation differentially 
modulates proliferation in adult oligodendrocyte progenitors. 
J Cell Biol 223(11):e202308064.
```

## Contact

Ipek Selcen, PhD  
[ipekselcen.github.io](https://ipekselcen.github.io) | [@ipekselcen](https://github.com/ipekselcen)
