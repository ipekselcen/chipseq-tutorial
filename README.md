# H4K8ac ChIP-seq in oligodendrocyte progenitor cells

This repository contains scripts I’m organizing around the H4K8ac ChIP-seq work in our study of adult and neonatal oligodendrocyte progenitor cells. The data are available through [GEO: GSE263808](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE263808).

The scripts are a work in progress, not a tested, ready-to-run pipeline. They are currently in the repository root, while `run_all.sh` expects them under `scripts/`, so the runner will not work as written.

## What’s in the repository

The numbered shell and R scripts cover reference setup, data download, quality control, trimming, alignment, filtering, peak calling, differential analysis, annotation, and plotting. `environment.yml` lists the software environment.

## Study

Dansu DK*, Selcen I*, et al. (2024). Histone H4 acetylation differentially modulates proliferation in adult oligodendrocyte progenitors. *Journal of Cell Biology* 223(11):e202308064. [doi:10.1083/jcb.202308064](https://doi.org/10.1083/jcb.202308064)

*Co-first authors.
