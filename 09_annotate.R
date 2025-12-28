#!/usr/bin/env Rscript
# 09_annotate.R - Annotate peaks with ChIPseeker

suppressPackageStartupMessages({
  library(ChIPseeker)
  library(TxDb.Mmusculus.UCSC.mm10.knownGene)
  library(org.Mm.eg.db)
  library(clusterProfiler)
  library(GenomicRanges)
})

cat("Annotating peaks...\n")

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
dir.create("results/annotation", showWarnings = FALSE)
dir.create("results/figures", showWarnings = FALSE)

# Load differential peaks
adult_peaks <- readPeakFile("results/differential/adult_enriched_macs2.bed")
neonatal_peaks <- readPeakFile("results/differential/neonatal_enriched_macs2.bed")

# Annotate to genes
adult_anno <- annotatePeak(adult_peaks, 
                          tssRegion = c(-3000, 3000),
                          TxDb = txdb,
                          annoDb = "org.Mm.eg.db")

neonatal_anno <- annotatePeak(neonatal_peaks,
                             tssRegion = c(-3000, 3000),
                             TxDb = txdb,
                             annoDb = "org.Mm.eg.db")

# Save annotations
write.csv(as.data.frame(adult_anno), 
          "results/annotation/adult_peaks_annotated.csv",
          row.names = FALSE)

write.csv(as.data.frame(neonatal_anno),
          "results/annotation/neonatal_peaks_annotated.csv",
          row.names = FALSE)

# Extract gene lists
adult_genes <- unique(as.data.frame(adult_anno)$geneId)
neonatal_genes <- unique(as.data.frame(neonatal_anno)$geneId)

adult_symbols <- mapIds(org.Mm.eg.db, 
                       keys = adult_genes,
                       column = "SYMBOL",
                       keytype = "ENTREZID")

neonatal_symbols <- mapIds(org.Mm.eg.db,
                          keys = neonatal_genes,
                          column = "SYMBOL",
                          keytype = "ENTREZID")

writeLines(adult_symbols, "results/annotation/adult_genes.txt")
writeLines(neonatal_symbols, "results/annotation/neonatal_genes.txt")

# GO enrichment
adult_go <- enrichGO(gene = adult_genes,
                    OrgDb = org.Mm.eg.db,
                    ont = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.05)

write.csv(as.data.frame(adult_go),
          "results/annotation/GO_enrichment_adult.csv",
          row.names = FALSE)

# Plots
pdf("results/figures/peak_annotation.pdf", width = 10, height = 6)
plotAnnoBar(list(Adult = adult_anno, Neonatal = neonatal_anno))
dev.off()

pdf("results/figures/distance_to_TSS.pdf", width = 10, height = 6)
plotDistToTSS(list(Adult = adult_anno, Neonatal = neonatal_anno))
dev.off()

if (nrow(as.data.frame(adult_go)) > 0) {
  pdf("results/figures/GO_enrichment.pdf", width = 10, height = 8)
  dotplot(adult_go, showCategory = 20)
  dev.off()
}

cat("Annotation complete\n")
cat(sprintf("Adult peaks annotated: %d genes\n", length(adult_genes)))
cat(sprintf("Neonatal peaks annotated: %d genes\n", length(neonatal_genes)))
