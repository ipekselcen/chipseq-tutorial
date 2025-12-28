#!/usr/bin/env Rscript
# 10_visualize.R - Create plots and figures

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(RColorBrewer)
})

cat("Creating visualizations...\n")

# Read peak counts
peak_files <- list.files("results/peaks", pattern = "broadPeak$", full.names = TRUE)
peak_counts <- data.frame(
  sample = basename(peak_files) %>% gsub("_peaks.broadPeak", "", .),
  peaks = sapply(peak_files, function(f) length(readLines(f)))
)

# Separate by condition
peak_counts$condition <- ifelse(grepl("aOPC", peak_counts$sample), "Adult", "Neonatal")

# Plot peak counts
pdf("results/figures/peak_counts.pdf", width = 8, height = 6)
ggplot(peak_counts, aes(x = sample, y = peaks, fill = condition)) +
  geom_bar(stat = "identity") +
  scale_fill_brewer(palette = "Set1") +
  labs(title = "H4K8ac Peaks per Sample",
       x = "Sample",
       y = "Number of Peaks") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

# Read differential results
adult_count <- length(readLines("results/differential/adult_enriched_macs2.bed"))
neonatal_count <- length(readLines("results/differential/neonatal_enriched_macs2.bed"))

diff_data <- data.frame(
  condition = c("Adult", "Neonatal"),
  sites = c(adult_count, neonatal_count)
)

# Plot differential sites
pdf("results/figures/differential_sites.pdf", width = 6, height = 6)
ggplot(diff_data, aes(x = condition, y = sites, fill = condition)) +
  geom_bar(stat = "identity") +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "Differential H4K8ac Sites",
       x = "Condition",
       y = "Number of Sites") +
  theme_minimal() +
  geom_text(aes(label = sites), vjust = -0.5)
dev.off()

cat("Visualizations complete\n")
cat("Plots saved to results/figures/\n")
