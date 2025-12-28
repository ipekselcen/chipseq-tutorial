#!/bin/bash
# 08_differential.sh - Differential binding analysis

set -euo pipefail

echo "Running differential binding analysis..."

mkdir -p results/differential

# Method 1: MACS2 bdgdiff
echo "Method 1: MACS2 bdgdiff..."

macs2 bdgdiff \
    --t1 results/peaks/aOPC_rep1_treat_pileup.bdg \
         results/peaks/aOPC_rep2_treat_pileup.bdg \
         results/peaks/aOPC_rep3_treat_pileup.bdg \
    --t2 results/peaks/nOPC_rep1_treat_pileup.bdg \
         results/peaks/nOPC_rep2_treat_pileup.bdg \
         results/peaks/nOPC_rep3_treat_pileup.bdg \
    --c1 results/peaks/aOPC_rep1_control_lambda.bdg \
    --c2 results/peaks/nOPC_rep1_control_lambda.bdg \
    -C 3.0 \
    -o results/differential/macs2_adult_vs_neonatal

# Rename outputs
mv results/differential/macs2_adult_vs_neonatal_c3.0_cond1.bed \
   results/differential/adult_enriched_macs2.bed
mv results/differential/macs2_adult_vs_neonatal_c3.0_cond2.bed \
   results/differential/neonatal_enriched_macs2.bed

# Method 2: diffReps
echo "Method 2: diffReps..."

diffReps.pl \
    --treatment results/aligned/aOPC_rep1_ChIP_filtered.bam \
                results/aligned/aOPC_rep2_ChIP_filtered.bam \
                results/aligned/aOPC_rep3_ChIP_filtered.bam \
    --control results/aligned/nOPC_rep1_ChIP_filtered.bam \
              results/aligned/nOPC_rep2_ChIP_filtered.bam \
              results/aligned/nOPC_rep3_ChIP_filtered.bam \
    --window 1000 \
    --step 100 \
    --pval 0.01 \
    --fdr 0.01 \
    --nproc 8 \
    --report results/differential/diffreps_results \
    2>&1 | tee results/differential/diffreps.log

# Extract significant sites
awk '$11 < 0.01 && $7 > 0' results/differential/diffreps_results.txt \
    > results/differential/adult_enriched_diffreps.txt
awk '$11 < 0.01 && $7 < 0' results/differential/diffreps_results.txt \
    > results/differential/neonatal_enriched_diffreps.txt

# Summary
echo ""
echo "Results:"
echo "Adult enriched (MACS2): $(wc -l < results/differential/adult_enriched_macs2.bed)"
echo "Neonatal enriched (MACS2): $(wc -l < results/differential/neonatal_enriched_macs2.bed)"
echo "Adult enriched (diffReps): $(wc -l < results/differential/adult_enriched_diffreps.txt)"
echo "Neonatal enriched (diffReps): $(wc -l < results/differential/neonatal_enriched_diffreps.txt)"

echo "Differential analysis complete"
