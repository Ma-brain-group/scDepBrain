#!/bin/bash
# Re-run MAGMA gene-based for HW (Howard 2019) with the same params as
# the EUR/DIV runs so all 3 are strictly comparable:
#   - window = 50,50
#   - gene-loc extendedMHCexcluded
#   - reuse mdd.location_v2 + mdd.Pval_v2 (already MAF>0.01, no MHC, has rsID)
#   - HW has no per-SNP N column — use single N=500199 (Howard 2019)
#SBATCH -J magma_hw_pp
#SBATCH -p defq
#SBATCH -A hpcusers
#SBATCH -c 2
#SBATCH --mem=32G
#SBATCH -t 10-00:00:00
#SBATCH -o /mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_MAGMA_PGC2025_20260617/logs/slurm-%j.out
#SBATCH -e /mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_MAGMA_PGC2025_20260617/logs/slurm-%j.err

set -euo pipefail
PROJ=/mnt/isilon/gandal_lab/mayl/07_scDepBrain
MAGMA_DIR=/mnt/isilon/gandal_lab/mayl/01_GWAS_tools/MAGMA
PP=$PROJ/results_MAGMA_PGC2025_20260617
OUT=$PP/EUR/..   # stage HW in its own subfolder

# Put HW outputs in EUR_DIV-sibling folder
HW_OUT=$PP/HW
mkdir -p $HW_OUT
SNP_LOC=$PROJ/results_MAGMA_MDD_20260531/mdd.location_v2
PVAL=$PROJ/results_MAGMA_MDD_20260531/mdd.Pval_v2
GENE_LOC=$PP/NCBI37.3.gene.loc.extendedMHCexcluded
ANNOT_PREFIX=$HW_OUT/hw.annot
GENES_PREFIX=$HW_OUT/hw

echo "[$(date)] MAGMA annotation (window=50,50, HW) ..."
$MAGMA_DIR/magma \
  --snp-loc $SNP_LOC \
  --annotate window=50,50 \
  --gene-loc $GENE_LOC \
  --out $ANNOT_PREFIX

echo "[$(date)] MAGMA gene-based (N=500199 single, HW) ..."
$MAGMA_DIR/magma \
  --bfile $MAGMA_DIR/g1000_eur/g1000_eur \
  --pval $PVAL N=500199 \
  --gene-annot $ANNOT_PREFIX.genes.annot \
  --out $GENES_PREFIX

echo "[$(date)] DONE — $GENES_PREFIX.genes.raw"
