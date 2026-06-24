#!/bin/bash
# ============================================================
# Step 3: MAGMA cell-type enrichment
# For each (GWAS × gene-set-file): run MAGMA gene-set analysis.
# 3 GWAS × 2 gene-sets = 6 outputs.
# Depends on: gene-based for all 3 GWAS + cell-type gene sets
# Runtime: ~5-10 min total (MAGMA gene-set is fast).
# ============================================================
#SBATCH -J magma_enrichment
#SBATCH -p defq
#SBATCH -A hpcusers
#SBATCH -c 2
#SBATCH --mem=16G
#SBATCH -t 02:00:00
#SBATCH -o /mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_MAGMA_PGC2025_20260617/logs/slurm-%j.out
#SBATCH -e /mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_MAGMA_PGC2025_20260617/logs/slurm-%j.err

set -euo pipefail
PROJ=/mnt/isilon/gandal_lab/mayl/07_scDepBrain
MAGMA_DIR=/mnt/isilon/gandal_lab/mayl/01_GWAS_tools/MAGMA
PP=$PROJ/results_MAGMA_PGC2025_20260617
ENR=$PP/enrichment

GENES_HW=$PP/HW/hw.genes.raw
GENES_EUR=$PP/EUR/pgc2025_eur.genes.raw
GENES_DIV=$PP/DIV/pgc2025_div.genes.raw

SET_ANNO=$PP/celltype_genesets/anno_top10.txt
SET_FINAL=$PP/celltype_genesets/final_anno_top10.txt

# Sanity check inputs
for f in $GENES_HW $GENES_EUR $GENES_DIV $SET_ANNO $SET_FINAL; do
  [[ -e "$f" ]] || { echo "ERROR missing: $f"; exit 1; }
done
echo "All inputs present."

# Enrichment runs
declare -A gwas_files=( [HW]=$GENES_HW [EUR]=$GENES_EUR [DIV]=$GENES_DIV )
declare -A set_files=( [anno]=$SET_ANNO [final_anno]=$SET_FINAL )

for tag in HW EUR DIV; do
  for setname in anno final_anno; do
    out=$ENR/${tag}_${setname}
    echo "[$(date)] $tag x $setname -> $out ..."
    $MAGMA_DIR/magma \
      --gene-results ${gwas_files[$tag]} \
      --set-annot   ${set_files[$setname]} \
      --out $out
  done
done

echo "[$(date)] DONE. Outputs:"
ls -lh $ENR/*.gsa.out 2>/dev/null
