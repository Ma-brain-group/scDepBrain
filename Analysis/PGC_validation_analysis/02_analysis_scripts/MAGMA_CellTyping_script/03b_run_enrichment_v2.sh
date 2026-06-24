#!/bin/bash
# v2 enrichment: include BOTH --gene-covar (Bryois standard) and --set-annot (sanity check)
#SBATCH -J magma_enrich_v2
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

declare -A gwas_files=( [HW]=$GENES_HW [EUR]=$GENES_EUR [DIV]=$GENES_DIV )
declare -a granularities=( anno final_anno )

for tag in HW EUR DIV; do
  for gran in "${granularities[@]}"; do
    # 1) gene-covar (Bryois MAGMA_CellTyping standard)
    out_cov=$ENR/${tag}_${gran}_covar_v2
    cov_file=$PP/celltype_genesets/${gran}_gene_covar_v2.txt
    echo "[$(date)] $tag x $gran [gene-covar] ..."
    $MAGMA_DIR/magma \
      --gene-results ${gwas_files[$tag]} \
      --gene-covar $cov_file \
      --out $out_cov

    # 2) set-annot top 10% (sanity / direct comparison)
    out_set=$ENR/${tag}_${gran}_set_v2
    set_file=$PP/celltype_genesets/${gran}_top10_v2.txt
    echo "[$(date)] $tag x $gran [set-annot top10] ..."
    $MAGMA_DIR/magma \
      --gene-results ${gwas_files[$tag]} \
      --set-annot   $set_file \
      --out $out_set
  done
done

echo "[$(date)] DONE. Outputs:"
ls -lh $ENR/*_v2.gsa.out 2>/dev/null
