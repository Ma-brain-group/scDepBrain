#!/bin/bash
# ============================================================
# MAGMA gene-based analysis for PGC 2025 EUR GWAS
# Following user's reference pipeline:
#   - window = 50,50 (±50 kb)
#   - --gene-loc extendedMHCexcluded
#   - --pval ncol=N (per-SNP NEFF)
#
# Inputs:  pgc-mdd2025_no23andMe_eur_v3-49-24-11.tsv.gz
# Outputs (in EUR/):
#   pgc2025_eur.snp_loc      (SNP CHR BP)
#   pgc2025_eur.pval         (SNP P N)   N = NEFF column
#   pgc2025_eur.genes.annot  (annotation step)
#   pgc2025_eur.genes.raw    (gene-based output, used downstream)
# ============================================================
#SBATCH -J magma_eur_pp
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
OUT=$PP/EUR
GENE_LOC=$PP/NCBI37.3.gene.loc.extendedMHCexcluded
GWAS_IN=$PROJ/MDD_2025_Cell_GWAS_summary/pgc-mdd2025_no23andMe_eur_v3-49-24-11.tsv.gz

SNP_LOC=$OUT/pgc2025_eur.snp_loc
PVAL=$OUT/pgc2025_eur.pval
ANNOT_PREFIX=$OUT/pgc2025_eur.annot
GENES_PREFIX=$OUT/pgc2025_eur

echo "[$(date)] Step 1: extract SNP loc + pval (with N=NEFF) from $GWAS_IN ..."
# Parse header to find column indices; write 2 files in single awk pass
zcat "$GWAS_IN" | awk -v OUT_LOC="$SNP_LOC" -v OUT_P="$PVAL" '
  BEGIN { OFS = " "
          print "SNP","CHR","BP"   > OUT_LOC
          print "SNP","P","N"      > OUT_P }
  /^##/ { next }
  /^#CHROM/ {
      for (i=1; i<=NF; i++) {
          c=$i; sub("#","",c)
          if (c=="CHROM") C=i
          if (c=="POS")   P=i
          if (c=="ID")    I=i
          if (c=="PVAL")  V=i
          if (c=="NEFF")  N=i
      }
      next
  }
  {
      chr=$C; pos=$P+0; rs=$I; pv=$V; nf=$N
      if (chr !~ /^[0-9]+$/) next
      chr_n = chr + 0
      if (chr_n < 1 || chr_n > 22) next
      if (rs == "." || rs == "") next
      if (pv == "" || pv == ".") next
      if (nf == "" || nf == ".") next
      # MAGMA chokes on p=0; floor to 1e-300
      if (pv+0 <= 0) pv = "1e-300"
      print rs, chr_n, pos >> OUT_LOC
      printf "%s %s %d\n", rs, pv, int(nf+0.5) >> OUT_P
  }'
echo "[$(date)] Wrote $(wc -l < $SNP_LOC) snp_loc rows, $(wc -l < $PVAL) pval rows"
echo "Head snp_loc:" && head -3 $SNP_LOC
echo "Head pval:"    && head -3 $PVAL

echo ""
echo "[$(date)] Step 2: MAGMA annotation (window=50,50) ..."
$MAGMA_DIR/magma \
  --snp-loc $SNP_LOC \
  --annotate window=50,50 \
  --gene-loc $GENE_LOC \
  --out $ANNOT_PREFIX

echo ""
echo "[$(date)] Step 3: MAGMA gene-based (ncol=N) ..."
$MAGMA_DIR/magma \
  --bfile $MAGMA_DIR/g1000_eur/g1000_eur \
  --pval $PVAL ncol=N \
  --gene-annot $ANNOT_PREFIX.genes.annot \
  --out $GENES_PREFIX

echo ""
echo "[$(date)] DONE — output: $GENES_PREFIX.genes.raw"
ls -lh $OUT
