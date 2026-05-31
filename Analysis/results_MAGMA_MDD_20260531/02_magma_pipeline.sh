#!/bin/bash
# ============================================================
# MAGMA pipeline for MDD GWAS (ieu-b-102, N=500,199)
#   Step 1: SNP -> gene annotation  (±20 kb window)
#   Step 2: gene-based association  (SNP-wise mean, 1000G EUR LD ref)
#   Step 3: KEGG pathway-based set analysis
#
# Fixes vs the original snippet:
#   - --bfile path:  g1000_eur/g1000_eur (NOT 1000G_data/g1000_eur)
#   - --pval N=500199 on same token group as the file (Howard 2019)
#   - --set-annot uses MSigDB 2025 KEGG legacy GMT (Entrez IDs),
#     because KEGG_for_MAGMA_annotated.txt does not exist
#   - --snp-loc uses freshly-extracted mdd.location_v2 (3 cols)
#   - Outputs go to dedicated results folder (project rule 5)
#
# Run as a SLURM job (project rule 2 — no heavy work on login node).
# Submit with:  sbatch 02_magma_pipeline.sh
# ============================================================
#SBATCH -J magma_mdd
#SBATCH -c 1
#SBATCH --mem=16G
#SBATCH -t 06:00:00
#SBATCH -o /mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_MAGMA_MDD_20260531/slurm-%j.out
#SBATCH -e /mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_MAGMA_MDD_20260531/slurm-%j.err
# NOTE: add  #SBATCH -p <partition>  and  #SBATCH -A <account>  lines
#       once you tell me your lab's SLURM partition/account names.

set -euo pipefail

# === Paths ===
PROJ=/mnt/isilon/gandal_lab/mayl/07_scDepBrain
MAGMA_DIR=/mnt/isilon/gandal_lab/mayl/01_GWAS_tools/MAGMA
DATA=$PROJ/MDD_2019_NN_ieu-b-102
OUT=$PROJ/results_MAGMA_MDD_20260531

MAGMA=$MAGMA_DIR/magma
SNPLOC=$OUT/mdd.location_v2                                       # produced by 01_extract_snp_loc.sh
PVAL=$DATA/mdd.Pval
GENELOC=$MAGMA_DIR/NCBI37.3.gene.loc
BFILE=$MAGMA_DIR/g1000_eur/g1000_eur
KEGG=$MAGMA_DIR/msigdb_v2025.1.Hs_geneset_GMTs/c2.cp.kegg_legacy.v2025.1.Hs.entrez.gmt
N_SAMPLES=500199

# === Sanity checks ===
for f in "$MAGMA" "$SNPLOC" "$PVAL" "$GENELOC" "${BFILE}.bed" "$KEGG"; do
    [[ -e "$f" ]] || { echo "ERROR: missing $f"; exit 1; }
done

# === Step 0 — auto-run extraction if v2 location is missing ===
if [[ ! -s "$SNPLOC" ]]; then
    echo "[$(date)] mdd.location_v2 not found, running 01_extract_snp_loc.sh ..."
    bash "$OUT/01_extract_snp_loc.sh"
fi

# === Step 1: MAGMA annotation (SNP -> gene, ±20 kb) ===
echo "[$(date)] === Step 1: SNP -> gene annotation ==="
"$MAGMA" \
    --snp-loc "$SNPLOC" \
    --annotate window=20,20 \
    --gene-loc "$GENELOC" \
    --out "$OUT/mdd.hg19_SNP_Gene_annotation"

# === Step 2: gene-based association analysis (N = Howard 2019 meta) ===
echo "[$(date)] === Step 2: gene-based analysis (N=$N_SAMPLES) ==="
"$MAGMA" \
    --bfile "$BFILE" \
    --pval "$PVAL" N=$N_SAMPLES \
    --gene-annot "$OUT/mdd.hg19_SNP_Gene_annotation.genes.annot" \
    --out "$OUT/mdd.hg19_SNP_Gene_Analysis_P"

# === Step 3: KEGG pathway-based set analysis ===
echo "[$(date)] === Step 3: KEGG (MSigDB 2025 legacy) set analysis ==="
"$MAGMA" \
    --gene-results "$OUT/mdd.hg19_SNP_Gene_Analysis_P.genes.raw" \
    --set-annot "$KEGG" \
    --out "$OUT/mdd.KEGG_MAGMA_Gene_set_results"

echo "[$(date)] === All MAGMA steps complete ==="
ls -lh "$OUT/"
