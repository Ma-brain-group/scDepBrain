#!/bin/bash
# all_gene control: make annot + LD score in ONE array task per chr
#SBATCH -J all_gene
#SBATCH -p defq
#SBATCH -A hpcusers
#SBATCH -c 1
#SBATCH --mem=8G
#SBATCH -t 02:00:00
#SBATCH --array=1-22
#SBATCH -o /mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_LDSC_PGC2025_20260617/logs/allgene-%A_%a.out
#SBATCH -e /mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_LDSC_PGC2025_20260617/logs/allgene-%A_%a.err

set -euo pipefail
source /home/may2/miniconda3/etc/profile.d/conda.sh
conda activate ldsc

PP=/mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_LDSC_PGC2025_20260617
LDSC=$PP/ldsc_software/ldsc.py
MAKE_ANNOT=$PP/ldsc_software/make_annot.py
GS_FILE=$PP/celltype_genesets/all_gene.geneset.txt
COORD=$PP/refs/GRCh37/SEG_resources/coords/gencode.v41lift37.gene_name.coords.tsv
PLINK_PREFIX=$PP/refs/GRCh37/1000G_EUR_Phase3_plink/1000G.EUR.QC
HM3=$PP/refs/GRCh37/SEG_resources/w_hm3.print.snplist

ANNOT_DIR=$PP/annot/all_gene
mkdir -p $ANNOT_DIR

chr=$SLURM_ARRAY_TASK_ID
ANNOT=$ANNOT_DIR/all_gene.${chr}.annot.gz
OUT_PREFIX=$ANNOT_DIR/all_gene.${chr}

if [[ ! -e $ANNOT ]]; then
    echo "[$(date)] make_annot chr $chr"
    python $MAKE_ANNOT \
        --gene-set-file $GS_FILE \
        --gene-coord-file $COORD \
        --windowsize 100000 \
        --bimfile ${PLINK_PREFIX}.${chr}.bim \
        --annot-file $ANNOT
fi

if [[ ! -e ${OUT_PREFIX}.l2.ldscore.gz ]]; then
    echo "[$(date)] LD score chr $chr"
    python $LDSC \
        --l2 \
        --bfile ${PLINK_PREFIX}.${chr} \
        --ld-wind-cm 1 \
        --annot $ANNOT \
        --thin-annot \
        --print-snps $HM3 \
        --out $OUT_PREFIX
fi
echo "[$(date)] done"
