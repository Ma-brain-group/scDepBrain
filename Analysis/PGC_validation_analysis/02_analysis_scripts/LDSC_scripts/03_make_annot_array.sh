#!/bin/bash
# Make annotation files for each (cell type, chr).
# Array dimension: N_celltypes * 22 = (28 * 22) = up to 616 tasks
# Each task reads gene set -> SNPs in ±100kb window -> annot.gz
#SBATCH -J ldsc_annot
#SBATCH -p defq
#SBATCH -A hpcusers
#SBATCH -c 1
#SBATCH --mem=4G
#SBATCH -t 00:30:00
#SBATCH --array=1-616%80
#SBATCH -o /mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_LDSC_PGC2025_20260617/logs/annot-%A_%a.out
#SBATCH -e /mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_LDSC_PGC2025_20260617/logs/annot-%A_%a.err

set -euo pipefail
source /home/may2/miniconda3/etc/profile.d/conda.sh
conda activate ldsc

PP=/mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_LDSC_PGC2025_20260617
GS_DIR=$PP/celltype_genesets
COORD=$PP/refs/GRCh37/SEG_resources/coords/gencode.v41lift37.gene_name.coords.tsv
PLINK_PREFIX=$PP/refs/GRCh37/1000G_EUR_Phase3_plink/1000G.EUR.QC
MAKE_ANNOT=$PP/ldsc_software/make_annot.py

# Build list of (geneset, chr)
mapfile -t GS_FILES < <(ls $GS_DIR/*.geneset.txt 2>/dev/null | sort)
N_GS=${#GS_FILES[@]}
N_CHR=22

# Task index → (gs_idx, chr)
task=$SLURM_ARRAY_TASK_ID
gs_idx=$(( (task - 1) / N_CHR ))
chr=$(( (task - 1) % N_CHR + 1 ))

if [[ $gs_idx -ge $N_GS ]]; then
    echo "Task $task: out of range (gs_idx=$gs_idx, N_GS=$N_GS) — skipping"
    exit 0
fi

GS_FILE=${GS_FILES[$gs_idx]}
GS_NAME=$(basename $GS_FILE .geneset.txt)
ANNOT_DIR=$PP/annot/$GS_NAME
mkdir -p $ANNOT_DIR
OUT_ANNOT=$ANNOT_DIR/${GS_NAME}.${chr}.annot.gz

if [[ -e $OUT_ANNOT ]]; then
    echo "[$(date)] $GS_NAME chr $chr already done, skip"
    exit 0
fi

echo "[$(date)] make_annot: $GS_NAME chr $chr"
python $MAKE_ANNOT \
    --gene-set-file $GS_FILE \
    --gene-coord-file $COORD \
    --windowsize 100000 \
    --bimfile ${PLINK_PREFIX}.${chr}.bim \
    --annot-file $OUT_ANNOT

echo "[$(date)] done"
