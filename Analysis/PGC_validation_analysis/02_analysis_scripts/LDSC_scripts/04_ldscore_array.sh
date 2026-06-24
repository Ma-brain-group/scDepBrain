#!/bin/bash
# Calculate cell-type-specific LD scores per (cell type, chr).
# Paper-method: --ld-wind-cm 1, --thin-annot, restricted to HM3 SNPs via --print-snps
#SBATCH -J ldsc_l2
#SBATCH -p defq
#SBATCH -A hpcusers
#SBATCH -c 1
#SBATCH --mem=8G
#SBATCH -t 02:00:00
#SBATCH --array=1-616%80
#SBATCH -o /mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_LDSC_PGC2025_20260617/logs/l2-%A_%a.out
#SBATCH -e /mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_LDSC_PGC2025_20260617/logs/l2-%A_%a.err

set -euo pipefail
source /home/may2/miniconda3/etc/profile.d/conda.sh
conda activate ldsc

PP=/mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_LDSC_PGC2025_20260617
PLINK_PREFIX=$PP/refs/GRCh37/1000G_EUR_Phase3_plink/1000G.EUR.QC
HM3=$PP/refs/GRCh37/SEG_resources/w_hm3.print.snplist
LDSC=$PP/ldsc_software/ldsc.py

mapfile -t GS_FILES < <(ls $PP/celltype_genesets/*.geneset.txt 2>/dev/null | sort)
N_GS=${#GS_FILES[@]}
N_CHR=22
task=$SLURM_ARRAY_TASK_ID
gs_idx=$(( (task - 1) / N_CHR ))
chr=$(( (task - 1) % N_CHR + 1 ))

if [[ $gs_idx -ge $N_GS ]]; then exit 0; fi

GS_FILE=${GS_FILES[$gs_idx]}
GS_NAME=$(basename $GS_FILE .geneset.txt)
ANNOT_DIR=$PP/annot/$GS_NAME
ANNOT=$ANNOT_DIR/${GS_NAME}.${chr}.annot.gz
OUT_PREFIX=$ANNOT_DIR/${GS_NAME}.${chr}

if [[ -e ${OUT_PREFIX}.l2.ldscore.gz ]]; then
    echo "$GS_NAME chr $chr already done"
    exit 0
fi

[[ -e $ANNOT ]] || { echo "ERROR: annot $ANNOT not found"; exit 1; }

echo "[$(date)] LD score: $GS_NAME chr $chr"
python $LDSC \
    --l2 \
    --bfile ${PLINK_PREFIX}.${chr} \
    --ld-wind-cm 1 \
    --annot $ANNOT \
    --thin-annot \
    --print-snps $HM3 \
    --out $OUT_PREFIX

echo "[$(date)] done"
