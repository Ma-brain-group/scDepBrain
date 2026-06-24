#!/bin/bash
# Run stratified LDSC: 3 GWAS × N cell types
# Each task: partition heritability with baseline_LD_v2.2 + cell_type annotation
#SBATCH -J sldsc
#SBATCH -p defq
#SBATCH -A hpcusers
#SBATCH -c 1
#SBATCH --mem=8G
#SBATCH -t 02:00:00
#SBATCH --array=1-84%30
#SBATCH -o /mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_LDSC_PGC2025_20260617/logs/sldsc-%A_%a.out
#SBATCH -e /mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_LDSC_PGC2025_20260617/logs/sldsc-%A_%a.err

set -euo pipefail
source /home/may2/miniconda3/etc/profile.d/conda.sh
conda activate ldsc

PP=/mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_LDSC_PGC2025_20260617
LDSC=$PP/ldsc_software/ldsc.py
BASELINE=$PP/refs/GRCh37/baseline_LD_v2.2/baselineLD.
FRQ=$PP/refs/GRCh37/1000G_Phase3_frq/1000G.EUR.QC.
WEIGHTS=$PP/refs/GRCh37/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.

GWAS_LIST=(HW EUR DIV)
mapfile -t GS_FILES < <(ls $PP/celltype_genesets/*.geneset.txt 2>/dev/null | sort)
N_GS=${#GS_FILES[@]}
N_GWAS=${#GWAS_LIST[@]}

task=$SLURM_ARRAY_TASK_ID
gwas_idx=$(( (task - 1) / N_GS ))
gs_idx=$(( (task - 1) % N_GS ))
if [[ $gwas_idx -ge $N_GWAS ]]; then exit 0; fi

GWAS=${GWAS_LIST[$gwas_idx]}
SUMSTATS=$PP/sumstats/${GWAS}.sumstats.gz
GS_FILE=${GS_FILES[$gs_idx]}
GS_NAME=$(basename $GS_FILE .geneset.txt)
ANNOT_PREFIX=$PP/annot/${GS_NAME}/${GS_NAME}.
OUT=$PP/enrichment/${GWAS}__${GS_NAME}

mkdir -p $(dirname $OUT)

if [[ -e ${OUT}.results ]]; then
    echo "$GWAS x $GS_NAME already done"
    exit 0
fi

echo "[$(date)] S-LDSC: $GWAS x $GS_NAME"
python $LDSC \
    --h2 $SUMSTATS \
    --ref-ld-chr ${BASELINE},${ANNOT_PREFIX} \
    --frqfile-chr $FRQ \
    --w-ld-chr $WEIGHTS \
    --overlap-annot \
    --print-coefficients \
    --out $OUT

echo "[$(date)] done"
