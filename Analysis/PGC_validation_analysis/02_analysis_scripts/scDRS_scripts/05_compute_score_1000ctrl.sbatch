#!/bin/bash
# 1000-ctrl version (re-run after 500 ctrl results in)
# Higher memory to handle 2x control matrix
#SBATCH -J scdrs_1k
#SBATCH -p defq
#SBATCH -A hpcusers
#SBATCH -c 8
#SBATCH --mem=600G
#SBATCH -t 10-00:00:00
#SBATCH --array=0-2
#SBATCH -o /mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_scDRS_20260618/logs/score-1k-%A_%a.out
#SBATCH -e /mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_scDRS_20260618/logs/score-1k-%A_%a.err

set -euo pipefail
source /home/may2/miniconda3/etc/profile.d/conda.sh
conda activate pyscenic

PP=/mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_scDRS_20260618
H5AD=$PP/data/MDD_scDRS.h5ad

# Array index → (trait, gs_file)
GS_FILES=("$PP/gs/HW_EUR.gs" "$PP/gs/HW_EUR.gs" "$PP/gs/DIV.gs")
# We split HW_EUR into 2 array tasks by trait via filter; for simplicity run both traits in one job (array_idx 0 = HW_EUR, 2 traits)
# Use single gs per task and separate output dir
OUT_DIRS=("$PP/scores_1k_HW_EUR" "$PP/scores_1k_HW_EUR" "$PP/scores_1k_DIV")
GS_KEYS=("HW_EUR" "HW_EUR" "DIV")

# Simpler: 2 jobs total — array idx 0=HW_EUR, 1=DIV (skip duplicate index 2)
case $SLURM_ARRAY_TASK_ID in
  0) GS=$PP/gs/HW_EUR.gs; OUT=$PP/scores_1k; LABEL=HW_EUR ;;
  1) GS=$PP/gs/DIV.gs;    OUT=$PP/scores_1k; LABEL=DIV ;;
  *) echo "array_idx out of range"; exit 0 ;;
esac

mkdir -p $OUT

echo "[$(date)] $LABEL: compute-score (n_ctrl=1000)"
scdrs compute-score \
    --h5ad-file $H5AD \
    --h5ad-species human \
    --gs-file $GS \
    --gs-species human \
    --out-folder $OUT \
    --flag-filter-data True \
    --flag-raw-count True \
    --n-ctrl 1000

echo "[$(date)] $LABEL DONE"
ls -lh $OUT/
