#!/bin/bash
#SBATCH -J scdrs_score
#SBATCH -p defq
#SBATCH -A hpcusers
#SBATCH -c 16
#SBATCH --mem=600G
#SBATCH -t 10-00:00:00
#SBATCH -o /mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_scDRS_20260618/logs/score-%j.out
#SBATCH -e /mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_scDRS_20260618/logs/score-%j.err

set -euo pipefail
source /home/may2/miniconda3/etc/profile.d/conda.sh
conda activate pyscenic

PP=/mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_scDRS_20260618
H5AD=$PP/data/MDD_scDRS.h5ad
GS=$PP/gs/HW_EUR.gs
OUT=$PP/scores
mkdir -p $OUT

# Optional covariate file (skip if absent — scdrs computes anyway)
COV_OPT=""
if [[ -e $PP/data/cov.tsv ]]; then COV_OPT="--cov-file $PP/data/cov.tsv"; fi

echo "[$(date)] scdrs compute-score (HW + EUR, top 1000 genes z-weighted)"
scdrs compute-score \
    --h5ad-file $H5AD \
    --h5ad-species human \
    --gs-file $GS \
    --gs-species human \
    --out-folder $OUT \
    $COV_OPT \
    --flag-filter-data True \
    --flag-raw-count True \
    --n-ctrl 500

echo "[$(date)] DONE"
ls -lh $OUT/
