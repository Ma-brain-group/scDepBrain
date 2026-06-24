#!/bin/bash
#SBATCH -J scdrs_down_DIV
#SBATCH -p defq
#SBATCH -A hpcusers
#SBATCH -c 8
#SBATCH --mem=400G
#SBATCH -t 04:00:00
#SBATCH -o /mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_scDRS_20260618/logs/down-DIV-%j.out
#SBATCH -e /mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_scDRS_20260618/logs/down-DIV-%j.err

set -euo pipefail
source /home/may2/miniconda3/etc/profile.d/conda.sh
conda activate pyscenic

PP=/mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_scDRS_20260618
H5AD=$PP/data/MDD_scDRS.h5ad
SCORE=$PP/scores/DIV.full_score.gz
OUT=$PP/downstream
mkdir -p $OUT

if [[ ! -e $SCORE ]]; then
    echo "ERROR: $SCORE missing"; exit 1
fi

echo "[$(date)] DIV perform-downstream (anno + final_anno)"
scdrs perform-downstream \
    --h5ad-file $H5AD \
    --score-file $SCORE \
    --out-folder $OUT \
    --group-analysis anno,final_anno \
    --flag-filter-data True \
    --flag-raw-count True

echo "[$(date)] DONE"
ls -lh $OUT/DIV*
