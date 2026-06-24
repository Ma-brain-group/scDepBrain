#!/bin/bash
#SBATCH -J scdrs_down
#SBATCH -p defq
#SBATCH -A hpcusers
#SBATCH -c 8
#SBATCH --mem=400G
#SBATCH -t 04:00:00
#SBATCH -o /mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_scDRS_20260618/logs/down-%j.out
#SBATCH -e /mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_scDRS_20260618/logs/down-%j.err

set -euo pipefail
source /home/may2/miniconda3/etc/profile.d/conda.sh
conda activate pyscenic

PP=/mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_scDRS_20260618
H5AD=$PP/data/MDD_scDRS.h5ad
SCORE_DIR=$PP/scores
OUT=$PP/downstream
mkdir -p $OUT

# Run per trait: HW + EUR; group analysis on anno (8 major) + final_anno (20 subtypes)
for TRAIT in HW EUR; do
    SCORE=$SCORE_DIR/${TRAIT}.full_score.gz
    if [[ ! -e $SCORE ]]; then
        echo "WARN: $SCORE missing, skip"; continue
    fi
    echo "[$(date)] $TRAIT perform-downstream (anno + final_anno)"
    scdrs perform-downstream \
        --h5ad-file $H5AD \
        --score-file $SCORE \
        --out-folder $OUT \
        --group-analysis anno,final_anno \
        --flag-filter-data True \
        --flag-raw-count True
done

echo "[$(date)] DONE"
ls -lh $OUT/
