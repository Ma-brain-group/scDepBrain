#!/bin/bash
#SBATCH -J scdrs_down_EUR
#SBATCH -p defq
#SBATCH -A hpcusers
#SBATCH -c 8
#SBATCH --mem=400G
#SBATCH -t 04:00:00
#SBATCH -o /mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_scDRS_20260618/logs/down-EUR-%j.out
#SBATCH -e /mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_scDRS_20260618/logs/down-EUR-%j.err

set -euo pipefail
source /home/may2/miniconda3/etc/profile.d/conda.sh
conda activate pyscenic

PP=/mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_scDRS_20260618
scdrs perform-downstream \
    --h5ad-file $PP/data/MDD_scDRS.h5ad \
    --score-file $PP/scores/EUR.full_score.gz \
    --out-folder $PP/downstream \
    --group-analysis anno,final_anno \
    --flag-filter-data True \
    --flag-raw-count True

echo "DONE"
ls -lh $PP/downstream/EUR*
