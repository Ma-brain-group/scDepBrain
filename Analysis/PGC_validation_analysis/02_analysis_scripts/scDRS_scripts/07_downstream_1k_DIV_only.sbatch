#!/bin/bash
#SBATCH -J down_1k_DIV
#SBATCH -p defq
#SBATCH -A hpcusers
#SBATCH -c 8
#SBATCH --mem=400G
#SBATCH -t 06:00:00
#SBATCH -o /mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_scDRS_20260618/logs/down-1k-DIV-%j.out
#SBATCH -e /mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_scDRS_20260618/logs/down-1k-DIV-%j.err

set -euo pipefail
source /home/may2/miniconda3/etc/profile.d/conda.sh
conda activate pyscenic

PP=/mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_scDRS_20260618
scdrs perform-downstream \
    --h5ad-file $PP/data/MDD_scDRS.h5ad \
    --score-file $PP/scores_1k/DIV.full_score.gz \
    --out-folder $PP/downstream_1k \
    --group-analysis anno,final_anno \
    --flag-filter-data True \
    --flag-raw-count True

echo "DONE"
ls -lh $PP/downstream_1k/DIV*
