#!/bin/bash
#SBATCH -J PGC_down
#SBATCH -p defq
#SBATCH -A hpcusers
#SBATCH -c 8
#SBATCH --mem=400G
#SBATCH -t 04:00:00
#SBATCH -o /mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_PGC295_geneset_20260618/logs/down-%j.out
#SBATCH -e /mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_PGC295_geneset_20260618/logs/down-%j.err

set -euo pipefail
source /home/may2/miniconda3/etc/profile.d/conda.sh
conda activate pyscenic

PP=/mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_PGC295_geneset_20260618
H5AD=/mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_scDRS_20260618/data/MDD_scDRS.h5ad
OUT=$PP/downstream
mkdir -p $OUT

for TAG in HW EUR DIV; do
    SCORE=$PP/scores/PGC295_${TAG}.full_score.gz
    [[ -e $SCORE ]] || { echo "WARN: $SCORE missing"; continue; }
    echo "[$(date)] $TAG: perform-downstream"
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
