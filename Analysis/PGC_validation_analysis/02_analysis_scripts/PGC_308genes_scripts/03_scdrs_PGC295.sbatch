#!/bin/bash
# scDRS on PGC 295 high-confidence MDD genes, weighted by MAGMA z-score
# Array task: 0=HW, 1=EUR, 2=DIV
#SBATCH -J PGC_scdrs
#SBATCH -p defq
#SBATCH -A hpcusers
#SBATCH -c 8
#SBATCH --mem=600G
#SBATCH -t 10-00:00:00
#SBATCH --array=0-2
#SBATCH -o /mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_PGC295_geneset_20260618/logs/scdrs-%A_%a.out
#SBATCH -e /mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_PGC295_geneset_20260618/logs/scdrs-%A_%a.err

set -euo pipefail
source /home/may2/miniconda3/etc/profile.d/conda.sh
conda activate pyscenic

PP=/mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_PGC295_geneset_20260618
H5AD=/mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_scDRS_20260618/data/MDD_scDRS.h5ad

TAGS=(HW EUR DIV)
TAG=${TAGS[$SLURM_ARRAY_TASK_ID]}
GS=$PP/gs/PGC295_${TAG}_weighted.gs
OUT=$PP/scores
mkdir -p $OUT

echo "[$(date)] $TAG: compute-score (PGC295 weighted, n_ctrl=500)"
scdrs compute-score \
    --h5ad-file $H5AD \
    --h5ad-species human \
    --gs-file $GS \
    --gs-species human \
    --out-folder $OUT \
    --flag-filter-data True \
    --flag-raw-count True \
    --n-ctrl 500

echo "[$(date)] $TAG DONE"
ls -lh $OUT/PGC295_${TAG}*
