#!/bin/bash
#SBATCH -J sol1_06_eur
#SBATCH -p defq
#SBATCH -A hpcusers
#SBATCH -c 4
#SBATCH --mem=400G
#SBATCH -t 10-00:00:00
#SBATCH --array=1-5
#SBATCH -o /mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_scPagwas_MDD2025_EUR_20260610/slurm-%A_%a.out
#SBATCH -e /mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_scPagwas_MDD2025_EUR_20260610/slurm-%A_%a.err

set -euo pipefail
cd /mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_scPagwas_MDD2025_EUR_20260610
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1
source /home/may2/miniconda3/etc/profile.d/conda.sh
conda activate pyscenic
echo "[$(date)] Running chunk $SLURM_ARRAY_TASK_ID on $(hostname)"
Rscript /mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_scPagwas_MDD2025_EUR_20260610/sol1_06_chunk_regression.R
echo "[$(date)] Chunk $SLURM_ARRAY_TASK_ID done"
