#!/bin/bash
#SBATCH -J div_pp_06
#SBATCH -p defq
#SBATCH -A hpcusers
#SBATCH -c 4
#SBATCH --mem=600G
#SBATCH -t 10-00:00:00
#SBATCH --array=1-5
#SBATCH -o /mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_scPagwas_DIV_paperparams_20260615/slurm-%A_%a.out
#SBATCH -e /mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_scPagwas_DIV_paperparams_20260615/slurm-%A_%a.err

set -euo pipefail
cd /mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_scPagwas_DIV_paperparams_20260615
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1
source /home/may2/miniconda3/etc/profile.d/conda.sh
conda activate pyscenic
echo "[$(date)] chunk $SLURM_ARRAY_TASK_ID on $(hostname)"
Rscript /mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_scPagwas_DIV_paperparams_20260615/sol1_06_chunk_regression.R
