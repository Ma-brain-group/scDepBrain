#!/bin/bash
# Submit: sbatch --dependency=afterok:<prep_jobid> 02_run_scPagwas.sbatch
#SBATCH -J scPagwas_EUR
#SBATCH -p defq
#SBATCH -A hpcusers
#SBATCH -c 16
#SBATCH --mem=600G
#SBATCH -t 10-00:00:00
#SBATCH -o /mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_scPagwas_MDD2025_EUR_20260610/slurm-%j.out
#SBATCH -e /mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_scPagwas_MDD2025_EUR_20260610/slurm-%j.err

set -euo pipefail
cd /mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_scPagwas_MDD2025_EUR_20260610
echo "[$(date)] Starting scPagwas (PGC 2025 EUR) on $(hostname)"
echo "[$(date)] Memory snapshot:" && free -g | head -3

# Avoid "two levels of parallelism" — see MDD_20260531/04_run_scPagwas.sbatch
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1
export GOTO_NUM_THREADS=1

source /home/may2/miniconda3/etc/profile.d/conda.sh
conda activate pyscenic
echo "[$(date)] Using R: $(which Rscript)"

Rscript /mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_scPagwas_MDD2025_EUR_20260610/02_run_scPagwas.R

echo "[$(date)] Done"
echo "[$(date)] Final memory:" && free -g | head -3
