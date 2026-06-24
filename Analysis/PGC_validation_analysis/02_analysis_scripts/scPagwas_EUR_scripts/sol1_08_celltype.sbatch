#!/bin/bash
#SBATCH -J sol1_08_eur_ct
#SBATCH -p defq
#SBATCH -A hpcusers
#SBATCH -c 4
#SBATCH --mem=400G
#SBATCH -t 10-00:00:00
#SBATCH -o /mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_scPagwas_MDD2025_EUR_20260610/slurm-%j.out
#SBATCH -e /mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_scPagwas_MDD2025_EUR_20260610/slurm-%j.err

set -euo pipefail
cd /mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_scPagwas_MDD2025_EUR_20260610
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1
source /home/may2/miniconda3/etc/profile.d/conda.sh
conda activate pyscenic
Rscript /mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_scPagwas_MDD2025_EUR_20260610/sol1_08_celltype.R
