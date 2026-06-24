#!/bin/bash
#SBATCH -J div_pp_07
#SBATCH -p defq
#SBATCH -A hpcusers
#SBATCH -c 4
#SBATCH --mem=800G
#SBATCH -t 10-00:00:00
#SBATCH -o /mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_scPagwas_DIV_paperparams_20260615/slurm-%j.out
#SBATCH -e /mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_scPagwas_DIV_paperparams_20260615/slurm-%j.err

set -euo pipefail
cd /mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_scPagwas_DIV_paperparams_20260615
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1
source /home/may2/miniconda3/etc/profile.d/conda.sh
conda activate pyscenic
Rscript /mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_scPagwas_DIV_paperparams_20260615/sol1_07_finalize_paper.R
