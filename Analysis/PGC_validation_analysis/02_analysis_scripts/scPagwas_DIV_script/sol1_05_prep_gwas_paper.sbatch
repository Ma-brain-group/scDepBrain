#!/bin/bash
#SBATCH -J div_pp_05
#SBATCH -p defq
#SBATCH -A hpcusers
#SBATCH -c 2
#SBATCH --mem=100G
#SBATCH -t 10-00:00:00
#SBATCH -o /mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_scPagwas_DIV_paperparams_20260615/slurm-%j.out
#SBATCH -e /mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_scPagwas_DIV_paperparams_20260615/slurm-%j.err

set -euo pipefail
cd /mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_scPagwas_DIV_paperparams_20260615
source /home/may2/miniconda3/etc/profile.d/conda.sh
conda activate pyscenic
Rscript /mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_scPagwas_DIV_paperparams_20260615/sol1_05_prep_gwas_paper.R
