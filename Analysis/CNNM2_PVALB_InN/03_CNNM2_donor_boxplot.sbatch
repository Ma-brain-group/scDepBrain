#!/bin/bash
#SBATCH -J CNNM2_donor_box
#SBATCH -p defq
#SBATCH -A hpcusers
#SBATCH -c 4
#SBATCH --mem=200G
#SBATCH -t 02:00:00
#SBATCH -o /mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_CNNM2_replication_20260619/logs/slurm-%j.out
#SBATCH -e /mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_CNNM2_replication_20260619/logs/slurm-%j.err

set -euo pipefail
source /home/may2/miniconda3/etc/profile.d/conda.sh
conda activate pyscenic
Rscript /mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_CNNM2_replication_20260619/scripts/03_CNNM2_donor_boxplot.R
