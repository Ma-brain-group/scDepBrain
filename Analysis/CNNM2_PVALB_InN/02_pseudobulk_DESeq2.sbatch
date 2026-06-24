#!/bin/bash
#SBATCH -J CNNM2_pseudobulk
#SBATCH -p defq
#SBATCH -A hpcusers
#SBATCH -c 8
#SBATCH --mem=300G
#SBATCH -t 06:00:00
#SBATCH -o /mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_CNNM2_replication_20260619/logs/slurm-%j.out
#SBATCH -e /mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_CNNM2_replication_20260619/logs/slurm-%j.err

set -euo pipefail
source /home/may2/miniconda3/etc/profile.d/conda.sh
conda activate pyscenic

# Ensure DESeq2 available (Matrix.utils dropped — use sparse matrix multiplication instead)
Rscript -e 'if (!"DESeq2" %in% rownames(installed.packages())) { BiocManager::install("DESeq2", update=FALSE, ask=FALSE) }' 2>&1 | tail -5

Rscript /mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_CNNM2_replication_20260619/scripts/02_pseudobulk_DESeq2.R
