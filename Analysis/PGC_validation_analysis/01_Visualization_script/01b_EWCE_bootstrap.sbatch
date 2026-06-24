#!/bin/bash
#SBATCH -J PGC_EWCE
#SBATCH -p defq
#SBATCH -A hpcusers
#SBATCH -c 4
#SBATCH --mem=400G
#SBATCH -t 04:00:00
#SBATCH -o /mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_PGC295_geneset_20260618/logs/EWCE-%j.out
#SBATCH -e /mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_PGC295_geneset_20260618/logs/EWCE-%j.err

set -euo pipefail
source /home/may2/miniconda3/etc/profile.d/conda.sh
conda activate pyscenic
Rscript /mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_PGC295_geneset_20260618/scripts/01b_EWCE_bootstrap.R
