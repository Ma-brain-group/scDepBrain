#!/bin/bash
#SBATCH -J PGC_modscore
#SBATCH -p defq
#SBATCH -A hpcusers
#SBATCH -c 4
#SBATCH --mem=400G
#SBATCH -t 04:00:00
#SBATCH -o /mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_PGC295_geneset_20260618/logs/modscore-%j.out
#SBATCH -e /mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_PGC295_geneset_20260618/logs/modscore-%j.err

set -euo pipefail
source /home/may2/miniconda3/etc/profile.d/conda.sh
conda activate pyscenic
Rscript /mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_PGC295_geneset_20260618/scripts/02_AddModuleScore.R
