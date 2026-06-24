#!/bin/bash
# S-LDSC v2: proper Finucane LDSC-SEG standard
# --ref-ld-chr baseline_LD_v2.2,all_gene,celltype
#SBATCH -J sldsc_v2
#SBATCH -p defq
#SBATCH -A hpcusers
#SBATCH -c 1
#SBATCH --mem=8G
#SBATCH -t 02:00:00
#SBATCH --array=1-56%30
#SBATCH -o /mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_LDSC_PGC2025_20260617/logs/sldsc-v2-%A_%a.out
#SBATCH -e /mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_LDSC_PGC2025_20260617/logs/sldsc-v2-%A_%a.err

set -euo pipefail
source /home/may2/miniconda3/etc/profile.d/conda.sh
conda activate ldsc

PP=/mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_LDSC_PGC2025_20260617
LDSC=$PP/ldsc_software/ldsc.py
BASELINE=$PP/refs/GRCh37/baseline_LD_v2.2/baselineLD.
ALL_GENE=$PP/annot/all_gene/all_gene.
FRQ=$PP/refs/GRCh37/1000G_Phase3_frq/1000G.EUR.QC.
WEIGHTS=$PP/refs/GRCh37/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.

declare -A GWAS=( [HW]=$PP/sumstats/HW.sumstats.gz \
                  [EUR]=$PP/from_miaot/sumstats/MDD_EUR_2025.hg38.sumstats.gz )
GWAS_KEYS=(HW EUR)

mapfile -t GS_FILES < <(ls $PP/celltype_genesets/major_*.geneset.txt $PP/celltype_genesets/subtype_*.geneset.txt 2>/dev/null | sort)
N_GS=${#GS_FILES[@]}

task=$SLURM_ARRAY_TASK_ID
gwas_idx=$(( (task - 1) / N_GS ))
gs_idx=$(( (task - 1) % N_GS ))
if [[ $gwas_idx -ge ${#GWAS_KEYS[@]} ]]; then exit 0; fi

KEY=${GWAS_KEYS[$gwas_idx]}
SUMSTATS=${GWAS[$KEY]}
GS_FILE=${GS_FILES[$gs_idx]}
GS_NAME=$(basename $GS_FILE .geneset.txt)
ANNOT_PREFIX=$PP/annot/${GS_NAME}/${GS_NAME}.
OUT=$PP/enrichment_v2/${KEY}__${GS_NAME}

mkdir -p $(dirname $OUT)

if [[ -e ${OUT}.results ]]; then
    echo "$KEY x $GS_NAME already done"; exit 0
fi

echo "[$(date)] $KEY x $GS_NAME (with all_gene control)"
python $LDSC \
    --h2 $SUMSTATS \
    --ref-ld-chr ${BASELINE},${ALL_GENE},${ANNOT_PREFIX} \
    --frqfile-chr $FRQ \
    --w-ld-chr $WEIGHTS \
    --overlap-annot \
    --print-coefficients \
    --out $OUT

echo "[$(date)] done"
