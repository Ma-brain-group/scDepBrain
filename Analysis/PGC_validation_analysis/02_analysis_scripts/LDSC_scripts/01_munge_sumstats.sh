#!/bin/bash
#SBATCH -J ldsc_munge
#SBATCH -p defq
#SBATCH -A hpcusers
#SBATCH -c 2
#SBATCH --mem=32G
#SBATCH -t 04:00:00
#SBATCH -o /mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_LDSC_PGC2025_20260617/logs/slurm-%j.out
#SBATCH -e /mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_LDSC_PGC2025_20260617/logs/slurm-%j.err

set -euo pipefail

source /home/may2/miniconda3/etc/profile.d/conda.sh
conda activate ldsc

PROJ=/mnt/isilon/gandal_lab/mayl/07_scDepBrain
LDSC_DIR=$PROJ/results_LDSC_PGC2025_20260617
LDSC_SOFT=$LDSC_DIR/ldsc_software
SUMOUT=$LDSC_DIR/sumstats
HM3=$LDSC_DIR/refs/GRCh37/w_hm3.snplist

mkdir -p $SUMOUT

# ============================================================
# Step 1A: HW (Howard 2019) — convert VCF -> tsv  (N=500199 fixed)
# ============================================================
HW_VCF=$PROJ/MDD_2019_NN_ieu-b-102/ieu-b-102.vcf
HW_TSV=$SUMOUT/HW_raw.tsv

if [[ ! -e $HW_TSV ]]; then
    echo "[$(date)] HW: parsing VCF -> tsv ..."
    awk 'BEGIN { OFS="\t"; print "SNP","A1","A2","BETA","SE","P","N" }
        /^#/ { next }
        {
            chrom=$1; pos=$2; rsid=$3; ref=$4; alt=$5
            if (chrom !~ /^([0-9]|1[0-9]|2[0-2])$/) next
            if (rsid !~ /^rs/) next
            n=split($9,  fmt, ":")
            split($10, vals, ":")
            es=""; se=""; lp=""
            for (i=1;i<=n;i++) {
                if (fmt[i]=="ES") es=vals[i]
                if (fmt[i]=="SE") se=vals[i]
                if (fmt[i]=="LP") lp=vals[i]
            }
            if (es=="" || se=="" || lp=="") next
            p = (lp+0 > 300) ? "1e-300" : sprintf("%.10g", 10^(-(lp+0)))
            print rsid, alt, ref, es, se, p, 500199
        }' $HW_VCF > $HW_TSV
    wc -l $HW_TSV
fi

# Munge
if [[ ! -e $SUMOUT/HW.sumstats.gz ]]; then
    echo "[$(date)] HW: munge_sumstats.py ..."
    python $LDSC_SOFT/munge_sumstats.py \
        --sumstats $HW_TSV \
        --out $SUMOUT/HW \
        --merge-alleles $HM3 \
        --signed-sumstats BETA,0
fi

# ============================================================
# Step 1B: PGC 2025 EUR — convert tsv.gz -> tsv  (N=NEFF per-SNP)
# ============================================================
EUR_GZ=$PROJ/MDD_2025_Cell_GWAS_summary/pgc-mdd2025_no23andMe_eur_v3-49-24-11.tsv.gz
EUR_TSV=$SUMOUT/EUR_raw.tsv

if [[ ! -e $EUR_TSV ]]; then
    echo "[$(date)] EUR: parsing tsv.gz -> tsv ..."
    zcat $EUR_GZ | awk 'BEGIN { OFS="\t" }
        /^##/ { next }
        /^#CHROM/ {
            for (i=1;i<=NF;i++) {
                c=$i; sub("#","",c)
                if (c=="CHROM") C=i; if (c=="POS") P=i; if (c=="ID") I=i
                if (c=="EA") EA=i; if (c=="NEA") NEA=i
                if (c=="BETA") B=i; if (c=="SE") S=i; if (c=="PVAL") V=i; if (c=="NEFF") N=i
            }
            print "SNP","A1","A2","BETA","SE","P","N"
            next
        }
        {
            chr=$C; if (chr !~ /^([0-9]|1[0-9]|2[0-2])$/) next
            rs=$I; if (rs !~ /^rs/) next
            print rs, $EA, $NEA, $B, $S, $V, int($N+0.5)
        }' > $EUR_TSV
    wc -l $EUR_TSV
fi

if [[ ! -e $SUMOUT/EUR.sumstats.gz ]]; then
    echo "[$(date)] EUR: munge_sumstats.py ..."
    python $LDSC_SOFT/munge_sumstats.py \
        --sumstats $EUR_TSV \
        --out $SUMOUT/EUR \
        --merge-alleles $HM3 \
        --signed-sumstats BETA,0
fi

# ============================================================
# Step 1C: PGC 2025 DIV (trans-anc) — same parser
# ============================================================
DIV_GZ=$PROJ/MDD_2025_Cell_GWAS_summary/pgc-mdd2025_no23andMe_div_v3-49-46-01.tsv.gz
DIV_TSV=$SUMOUT/DIV_raw.tsv

if [[ ! -e $DIV_TSV ]]; then
    echo "[$(date)] DIV: parsing tsv.gz -> tsv ..."
    zcat $DIV_GZ | awk 'BEGIN { OFS="\t" }
        /^##/ { next }
        /^#CHROM/ {
            for (i=1;i<=NF;i++) {
                c=$i; sub("#","",c)
                if (c=="CHROM") C=i; if (c=="POS") P=i; if (c=="ID") I=i
                if (c=="EA") EA=i; if (c=="NEA") NEA=i
                if (c=="BETA") B=i; if (c=="SE") S=i; if (c=="PVAL") V=i; if (c=="NEFF") N=i
            }
            print "SNP","A1","A2","BETA","SE","P","N"
            next
        }
        {
            chr=$C; if (chr !~ /^([0-9]|1[0-9]|2[0-2])$/) next
            rs=$I; if (rs !~ /^rs/) next
            print rs, $EA, $NEA, $B, $S, $V, int($N+0.5)
        }' > $DIV_TSV
    wc -l $DIV_TSV
fi

if [[ ! -e $SUMOUT/DIV.sumstats.gz ]]; then
    echo "[$(date)] DIV: munge_sumstats.py ..."
    python $LDSC_SOFT/munge_sumstats.py \
        --sumstats $DIV_TSV \
        --out $SUMOUT/DIV \
        --merge-alleles $HM3 \
        --signed-sumstats BETA,0
fi

echo ""
echo "[$(date)] DONE — munged sumstats:"
ls -lh $SUMOUT/*.sumstats.gz
