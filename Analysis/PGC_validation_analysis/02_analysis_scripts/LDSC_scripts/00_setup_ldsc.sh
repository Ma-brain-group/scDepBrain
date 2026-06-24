#!/bin/bash
# ============================================================
# Step 0: Setup LDSC environment + clone repo + download refs
# Run on login node (network access required for clone + downloads).
# Total runtime ~15-30 min (mostly download time).
# Reference files go to: results_LDSC_PGC2025_20260617/refs/
# ============================================================
set -euo pipefail

PROJ=/mnt/isilon/gandal_lab/mayl/07_scDepBrain
LDSC_DIR=$PROJ/results_LDSC_PGC2025_20260617
REFS=$LDSC_DIR/refs

# ---- Part A: conda env + LDSC clone ----
echo "[$(date)] Part A: setup conda env 'ldsc' + clone LDSC repo"
source /home/may2/miniconda3/etc/profile.d/conda.sh

if ! conda env list | awk '{print $1}' | grep -qx ldsc; then
    echo "Creating conda env 'ldsc' with Python 2.7 ..."
    conda create -n ldsc -c bioconda -c conda-forge python=2.7 \
        bitarray=0.8 nose=1.3 'pybedtools' pandas=0.20 numpy=1.16 scipy=1.2 -y
fi

LDSC_REPO=$LDSC_DIR/ldsc_software
if [[ ! -d $LDSC_REPO/.git ]]; then
    echo "Cloning LDSC from GitHub ..."
    git clone https://github.com/bulik/ldsc.git $LDSC_REPO
fi

# Quick smoke test
conda activate ldsc
python $LDSC_REPO/ldsc.py -h | head -3
echo "LDSC installed at $LDSC_REPO"

# ---- Part B: download reference files ----
echo ""
echo "[$(date)] Part B: download LDSC reference files into $REFS"
cd $REFS

# All from Broad alkesgroup public LDSCORE server.
BASE=https://alkesgroup.broadinstitute.org/LDSCORE

declare -a urls=(
    "$BASE/1000G_Phase3_plinkfiles.tgz"          # 1000G EUR plink files (for LD score computation)
    "$BASE/1000G_Phase3_baselineLD_v2.2_ldscores.tgz"   # baseline LD scores
    "$BASE/1000G_Phase3_frq.tgz"                 # MAF info
    "$BASE/1000G_Phase3_weights_hm3_no_MHC.tgz"  # weights, HapMap3 SNPs, no MHC
    "$BASE/w_hm3.snplist.bz2"                    # HapMap3 SNP list (for munge_sumstats)
    "$BASE/eur_w_ld_chr.tgz"                     # EUR LD weights (for h2)
)

for u in "${urls[@]}"; do
    f=$(basename $u)
    if [[ ! -e $REFS/$f && ! -d $REFS/${f%.tgz} && ! -e $REFS/${f%.bz2} ]]; then
        echo "[$(date)] downloading $f ..."
        wget --quiet --show-progress -O $REFS/$f $u
    else
        echo "[$(date)] $f already present, skipping"
    fi
done

# ---- Extract archives ----
echo ""
echo "[$(date)] Extracting archives ..."
cd $REFS
for tar in *.tgz; do
    [[ -e $tar ]] || continue
    outdir=${tar%.tgz}
    if [[ ! -d $outdir ]]; then
        echo "  tar -xzf $tar"
        tar -xzf $tar
    fi
done

# bz2 single file
[[ -e w_hm3.snplist.bz2 && ! -e w_hm3.snplist ]] && bunzip2 -k w_hm3.snplist.bz2

echo ""
echo "[$(date)] DONE — references in $REFS"
ls -lah $REFS
