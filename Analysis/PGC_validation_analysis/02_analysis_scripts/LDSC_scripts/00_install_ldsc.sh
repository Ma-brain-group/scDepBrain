#!/bin/bash
# Step 0: Install LDSC + conda env (login node, network)
# Note: LDSC uses Python 2.7. We create a dedicated conda env.
set -euo pipefail

PROJ=/mnt/isilon/gandal_lab/mayl/07_scDepBrain
LDSC_DIR=$PROJ/results_LDSC_PGC2025_20260617
LDSC_REPO=$LDSC_DIR/ldsc_software

# Clone LDSC if not already
if [[ ! -d $LDSC_REPO/.git ]]; then
    echo "[$(date)] Cloning LDSC from GitHub ..."
    git clone https://github.com/bulik/ldsc.git $LDSC_REPO
fi

# Create conda env from environment.yml (gives correct deps)
source /home/may2/miniconda3/etc/profile.d/conda.sh
if ! conda env list | awk '{print $1}' | grep -qx ldsc; then
    echo "[$(date)] Creating conda env 'ldsc' from environment.yml ..."
    conda env create --file $LDSC_REPO/environment.yml
fi

# Test
echo "[$(date)] Test ldsc.py:"
conda activate ldsc
python $LDSC_REPO/ldsc.py -h 2>&1 | head -3
python $LDSC_REPO/munge_sumstats.py -h 2>&1 | head -3
echo "[$(date)] LDSC ready."
