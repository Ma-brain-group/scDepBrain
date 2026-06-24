#!/bin/bash
# Build GENE/HW/EUR wide pval TSV, then run scdrs munge-gs
# scdrs format: row=gene, col=trait
set -euo pipefail
source /home/may2/miniconda3/etc/profile.d/conda.sh
conda activate pyscenic

PP=/mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_scDRS_20260618
MAGMA_DIR=/mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_MAGMA_PGC2025_20260617
GENE_LOC=/mnt/isilon/gandal_lab/mayl/01_GWAS_tools/MAGMA/NCBI37.3.gene.loc

mkdir -p $PP/gs $PP/data

PVAL_WIDE=$PP/data/HW_EUR_gene_pval.tsv

echo "[$(date)] Build wide-format pval TSV: GENE  HW  EUR"
python <<'EOF'
import pandas as pd
gloc = pd.read_csv("/mnt/isilon/gandal_lab/mayl/01_GWAS_tools/MAGMA/NCBI37.3.gene.loc",
                   sep=r"\s+", header=None,
                   names=["entrez","chr","start","end","strand","symbol"],
                   dtype={"entrez":str})
sym = dict(zip(gloc["entrez"], gloc["symbol"]))

def load(path, tag):
    d = pd.read_csv(path, sep=r"\s+", dtype={"GENE":str})
    d["SYMBOL"] = d["GENE"].map(sym)
    d = d.dropna(subset=["SYMBOL"])
    d = d[["SYMBOL","P"]].rename(columns={"SYMBOL":"GENE","P":tag})
    # Dedup symbol (keep min p)
    d = d.sort_values(tag).drop_duplicates("GENE", keep="first")
    return d.set_index("GENE")

hw  = load("/mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_MAGMA_PGC2025_20260617/HW/hw.genes.out", "HW")
eur = load("/mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_MAGMA_PGC2025_20260617/EUR/pgc2025_eur.genes.out", "EUR")
merged = hw.join(eur, how="outer")
merged.index.name = "GENE"
merged.to_csv("/mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_scDRS_20260618/data/HW_EUR_gene_pval.tsv",
              sep="\t")
print(f"Wide pval: {merged.shape[0]} genes  HW non-NA: {merged['HW'].notna().sum()}  EUR non-NA: {merged['EUR'].notna().sum()}")
EOF
echo ""
head -2 $PVAL_WIDE
echo ""

echo "[$(date)] scdrs munge-gs (top 1000 z-score weighted)"
scdrs munge-gs \
    --out-file $PP/gs/HW_EUR.gs \
    --pval-file $PVAL_WIDE \
    --weight zscore --n-max 1000

echo ""
echo "==Resulting gs file=="
ls -lh $PP/gs/
cat $PP/gs/HW_EUR.gs | head -3 | awk '{print substr($0,1,200)}'
