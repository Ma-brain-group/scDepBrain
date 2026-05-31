#!/bin/bash
# ============================================================
# Extract MAGMA-format SNP location file from ieu-b-102 VCF.
#
# Why this script exists:
#   The original `mdd.location` is a 2-column file where chr & bp
#   were concatenated without a separator (e.g. "1548810499" — chr1
#   pos 548810499 or chr15 pos 48810499? ambiguous). MAGMA's
#   --snp-loc requires three separate columns: SNP CHR BP.
#   The source VCF (ieu-b-102.vcf, GRCh37) already has them split,
#   so we re-extract from there.
#
# Output: $OUT/mdd.location_v2  (space-separated: SNP CHR BP, no header)
#
# Filters:
#   - Keep autosomes only (chr 1-22) — standard for GWAS gene-based
#   - Drop rows with missing rsID ("." in VCF ID field)
#
# Runtime: ~3-5 min on a single core (805 MB plain-text VCF, awk pass)
# ============================================================
set -euo pipefail

PROJ=/mnt/isilon/gandal_lab/mayl/07_scDepBrain
DATA=$PROJ/MDD_2019_NN_ieu-b-102
OUT=$PROJ/results_MAGMA_MDD_20260531

VCF=$DATA/ieu-b-102.vcf
SNPLOC=$OUT/mdd.location_v2

echo "[$(date)] Extracting SNP/CHR/BP from $VCF ..."
awk 'BEGIN{OFS=" "} \
     !/^#/ && $1 ~ /^[0-9]+$/ && $3 != "." {print $3, $1, $2}' \
     "$VCF" > "$SNPLOC"

echo "[$(date)] Done."
echo "Output: $SNPLOC"
echo "Rows:   $(wc -l < "$SNPLOC")"
echo "Head:"
head -3 "$SNPLOC"
echo "Tail:"
tail -3 "$SNPLOC"
echo "Unique CHR values:"
awk '{print $2}' "$SNPLOC" | sort -u | tr '\n' ' '
echo ""
