#!/bin/bash
# ============================================================
# Convert PGC MDD 2025 (Trans-ancestry (DIV, multi-ancestry), no 23andMe) GWAS sumstats
# from 17-col PGC tsv.gz to scPagwas 6-col space-delim format.
#
# Source columns (PGC):
#   #CHROM POS ID EA NEA BETA SE PVAL FCAS FCON IMPINFO NEFF NCAS NCON HETI HETDF HETPVAL
# Output columns (scPagwas):
#   chrom pos rsid se beta maf
#
# Manuscript Methods filters:
#   - autosomes only (chr 1-22)
#   - require rsID (drop ".")
#   - exclude extended MHC (chr6: 25-35 Mb)
#   - MAF > 0.01  (computed as sample-size-weighted EAF then min(EAF,1-EAF))
#
# Runtime: ~3-5 min single-core awk pass on 223 MB gzipped (~1.5 GB raw)
# ============================================================
#SBATCH -J prep_div
#SBATCH -p defq
#SBATCH -A hpcusers
#SBATCH -c 1
#SBATCH --mem=4G
#SBATCH -t 00:30:00
#SBATCH -o /mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_scPagwas_MDD2025_TransAnc_20260610/slurm-%j.out
#SBATCH -e /mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_scPagwas_MDD2025_TransAnc_20260610/slurm-%j.err

set -euo pipefail
PROJ=/mnt/isilon/gandal_lab/mayl/07_scDepBrain
DATA=$PROJ/MDD_2025_Cell_GWAS_summary
OUT=$PROJ/results_scPagwas_MDD2025_TransAnc_20260610
GWAS_IN=$DATA/pgc-mdd2025_no23andMe_div_v3-49-46-01.tsv.gz
GWAS_OUT=$OUT/pgc-mdd2025_div_scPagwas_input.txt

echo "[$(date)] Converting $GWAS_IN -> $GWAS_OUT ..."

zcat "$GWAS_IN" | awk '
  BEGIN { OFS = " "; print "chrom","pos","rsid","se","beta","maf" }
  /^##/ { next }
  /^#CHROM/ {
      for (i=1; i<=NF; i++) {
          c = $i; sub("#","",c)
          if (c=="CHROM") C=i
          if (c=="POS")   P=i
          if (c=="ID")    I=i
          if (c=="BETA")  B=i
          if (c=="SE")    S=i
          if (c=="FCAS")  FA=i
          if (c=="FCON")  FC=i
          if (c=="NCAS")  Nca=i
          if (c=="NCON")  Nco=i
      }
      next
  }
  {
      chr = $C; pos = $P + 0; rs = $I
      if (chr !~ /^[0-9]+$/) next
      chr_n = chr + 0
      if (chr_n < 1 || chr_n > 22) next
      if (rs == "." || rs == "") next
      if (chr_n == 6 && pos >= 25000000 && pos <= 35000000) next

      beta = $B + 0; se = $S + 0
      fcas = $FA + 0; fcon = $FC + 0
      nca  = $Nca + 0; nco = $Nco + 0
      if (se <= 0 || (nca + nco) <= 0) next

      eaf = (fcas*nca + fcon*nco) / (nca + nco)
      maf = (eaf < 0.5) ? eaf : (1.0 - eaf)
      if (maf <= 0.01) next

      print chr_n, pos, rs, se, beta, maf
  }' > "$GWAS_OUT"

echo "[$(date)] Done"
echo "Rows:" && wc -l "$GWAS_OUT"
echo "Head:" && head -5 "$GWAS_OUT"
echo "CHR distribution:" && awk 'NR>1 {print $1}' "$GWAS_OUT" | sort -n | uniq -c
echo "MAF summary:" && awk 'NR>1 {sum+=$6; if ($6<min || min==0) min=$6; if ($6>max) max=$6; n++}
     END {print "  n=", n, "min=", min, "max=", max, "mean=", sum/n}' "$GWAS_OUT"
