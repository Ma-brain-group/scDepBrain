# ============================================================
# Solution 1 — Step 1.3: GWAS_summary_input + SnpToGene
#
# Reads the 6-col scPagwas-format GWAS file (already prepared),
# runs GWAS_summary_input + SnpToGene, saves to gwas_pagwas.RData.
# This step is per-GWAS, runs once.
# ============================================================
suppressPackageStartupMessages({
  library(scPagwas); library(bigreadr)
})

OUT     <- "/mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_scPagwas_MDD2025_EUR_20260610"
GWAS_IN <- file.path(OUT, "pgc-mdd2025_eur_scPagwas_input.txt")
OUT_RDATA <- file.path(OUT, "gwas_pagwas.RData")

cat("[", as.character(Sys.time()), "] Reading GWAS ...\n")
gwas_data <- bigreadr::fread2(GWAS_IN)
cat("Rows:", nrow(gwas_data), "Cols:", paste(colnames(gwas_data), collapse=" "), "\n")

Pagwas <- list()
Pagwas <- scPagwas::GWAS_summary_input(Pagwas = Pagwas, gwas_data = gwas_data)

cat("[", as.character(Sys.time()), "] SnpToGene (marg = ±20 kb per manuscript) ...\n")
Pagwas$snp_gene_df <- scPagwas::SnpToGene(
  gwas_data = Pagwas$gwas_data,
  block_annotation = block_annotation,
  marg = 20000
)

# Also pre-attach rawPathway_list (used by Link_pathway_blocks_gwas)
Pagwas$rawPathway_list <- Genes_by_pathway_kegg

cat("Saving to", OUT_RDATA, "\n")
save(Pagwas, file = OUT_RDATA, compress = TRUE)
cat("DONE\n")
