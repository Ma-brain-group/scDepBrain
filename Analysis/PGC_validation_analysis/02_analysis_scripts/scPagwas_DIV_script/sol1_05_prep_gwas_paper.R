# ============================================================
# HW paper-params version: GWAS_summary_input + SnpToGene
# with maf_filter = 0.1 (paper) and marg = 10000 (paper).
# ============================================================
suppressPackageStartupMessages({
  library(scPagwas); library(bigreadr)
})

OUT     <- "/mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_scPagwas_DIV_paperparams_20260615"
GWAS_IN <- file.path(OUT, "pgc-mdd2025_div_scPagwas_input.txt")
OUT_RDATA <- file.path(OUT, "gwas_pagwas.RData")

cat("[", as.character(Sys.time()), "] Reading GWAS ...\n")
gwas_data <- bigreadr::fread2(GWAS_IN)
cat("Input rows:", nrow(gwas_data), "\n")

Pagwas <- list()
Pagwas <- scPagwas::GWAS_summary_input(
  Pagwas = Pagwas,
  gwas_data = gwas_data,
  maf_filter = 0.1          # PAPER PARAM
)
cat("After maf>0.1 filter, gwas_data rows:", nrow(Pagwas$gwas_data), "\n")

cat("[", as.character(Sys.time()), "] SnpToGene (marg = 10000, PAPER PARAM) ...\n")
Pagwas$snp_gene_df <- scPagwas::SnpToGene(
  gwas_data = Pagwas$gwas_data,
  block_annotation = block_annotation,
  marg = 10000              # PAPER PARAM
)

Pagwas$rawPathway_list <- Genes_by_pathway_kegg

cat("Saving to", OUT_RDATA, "\n")
save(Pagwas, file = OUT_RDATA, compress = TRUE)
cat("DONE\n")
