# ============================================================
# HW paper-params: Step 1.8 celltype-only analysis
# Pagwas_perform_regression + Boot_evaluate matching paper code
# ============================================================
suppressPackageStartupMessages({
  library(Seurat); library(scPagwas); library(Matrix)
})

if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
  RhpcBLASctl::blas_set_num_threads(1)
}
options(default.nproc.blas = 1)
options(bigstatsr.check.parallel.blas = FALSE)
suppressMessages({
  utils::assignInNamespace(".Deprecate",
                           function(...) invisible(NULL),
                           ns = "SeuratObject")
})

SHARED <- "/mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_scPagwas_shared_20260612"
OUT    <- "/mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_scPagwas_DIV_paperparams_20260615"
BACKINGPATH <- file.path(OUT, "sol1_celltype_temp")
dir.create(BACKINGPATH, showWarnings = FALSE, recursive = TRUE)

cat("[", as.character(Sys.time()), "] Loading shared pagwas_preprocessed.RData ...\n")
load(file.path(SHARED, "pagwas_preprocessed.RData"))
rm(Single_data); gc()
Pagwas_full <- Pagwas
rm(Pagwas)

cat("[", as.character(Sys.time()), "] Loading paper-params gwas_pagwas.RData ...\n")
load(file.path(OUT, "gwas_pagwas.RData"))
Pagwas <- c(Pagwas, Pagwas_full)
rm(Pagwas_full); gc()

cat("[", as.character(Sys.time()), "] Pathway_annotation_input ...\n")
Pagwas <- scPagwas::Pathway_annotation_input(
  Pagwas = Pagwas,
  block_annotation = block_annotation
)

cat("[", as.character(Sys.time()), "] Link_pathway_blocks_gwas (celltype=TRUE) ...\n")
Pagwas <- scPagwas::Link_pathway_blocks_gwas(
  Pagwas      = Pagwas,
  chrom_ld    = chrom_ld,
  singlecell  = FALSE,
  celltype    = TRUE,
  backingpath = BACKINGPATH,
  n.cores     = 1
)

cat("[", as.character(Sys.time()), "] Pagwas_perform_regression ...\n")
Pagwas$lm_results <- scPagwas::Pagwas_perform_regression(
  Pathway_ld_gwas_data = Pagwas$Pathway_ld_gwas_data
)
cat("\n=== lm_results ===\n"); print(Pagwas$lm_results)

cat("\n[", as.character(Sys.time()), "] Boot_evaluate (200 iters, part=0.5) ...\n")
Pagwas <- scPagwas::Boot_evaluate(Pagwas, bootstrap_iters = 200, part = 0.5)
Pagwas$Pathway_ld_gwas_data <- NULL

cat("\n=== bootstrap_results ===\n")
print(Pagwas$bootstrap_results)

write.csv(Pagwas$bootstrap_results,
          file = file.path(OUT, "cell_type_pvalues_celltype_method_paperparams.csv"),
          row.names = FALSE)

save(Pagwas, file = file.path(OUT, "Pagwas_celltype_results_paperparams.RData"), compress = TRUE)

unlink(BACKINGPATH, recursive = TRUE)
cat("\n=== DONE (celltype paper-params) ===\n")
