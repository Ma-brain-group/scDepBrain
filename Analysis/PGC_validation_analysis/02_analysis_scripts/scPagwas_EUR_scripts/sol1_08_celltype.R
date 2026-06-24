# ============================================================
# Solution 1 — Step 1.8 (NEW): Cell-type-level regression + bootstrap
#
# Matches the user's STANDARD scPagwas pipeline:
#   3.7 Pagwas_perform_regression  (cell-type beta estimates)
#   Boot_evaluate (iters=200, part=0.5) (cell-type bootstrap p-values)
#
# Inputs (reused from earlier steps):
#   shared pagwas_preprocessed.RData (data_mat + pca_scCell_mat +
#     merge_scexpr + pca_cell_df + Celltype_anno + Pathway_list +
#     rawPathway_list + VariableFeatures)
#   gwas_pagwas.RData (gwas_data + snp_gene_df + rawPathway_list)
#
# Outputs (NEW, complementing earlier v2 outputs):
#   Pagwas_celltype_results.RData  (lm_results + bootstrap_results)
#   cell_type_pvalues_celltype_method.csv  (one row per cell type)
#
# Memory: peak ~150-300 GB (celltype mode is lighter than singlecell)
# Runtime: Link_pathway celltype-only ~1-3h + bootstrap ~10-20 min
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
OUT    <- "/mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_scPagwas_MDD2025_EUR_20260610"
BACKINGPATH <- file.path(OUT, "sol1_celltype_temp")
dir.create(BACKINGPATH, showWarnings = FALSE, recursive = TRUE)

# === Load shared preprocessed Pagwas (full data_mat + pca_scCell_mat + ...) ===
cat("[", as.character(Sys.time()), "] Loading shared pagwas_preprocessed.RData ...\n")
load(file.path(SHARED, "pagwas_preprocessed.RData"))   # brings Pagwas + Single_data
rm(Single_data); gc()   # not needed for cell-type analysis
Pagwas_full <- Pagwas
rm(Pagwas)

# === Load GWAS Pagwas ===
cat("[", as.character(Sys.time()), "] Loading gwas_pagwas.RData ...\n")
load(file.path(OUT, "gwas_pagwas.RData"))   # brings 'Pagwas' (gwas-side)
Pagwas <- c(Pagwas, Pagwas_full)
rm(Pagwas_full); gc()
cat("Pagwas merged. Names:", paste(names(Pagwas), collapse = ", "), "\n")

# === Pathway_annotation_input ===
cat("[", as.character(Sys.time()), "] Pathway_annotation_input ...\n")
Pagwas <- scPagwas::Pathway_annotation_input(
  Pagwas = Pagwas,
  block_annotation = block_annotation
)

# === Link_pathway_blocks_gwas in CELLTYPE-ONLY mode (no singlecell -> fast) ===
cat("[", as.character(Sys.time()), "] Link_pathway_blocks_gwas (celltype=TRUE, singlecell=FALSE) ...\n")
Pagwas <- scPagwas::Link_pathway_blocks_gwas(
  Pagwas      = Pagwas,
  chrom_ld    = chrom_ld,
  singlecell  = FALSE,
  celltype    = TRUE,
  backingpath = BACKINGPATH,
  n.cores     = 1
)
cat("Pathway_ld_gwas_data length:", length(Pagwas$Pathway_ld_gwas_data), "\n")

# === Pagwas_perform_regression for cell types ===
cat("[", as.character(Sys.time()), "] Pagwas_perform_regression ...\n")
Pagwas$lm_results <- scPagwas::Pagwas_perform_regression(
  Pathway_ld_gwas_data = Pagwas$Pathway_ld_gwas_data
)
cat("\n=== Cell-type regression beta (lm_results) ===\n")
print(Pagwas$lm_results)

# === Boot_evaluate bootstrap p-values (200 iters, 50% subsample) ===
cat("\n[", as.character(Sys.time()), "] Boot_evaluate (200 iters, part=0.5) ...\n")
Pagwas <- scPagwas::Boot_evaluate(Pagwas, bootstrap_iters = 200, part = 0.5)

# Free large object before save
Pagwas$Pathway_ld_gwas_data <- NULL

# === Save and report ===
cat("\n=== Bootstrap results ===\n")
print(Pagwas$bootstrap_results)

ct_df <- if (!is.null(Pagwas$bootstrap_results)) {
  Pagwas$bootstrap_results
} else {
  data.frame(annotation = names(Pagwas$lm_results$bv),
             bv = Pagwas$lm_results$bv)
}
write.csv(ct_df,
          file = file.path(OUT, "cell_type_pvalues_celltype_method.csv"),
          row.names = FALSE)
cat("Saved cell_type_pvalues_celltype_method.csv\n")

# Save the full result list for downstream / debugging
save(Pagwas, file = file.path(OUT, "Pagwas_celltype_results.RData"), compress = TRUE)
cat("Saved Pagwas_celltype_results.RData\n")

# Cleanup backing files
unlink(BACKINGPATH, recursive = TRUE)

cat("\n=== DONE (celltype standard pipeline) ===\n")
