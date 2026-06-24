# ============================================================
# Run scPagwas: PGC MDD 2025 EUR ancestry GWAS × MDD single-cell rds.
# Uses same compat patches as Howard 2019 validation run:
#   - SeuratObject:::.Deprecate -> no-op  (allow GetAssayData(slot=))
#   - Assay5 -> v3 Assay downgrade  (if needed)
#   - Drop SCT  (save memory)
#   - Relative output.dirs  (avoid malformed "./" + abs concat)
#   - Pre-create scPagwas_cache & temp dirs
# ============================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(scPagwas)
})

# Avoid bigparallelr "Two levels of parallelism" (see MDD_20260531 script)
if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
  RhpcBLASctl::blas_set_num_threads(1)
}
options(default.nproc.blas = 1)
options(bigstatsr.check.parallel.blas = FALSE)

# Compat patch (see 04_run_scPagwas.R in MDD_20260531 for full rationale)
suppressMessages({
  utils::assignInNamespace(".Deprecate",
                           function(...) invisible(NULL),
                           ns = "SeuratObject")
})
cat("Patched SeuratObject:::.Deprecate -> no-op.\n")

# ---- Paths ----
RDS_IN  <- "/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/integrated_heathy_brain_cells_short_read/MDD_singlecell_data_reannotation.rds"
OUT     <- "/mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_scPagwas_MDD2025_TransAnc_20260610"
GWAS    <- file.path(OUT, "pgc-mdd2025_div_scPagwas_input.txt")
RDS_OUT <- file.path(OUT, "MDD_scPagwas_PGC2025_TransAnc_v2026.rds")
PREFIX  <- "MDD_PGC2025_TransAnc_v2026"

setwd(OUT)
OUT_DIR <- "scPagwas_outputs"
dir.create(OUT_DIR,                              showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUT_DIR, "scPagwas_cache"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUT_DIR, "temp"),           showWarnings = FALSE, recursive = TRUE)

# ---- Load Seurat ----
cat("[", as.character(Sys.time()), "] Loading 24GB rds ...\n")
t0 <- Sys.time()
obj <- readRDS(RDS_IN)
cat(sprintf("Loaded in %.1f min\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))

# ---- Prepare object ----
DefaultAssay(obj) <- "RNA"

# Use MAJOR cell type column `anno` (8 categories) — matches published figure.
Idents(obj) <- obj$anno
cat("Idents() set from $anno (major types):\n")
print(table(Idents(obj)))

if (inherits(obj[["RNA"]], "Assay5")) {
  cat("Downgrading Assay5 -> v3 Assay ...\n")
  obj <- JoinLayers(obj, assay = "RNA")
  obj[["RNA"]] <- as(obj[["RNA"]], "Assay")
}
if ("SCT" %in% Assays(obj)) {
  cat("Dropping SCT assay\n")
  obj[["SCT"]] <- NULL
  gc()
}

# ---- Preserve OLD scPagwas (Howard 2019) columns by renaming ----
old_cols <- c("scPagwas.gPAS.score", "scPagwas.TRS.Score", "Random_Correct_BG_adjp")
for (col in old_cols) {
  if (col %in% colnames(obj@meta.data)) {
    new_name <- paste0(col, "_Howard2019_OLD")
    obj@meta.data[[new_name]] <- obj@meta.data[[col]]
    obj@meta.data[[col]] <- NULL
    cat(sprintf("Renamed %s -> %s (preserved)\n", col, new_name))
  }
}

# ---- Run scPagwas ----
cat("[", as.character(Sys.time()), "] Starting scPagwas_main2 (PGC 2025 Trans-ancestry) ...\n")
t1 <- Sys.time()
Pagwas <- scPagwas_main2(
  Pagwas           = NULL,
  gwas_data        = GWAS,
  Single_data      = obj,
  output.prefix    = PREFIX,
  output.dirs      = OUT_DIR,
  Pathway_list     = Genes_by_pathway_kegg,
  assay            = "RNA",
  block_annotation = block_annotation,
  chrom_ld         = chrom_ld,
  singlecell       = TRUE,
  celltype         = TRUE,
  iters_singlecell = 100,
  iters_celltype   = 200,
  marg             = 20000,        # ±20 kb (manuscript)
  maf_filter       = 0.01,
  min_clustercells = 10,
  min.pathway.size = 5,
  max.pathway.size = 300,
  n_topgenes       = 1000,
  n.cores          = 16,
  remove_outlier   = TRUE,
  seurat_return    = TRUE
)
cat(sprintf("scPagwas_main2 finished in %.1f min\n",
            as.numeric(difftime(Sys.time(), t1, units = "mins"))))

# ---- Save result Seurat ----
cat("Saving result to", RDS_OUT, "...\n")
saveRDS(Pagwas, file = RDS_OUT)

# ---- Cell-type summary ----
md <- Pagwas@meta.data
cat("\n=== Cell-type level TRS (mean per cell type) ===\n")
if ("scPagwas.TRS.Score" %in% colnames(md)) {
  ct_summary <- aggregate(scPagwas.TRS.Score ~ Idents(Pagwas),
                          data = md, FUN = function(x) c(mean = mean(x), n = length(x)))
  print(ct_summary)
}

cat("\n=== DONE ===\n")
