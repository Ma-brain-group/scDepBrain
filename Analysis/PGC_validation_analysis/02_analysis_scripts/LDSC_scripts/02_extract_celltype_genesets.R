# Extract cell-type gene sets for LDSC-SEG:
#   For each cell type (anno + final_anno), top 10% genes by specificity (Bryois)
#   Output one .geneset.txt per cell type (gene symbols, one per line)
suppressPackageStartupMessages({
  library(Seurat); library(Matrix)
})
utils::assignInNamespace(".Deprecate", function(...) invisible(NULL),
                         ns = "SeuratObject")

RDS <- "/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/integrated_heathy_brain_cells_short_read/MDD_singlecell_data_reannotation.rds"
OUT <- "/mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_LDSC_PGC2025_20260617/celltype_genesets"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

cat("[", as.character(Sys.time()), "] Loading 24GB rds ...\n")
t0 <- Sys.time()
obj <- readRDS(RDS)
DefaultAssay(obj) <- "RNA"
if (inherits(obj[["RNA"]], "Assay5")) {
  obj <- JoinLayers(obj, assay = "RNA")
  obj[["RNA"]] <- as(obj[["RNA"]], "Assay")
}
if ("SCT" %in% Assays(obj)) { obj[["SCT"]] <- NULL; gc() }
cat(sprintf("Loaded in %.1f min\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))

make_sets <- function(group_var, prefix) {
  cat("\n=== ", prefix, " (group:", group_var, ") ===\n")
  avg <- AverageExpression(obj, assays = "RNA", group.by = group_var,
                           return.seurat = FALSE)$RNA
  cat("avg dim:", paste(dim(avg), collapse = " x "), "\n")

  # Filter zero-sum
  keep <- rowSums(avg) > 0
  avg <- avg[keep, , drop = FALSE]
  spec <- sweep(avg, 1, rowSums(avg), "/")

  n_per <- floor(nrow(spec) * 0.10)
  cat("Top 10% genes per cell type:", n_per, "\n")
  cell_types <- colnames(spec)
  for (ct in cell_types) {
    top_idx <- order(spec[, ct], decreasing = TRUE)[1:n_per]
    top_genes <- rownames(spec)[top_idx]
    clean_ct <- gsub("[ /]+", "_", ct)
    fname <- file.path(OUT, paste0(prefix, "_", clean_ct, ".geneset.txt"))
    writeLines(top_genes, fname)
    cat("  ", clean_ct, ":", length(top_genes), "genes ->", basename(fname), "\n")
  }
}

make_sets("anno",       "major")
make_sets("final_anno", "subtype")

cat("\nAll gene sets in:", OUT, "\n")
cat("\n=== DONE ===\n")
