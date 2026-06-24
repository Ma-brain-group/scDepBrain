# =================================================================
# EWCE-style bootstrap: 295 PGC genes vs cell-type specificity matrix
# Standard method (Skene 2016) for testing if a pre-defined gene list
# is enriched for cell-type-specific expression.
# Conceptually equivalent to MAGMA CellTyping, but for a fixed gene list.
# =================================================================
suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(dplyr)
})
utils::assignInNamespace(".Deprecate", function(...) invisible(NULL),
                         ns = "SeuratObject")

RDS <- "/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/integrated_heathy_brain_cells_short_read/MDD_singlecell_data_reannotation.rds"
PGC_FILE <- "/mnt/isilon/gandal_lab/mayl/07_scDepBrain/PGC_308_high_confidence_genes.txt"
OUT <- "/mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_PGC295_geneset_20260618"
DAT <- file.path(OUT, "data")
N_BOOT <- 10000          # bootstrap replicates
set.seed(42)

pgc <- readLines(PGC_FILE); pgc <- pgc[nchar(trimws(pgc)) > 0]; pgc <- unique(pgc)
cat("PGC genes:", length(pgc), "\n")

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

# ------------- specificity matrix (Bryois: avg / sum) -------------
compute_spec <- function(obj, group_var) {
  avg <- AverageExpression(obj, assays = "RNA", group.by = group_var,
                           return.seurat = FALSE)$RNA
  avg <- avg[rowSums(avg) > 0, , drop = FALSE]
  spec <- sweep(avg, 1, rowSums(avg), "/")
  spec
}

# ------------- EWCE bootstrap -------------
ewce_bootstrap <- function(spec, gene_list, n_boot = 10000) {
  bg_genes <- rownames(spec)
  pgc_in <- intersect(gene_list, bg_genes)
  k <- length(pgc_in)
  cat("  Gene list in spec matrix:", k, "/", length(gene_list), "\n")

  # Observed mean specificity for PGC genes
  obs <- colMeans(spec[pgc_in, , drop = FALSE])

  # Bootstrap null: random samples of same size from background
  null <- matrix(0, nrow = n_boot, ncol = ncol(spec))
  colnames(null) <- colnames(spec)
  for (i in seq_len(n_boot)) {
    rand <- sample(bg_genes, k, replace = FALSE)
    null[i, ] <- colMeans(spec[rand, , drop = FALSE])
  }

  # One-sided P (observed > null) + fold-change vs null mean
  p_val <- sapply(seq_along(obs), function(j)
    (sum(null[, j] >= obs[j]) + 1) / (n_boot + 1))
  fc <- obs / colMeans(null)
  z <- (obs - colMeans(null)) / apply(null, 2, sd)

  data.frame(
    cell_type = colnames(spec),
    n_PGC = k,
    obs_mean_spec = round(obs, 5),
    null_mean_spec = round(colMeans(null), 5),
    fold_change = round(fc, 2),
    bootstrap_z = round(z, 2),
    P = signif(p_val, 3)
  )
}

# ----- MAJOR (8 cell types) -----
cat("\n[", as.character(Sys.time()), "] Major specificity ...\n")
spec_major <- compute_spec(obj, "anno")
cat("Major spec dim:", paste(dim(spec_major), collapse = " x "), "\n")
res_major <- ewce_bootstrap(spec_major, pgc, N_BOOT)
res_major$FDR <- signif(p.adjust(res_major$P, "BH"), 3)
res_major <- res_major %>% arrange(P)

cat("\n=== MAJOR (8) EWCE bootstrap ===\n")
print(res_major, row.names = FALSE)
write.csv(res_major, file.path(DAT, "EWCE_bootstrap_major.csv"), row.names = FALSE)

# ----- SUBTYPE (20) -----
cat("\n[", as.character(Sys.time()), "] Subtype specificity ...\n")
spec_sub <- compute_spec(obj, "final_anno")
cat("Subtype spec dim:", paste(dim(spec_sub), collapse = " x "), "\n")
res_sub <- ewce_bootstrap(spec_sub, pgc, N_BOOT)
res_sub$FDR <- signif(p.adjust(res_sub$P, "BH"), 3)
res_sub <- res_sub %>% arrange(P)

cat("\n=== SUBTYPE (20) EWCE bootstrap ===\n")
print(res_sub, row.names = FALSE)
write.csv(res_sub, file.path(DAT, "EWCE_bootstrap_subtype.csv"), row.names = FALSE)

cat("\n[", as.character(Sys.time()), "] DONE\n")
