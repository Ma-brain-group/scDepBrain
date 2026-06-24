# ============================================================
# Extract cell-type-specific gene sets for MAGMA enrichment
# Bryois et al. 2020 Nat Genet method:
#   1. Average expression per cell type
#   2. Specificity = avg_expr_in_type / sum(avg_expr_across_types)
#   3. Top 10% genes by specificity per cell type
#   4. Convert gene symbol -> Entrez (match NCBI37.3.gene.loc)
# Output 2 set-annot files in MAGMA format:
#   anno_top10.txt        — 8 major cell types
#   final_anno_top10.txt  — 20 subtypes
# Format: <set_name> <gene1> <gene2> ... <geneN>   (space-delim, one set per row)
# ============================================================
suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(dplyr); library(tibble)
})

suppressMessages({
  utils::assignInNamespace(".Deprecate", function(...) invisible(NULL),
                           ns = "SeuratObject")
})

RDS    <- "/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/integrated_heathy_brain_cells_short_read/MDD_singlecell_data_reannotation.rds"
OUT    <- "/mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_MAGMA_PGC2025_20260617/celltype_genesets"
GENELOC <- "/mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_MAGMA_PGC2025_20260617/NCBI37.3.gene.loc.extendedMHCexcluded"

# Read gene-loc to build symbol -> entrez mapping
cat("[", as.character(Sys.time()), "] Reading gene-loc ...\n")
gl <- read.table(GENELOC, header = FALSE, stringsAsFactors = FALSE,
                 col.names = c("entrez", "chr", "start", "end", "strand", "symbol"))
cat("gene-loc rows:", nrow(gl), "\n")
sym2ent <- setNames(as.character(gl$entrez), gl$symbol)

cat("[", as.character(Sys.time()), "] Loading 24 GB rds ...\n")
t0 <- Sys.time()
obj <- readRDS(RDS)
DefaultAssay(obj) <- "RNA"
if (inherits(obj[["RNA"]], "Assay5")) {
  obj <- JoinLayers(obj, assay = "RNA")
  obj[["RNA"]] <- as(obj[["RNA"]], "Assay")
}
cat(sprintf("Loaded in %.1f min. dim: %d x %d\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins")),
            nrow(obj), ncol(obj)))

# Drop SCT to save memory
if ("SCT" %in% Assays(obj)) {
  obj[["SCT"]] <- NULL; gc()
}

# Use RNA data slot (log-normalized) for specificity
make_set_file <- function(group_var, out_file, label) {
  cat("\n========== ", label, " ==========\n", sep = "")
  Idents(obj) <- obj[[group_var]][[1]]
  cell_types <- levels(droplevels(factor(Idents(obj))))
  cat("Cell types (", length(cell_types), "):\n  ", paste(cell_types, collapse = ", "), "\n")

  # Average expression per cell type (rows = genes, cols = cell types)
  cat("[", as.character(Sys.time()), "] AggregateExpression (", label, ") ...\n")
  avg <- AggregateExpression(obj, assays = "RNA", group.by = group_var,
                             return.seurat = FALSE)$RNA
  # AggregateExpression returns matrix; row = gene symbols, col = cell types
  # Sometimes column names get sanitized (spaces/special chars -> dots)
  cat("avg dim:", paste(dim(avg), collapse = " x "), "\n")
  cat("avg colnames:", paste(colnames(avg), collapse = ", "), "\n")

  # Normalize to per-cell-type total to compute specificity
  # specificity[gene, ct] = expr[gene, ct] / sum_over_ct(expr[gene, ct])
  rowtot <- rowSums(avg)
  # avoid 0/0 — drop genes with zero total expression
  keep <- rowtot > 0
  avg2 <- avg[keep, , drop = FALSE]
  rowtot2 <- rowtot[keep]
  spec <- sweep(avg2, 1, rowtot2, "/")
  cat("specificity matrix dim:", paste(dim(spec), collapse = " x "), "\n")

  # Per cell type, top 10% genes by specificity, restricted to genes with
  # an Entrez mapping in the (extended-MHC-excluded) gene-loc
  n_per_ct <- floor(nrow(spec) * 0.10)
  cat("Top 10% genes per cell type: keeping", n_per_ct, "/ ", nrow(spec), "\n")

  lines <- character()
  for (ct in colnames(spec)) {
    top_idx <- order(spec[, ct], decreasing = TRUE)[1:n_per_ct]
    top_sym <- rownames(spec)[top_idx]
    top_ent <- sym2ent[top_sym]
    top_ent <- unique(top_ent[!is.na(top_ent) & top_ent != ""])
    cat(sprintf("  %s : %d top symbols, %d unique Entrez after mapping\n",
                ct, length(top_sym), length(top_ent)))
    set_name <- gsub("[ /]+", "_", ct)
    lines <- c(lines, paste(c(set_name, top_ent), collapse = " "))
  }
  writeLines(lines, out_file)
  cat("Wrote", length(lines), "sets to", out_file, "\n")
}

dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
make_set_file("anno",       file.path(OUT, "anno_top10.txt"),       "anno (8 major types)")
make_set_file("final_anno", file.path(OUT, "final_anno_top10.txt"), "final_anno (20 subtypes)")

cat("\n=== DONE ===\n")
