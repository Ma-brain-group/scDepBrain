# ============================================================
# Cell-type gene sets v2: FIX two bugs
#   Bug 1: use AverageExpression (mean per ct), not AggregateExpression (sum)
#   Bug 2: also output --gene-covar file (continuous specificity, Bryois 2020 / MAGMA_CellTyping standard)
#
# Outputs in celltype_genesets/:
#   anno_top10_v2.txt              — set-annot, top 10% by specificity, 8 major types
#   anno_gene_covar_v2.txt         — gene-covar, continuous specificity, 8 major types
#   final_anno_top10_v2.txt        — set-annot, 20 subtypes
#   final_anno_gene_covar_v2.txt   — gene-covar, 20 subtypes
#
# Specificity formula (Bryois 2020):
#   spec[g, ct] = avg_expr[g, ct] / sum_over_ct(avg_expr[g, ct])
#
# gene_covar.txt format expected by MAGMA:
#   GENE celltype1 celltype2 ... celltypeN
#   ENTREZ1 0.3 0.1 ...
# (space-delimited, header line with GENE + cell type names)
# ============================================================
suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(dplyr)
})
suppressMessages({
  utils::assignInNamespace(".Deprecate", function(...) invisible(NULL),
                           ns = "SeuratObject")
})

RDS    <- "/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/integrated_heathy_brain_cells_short_read/MDD_singlecell_data_reannotation.rds"
OUT    <- "/mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_MAGMA_PGC2025_20260617/celltype_genesets"
GENELOC <- "/mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_MAGMA_PGC2025_20260617/NCBI37.3.gene.loc.extendedMHCexcluded"

dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# Symbol -> Entrez via gene-loc
gl <- read.table(GENELOC, header = FALSE, stringsAsFactors = FALSE,
                 col.names = c("entrez", "chr", "start", "end", "strand", "symbol"))
sym2ent <- setNames(as.character(gl$entrez), gl$symbol)

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

make_set_files <- function(group_var, prefix, label) {
  cat("\n========== ", label, " ==========\n", sep = "")
  Idents(obj) <- obj[[group_var]][[1]]
  cts <- levels(droplevels(factor(Idents(obj))))
  cat("Cell types (", length(cts), "):\n  ", paste(cts, collapse = ", "), "\n")

  # *** FIX BUG 1: AverageExpression (mean), not AggregateExpression (sum) ***
  cat("[", as.character(Sys.time()), "] AverageExpression — mean per cell type ...\n")
  avg <- AverageExpression(obj, assays = "RNA", group.by = group_var,
                           return.seurat = FALSE)$RNA
  cat("avg dim:", paste(dim(avg), collapse = " x "), "\n")
  cat("avg colnames:", paste(colnames(avg), collapse = ", "), "\n")

  # Filter to genes mappable to Entrez (extended-MHC-excluded)
  keep_sym <- intersect(rownames(avg), names(sym2ent))
  cat("Genes mappable to Entrez:", length(keep_sym), "/", nrow(avg), "\n")
  avg <- avg[keep_sym, , drop = FALSE]
  # Drop genes with 0 sum across all cell types
  rowtot <- rowSums(avg)
  avg <- avg[rowtot > 0, , drop = FALSE]
  cat("After zero-sum filter:", nrow(avg), "\n")

  # Specificity
  rowtot <- rowSums(avg)
  spec <- sweep(avg, 1, rowtot, "/")

  # Convert rownames to Entrez
  ent <- sym2ent[rownames(spec)]
  # Dedup: if multiple symbols map to same entrez, aggregate by mean
  spec_df <- as.data.frame(as.matrix(spec))
  spec_df$entrez <- ent
  spec_df <- spec_df %>%
    group_by(entrez) %>%
    summarise(across(everything(), mean), .groups = "drop")
  ent_vec <- spec_df$entrez
  spec_mat <- as.matrix(spec_df[, -1])
  rownames(spec_mat) <- ent_vec
  cat("After Entrez dedup:", nrow(spec_mat), "unique Entrez genes\n")

  # ----- Output 1: gene_covar (CONTINUOUS, Bryois MAGMA_CellTyping standard) -----
  covar_file <- file.path(OUT, paste0(prefix, "_gene_covar_v2.txt"))
  # MAGMA wants: GENE celltype1 celltype2 ...  (whitespace; first column "GENE")
  # Replace spaces/slashes in cell type names for header
  ct_clean <- gsub("[ /]+", "_", colnames(spec_mat))
  colnames(spec_mat) <- ct_clean
  out_df <- data.frame(GENE = ent_vec, spec_mat, check.names = FALSE)
  write.table(out_df, file = covar_file, quote = FALSE, sep = " ",
              row.names = FALSE)
  cat("Wrote", nrow(out_df), "rows x", ncol(spec_mat), "cell types to", covar_file, "\n")

  # ----- Output 2: top 10% set-annot (for comparison/sanity) -----
  set_file <- file.path(OUT, paste0(prefix, "_top10_v2.txt"))
  n_per <- floor(nrow(spec_mat) * 0.10)
  cat("Top 10% genes/ct:", n_per, "\n")
  lines <- character()
  for (ct in colnames(spec_mat)) {
    top_idx <- order(spec_mat[, ct], decreasing = TRUE)[1:n_per]
    top_ent <- rownames(spec_mat)[top_idx]
    lines <- c(lines, paste(c(ct, top_ent), collapse = " "))
  }
  writeLines(lines, set_file)
  cat("Wrote", length(lines), "sets to", set_file, "\n")
}

make_set_files("anno",       "anno",       "anno (8 major)")
make_set_files("final_anno", "final_anno", "final_anno (20 subtypes)")

cat("\n=== DONE (v2) ===\n")
