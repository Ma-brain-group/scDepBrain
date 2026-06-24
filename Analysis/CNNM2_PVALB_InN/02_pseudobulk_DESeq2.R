# =================================================================
# Pseudobulk DESeq2 DE — CNNM2 MDD vs Control
# Standard donor-level approach per Squair et al. 2021 (Nat Commun)
# Aggregate raw counts per (donor × cell_type), then DESeq2
# Compares with previous per-cell Wilcoxon (which is pseudoreplication)
#
# Levels:
#   (1) 8 broad cell types from 2023_NC.rds
#   (2) Inhibitory subtypes from NC_inhibitory_neurons_reannotation.rds
#   (3) Excitatory subtypes from NC_Excitatory_neuron_subtypes.rds
# =================================================================
suppressPackageStartupMessages({
  library(Seurat); library(DESeq2); library(Matrix); library(dplyr)
  library(ggplot2); library(ggpubr); library(data.table)
})
utils::assignInNamespace(".Deprecate", function(...) invisible(NULL),
                         ns = "SeuratObject")

DATA_DIR <- "/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/replication_4_MDD_brain_cells"
OUT <- "/mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_CNNM2_replication_20260619"
FIG <- file.path(OUT, "figures"); DAT <- file.path(OUT, "data")

GENE <- "CNNM2"

# =================================================================
# Helper: pseudobulk + DESeq2 per cell type
# =================================================================
pseudobulk_de <- function(obj, ct_col, min_cells = 10, covar_cols = NULL) {
  # Get raw counts
  counts <- GetAssayData(obj, assay = "RNA", slot = "counts")
  md <- obj@meta.data
  md$cell_id <- rownames(md)

  # Define grouping: donor × cell_type
  if (!"sample" %in% colnames(md)) stop("'sample' column missing")
  if (!"phenotype" %in% colnames(md)) stop("'phenotype' column missing")
  if (!ct_col %in% colnames(md)) stop(ct_col, " column missing")
  md$grp <- paste(md$sample, md[[ct_col]], sep = "__")

  # Filter to phenotype == Case or Control
  keep <- md$phenotype %in% c("Case","Control")
  md <- md[keep, ]; counts <- counts[, rownames(md)]

  # Filter (sample × cell_type) with at least min_cells
  cell_counts <- table(md$grp)
  good_grp <- names(cell_counts)[cell_counts >= min_cells]
  md <- md[md$grp %in% good_grp, ]; counts <- counts[, rownames(md)]
  cat("Cells after filter:", ncol(counts), "  Groups:", length(good_grp), "\n")

  # Aggregate counts: sum within each grp (pseudobulk) via sparse matrix multiplication
  grp_factor <- factor(md$grp)
  # group indicator matrix: n_cells x n_groups
  ind_mat <- Matrix::sparse.model.matrix(~ 0 + grp_factor)
  colnames(ind_mat) <- gsub("^grp_factor", "", colnames(ind_mat))
  # counts (genes x cells) %*% ind (cells x groups) = genes x groups
  pseudo <- as.matrix(counts %*% ind_mat)
  colnames(pseudo) <- colnames(ind_mat)
  cat("Pseudobulk matrix:", paste(dim(pseudo), collapse = " x "), "\n")

  # Sample metadata at pseudobulk level
  smd <- md %>% select(grp, sample, phenotype, any_of(c("sex","batch","Chemistry","age")), all_of(ct_col)) %>%
    distinct(grp, .keep_all = TRUE)
  rownames(smd) <- smd$grp
  smd <- smd[colnames(pseudo), , drop = FALSE]
  smd$phenotype <- factor(smd$phenotype, levels = c("Control","Case"))

  # Test per cell type
  cts <- unique(smd[[ct_col]])
  res_list <- lapply(cts, function(c) {
    keep <- smd[[ct_col]] == c
    s <- smd[keep, , drop = FALSE]
    if (length(unique(s$phenotype)) < 2) return(NULL)
    if (min(table(s$phenotype)) < 3) {
      cat("Skip", c, ": min donor < 3\n"); return(NULL)
    }
    p <- pseudo[, keep, drop = FALSE]
    # Filter low-count genes within this cell type
    p <- p[rowSums(p) >= 10, , drop = FALSE]

    # Build design (covariates if available + balanced)
    design <- ~ phenotype
    used_covars <- c()
    if ("sex" %in% colnames(s) && length(unique(s$sex)) > 1 &&
        all(table(s$sex, s$phenotype) > 0)) {
      design <- update(design, ~ sex + .)
      used_covars <- c(used_covars, "sex")
    }

    dds <- tryCatch({
      DESeqDataSetFromMatrix(countData = p, colData = s, design = design)
    }, error = function(e) {cat("DESeq error:", c, conditionMessage(e), "\n"); NULL})
    if (is.null(dds)) return(NULL)
    dds <- DESeq(dds, quiet = TRUE)
    res <- results(dds, contrast = c("phenotype","Case","Control"))
    cnnm2 <- if (GENE %in% rownames(res)) as.data.frame(res[GENE, ]) else NULL
    if (!is.null(cnnm2)) {
      cnnm2$cell_type <- c
      cnnm2$n_donors_case <- sum(s$phenotype == "Case")
      cnnm2$n_donors_ctrl <- sum(s$phenotype == "Control")
      cnnm2$covariates <- paste(used_covars, collapse=",")
    }
    cnnm2
  })
  res_df <- do.call(rbind, res_list)
  return(res_df)
}

# =================================================================
# Step 1: All broad cell types from 2023_NC.rds
# =================================================================
cat("\n[", as.character(Sys.time()), "] Loading 2023_NC.rds ...\n")
obj <- readRDS(file.path(DATA_DIR, "2023_NC.rds"))
DefaultAssay(obj) <- "RNA"
if (inherits(obj[["RNA"]], "Assay5")) {
  obj <- JoinLayers(obj, assay = "RNA")
  obj[["RNA"]] <- as(obj[["RNA"]], "Assay")
}
cat("Cells:", ncol(obj),
    "  Donors (samples):", length(unique(obj$sample)),
    "  Phenotype levels:", paste(unique(obj$phenotype), collapse = ", "), "\n")

# Confirm n_donors per phenotype
d_pheno <- obj@meta.data %>% select(sample, phenotype) %>% distinct()
cat("Donor count per phenotype:\n")
print(table(d_pheno$phenotype))

# DE for broad cell types
res_broad <- pseudobulk_de(obj, ct_col = "broad_cell_type")
res_broad$level <- "broad_cell_type"
write.csv(res_broad, file.path(DAT, "pseudobulk_CNNM2_broad.csv"), row.names = FALSE)
rm(obj); gc()

# =================================================================
# Step 2: Inhibitory subtypes
# =================================================================
cat("\n[", as.character(Sys.time()), "] Loading NC_inhibitory ...\n")
obj_in <- readRDS(file.path(DATA_DIR, "NC_inhibitory_neurons_reannotation.rds"))
DefaultAssay(obj_in) <- "RNA"
if (inherits(obj_in[["RNA"]], "Assay5")) {
  obj_in <- JoinLayers(obj_in, assay = "RNA")
  obj_in[["RNA"]] <- as(obj_in[["RNA"]], "Assay")
}
cat("Inhibitory cells:", ncol(obj_in), "\n")
res_in <- pseudobulk_de(obj_in, ct_col = "RNA_snn_res.0.1_cell_type")
res_in$level <- "Inhibitory_subtype"
write.csv(res_in, file.path(DAT, "pseudobulk_CNNM2_inhibitory.csv"), row.names = FALSE)
rm(obj_in); gc()

# =================================================================
# Step 3: Excitatory subtypes
# =================================================================
cat("\n[", as.character(Sys.time()), "] Loading NC_Excitatory ...\n")
obj_ex <- readRDS(file.path(DATA_DIR, "NC_Excitatory_neuron_subtypes.rds"))
DefaultAssay(obj_ex) <- "RNA"
if (inherits(obj_ex[["RNA"]], "Assay5")) {
  obj_ex <- JoinLayers(obj_ex, assay = "RNA")
  obj_ex[["RNA"]] <- as(obj_ex[["RNA"]], "Assay")
}
ex_subtype_col <- grep("Excitatory.*subtype|subtype|cell_type",
                       colnames(obj_ex@meta.data),
                       value = TRUE, ignore.case = TRUE)
use_col <- ex_subtype_col[1]
for (c in ex_subtype_col) {
  v <- as.character(obj_ex@meta.data[[c]])
  if (any(grepl("^Ex-", v))) { use_col <- c; break }
}
cat("Using col:", use_col, "\n")
cat("Excitatory cells:", ncol(obj_ex), "\n")
res_ex <- pseudobulk_de(obj_ex, ct_col = use_col)
res_ex$level <- "Excitatory_subtype"
write.csv(res_ex, file.path(DAT, "pseudobulk_CNNM2_excitatory.csv"), row.names = FALSE)
rm(obj_ex); gc()

# =================================================================
# Combine + plot
# =================================================================
all_res <- bind_rows(res_broad, res_in, res_ex)
all_res$FDR <- p.adjust(all_res$pvalue, method = "BH")
all_res <- all_res %>%
  arrange(pvalue) %>%
  mutate(log2FoldChange = round(log2FoldChange, 3),
         pvalue = signif(pvalue, 3),
         padj = signif(padj, 3),
         FDR = signif(FDR, 3))

write.csv(all_res, file.path(DAT, "pseudobulk_CNNM2_all_levels.csv"),
          row.names = FALSE)
options(width = 200)
cat("\n========= Pseudobulk DESeq2 CNNM2 — all levels =========\n")
print(all_res, row.names = FALSE)

# Forest plot of effect sizes
all_res$cell_label <- paste0(all_res$cell_type, " (n=",
                              all_res$n_donors_case, "/", all_res$n_donors_ctrl, ")")
all_res$cell_label <- factor(all_res$cell_label,
                              levels = all_res$cell_label[order(all_res$log2FoldChange)])
all_res$sig <- ifelse(all_res$pvalue < 0.05, "P<0.05",
                ifelse(all_res$pvalue < 0.1,  "P<0.1", "NS"))

p_forest <- ggplot(all_res, aes(x = log2FoldChange, y = cell_label,
                                  color = sig)) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = log2FoldChange - 1.96 * lfcSE,
                     xmax = log2FoldChange + 1.96 * lfcSE), height = 0.3) +
  geom_vline(xintercept = 0, linetype = 2, color = "gray60") +
  facet_grid(level ~ ., scales = "free_y", space = "free_y") +
  scale_color_manual(values = c("P<0.05" = "#cc3a36", "P<0.1" = "#f4a261",
                                 "NS" = "gray60")) +
  labs(x = "log2 fold change (Case vs Control)",
       y = NULL,
       title = "CNNM2 pseudobulk DESeq2 — donor-level DE",
       subtitle = "Cell type (n_case_donors / n_control_donors) | error bars = 95% CI") +
  theme_bw(base_size = 11) +
  theme(strip.text = element_text(face = "bold"),
        legend.position = "top")

ggsave(file.path(FIG, "CNNM2_pseudobulk_forest.pdf"), p_forest,
       width = 9, height = 7)
ggsave(file.path(FIG, "CNNM2_pseudobulk_forest.png"), p_forest,
       width = 9, height = 7, dpi = 200, bg = "white")

cat("\n[", as.character(Sys.time()), "] DONE\n")
