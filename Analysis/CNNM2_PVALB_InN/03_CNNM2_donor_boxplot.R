# =================================================================
# Donor-level boxplot of CNNM2 pseudobulk expression
# 3 levels: All cells / Inhibitory neurons / In_PVALB subtype
# Each point = 1 donor (37 Case + 35 Control)
# DESeq2 P value annotated on each panel
# =================================================================
suppressPackageStartupMessages({
  library(Seurat); library(DESeq2); library(Matrix); library(dplyr)
  library(ggplot2); library(ggpubr); library(patchwork)
})
utils::assignInNamespace(".Deprecate", function(...) invisible(NULL),
                         ns = "SeuratObject")

DATA_DIR <- "/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/replication_4_MDD_brain_cells"
OUT <- "/mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_CNNM2_replication_20260619"
FIG <- file.path(OUT, "figures"); DAT <- file.path(OUT, "data")

GENE <- "CNNM2"

# =================================================================
# Helper: build pseudobulk + DESeq2 + return per-donor CNNM2 values
# =================================================================
pseudobulk_cnnm2 <- function(obj, cell_filter = NULL, label = "X") {
  cat("\n[", as.character(Sys.time()), "] Processing:", label, "\n")
  md <- obj@meta.data
  md$cell_id <- rownames(md)

  if (!is.null(cell_filter)) {
    keep_idx <- cell_filter(md)
    md <- md[keep_idx, ]
    cat("Cells after filter:", nrow(md), "\n")
  }

  # Filter to Case/Control
  keep <- md$phenotype %in% c("Case","Control")
  md <- md[keep, ]
  cat("Cells with Case/Control phenotype:", nrow(md), "\n")

  # Get raw counts for these cells
  counts <- GetAssayData(obj, assay = "RNA", slot = "counts")[, md$cell_id]

  # Aggregate per donor (sample)
  donors <- factor(md$sample)
  ind <- Matrix::sparse.model.matrix(~ 0 + donors)
  colnames(ind) <- gsub("^donors", "", colnames(ind))
  pseudo <- as.matrix(counts %*% ind)
  cat("Pseudobulk:", paste(dim(pseudo), collapse = " x "), "\n")

  # Donor-level metadata
  smd <- md %>% select(sample, phenotype, any_of(c("sex","batch","age"))) %>%
    distinct() %>% as.data.frame()
  rownames(smd) <- smd$sample
  smd <- smd[colnames(pseudo), , drop = FALSE]
  smd$phenotype <- factor(smd$phenotype, levels = c("Control","Case"))

  # Min donors per phenotype
  cat("Donors per phenotype:\n"); print(table(smd$phenotype))

  # Filter low-count genes for DESeq2 stability
  pseudo_f <- pseudo[rowSums(pseudo) >= 10, ]

  # Build DESeq2 with phenotype + sex
  design <- ~ phenotype
  if ("sex" %in% colnames(smd) && length(unique(smd$sex)) > 1 &&
      all(table(smd$sex, smd$phenotype) > 0)) {
    design <- ~ sex + phenotype
  }
  dds <- DESeqDataSetFromMatrix(countData = pseudo_f, colData = smd, design = design)
  dds <- DESeq(dds, quiet = TRUE)

  # CNNM2 result
  res <- results(dds, contrast = c("phenotype","Case","Control"))
  cnnm2_res <- if (GENE %in% rownames(res)) as.data.frame(res[GENE, ]) else NULL
  if (is.null(cnnm2_res)) stop("CNNM2 missing from DE result")
  cat("CNNM2 DE: log2FC =", round(cnnm2_res$log2FoldChange, 3),
      ", P =", signif(cnnm2_res$pvalue, 3), "\n")

  # Get normalized CNNM2 per donor (size-factor normalized then log2(x+1))
  norm_counts <- counts(dds, normalized = TRUE)
  cnnm2_per_donor <- log2(norm_counts[GENE, ] + 1)
  plot_df <- data.frame(
    donor = colnames(pseudo),
    phenotype = smd$phenotype,
    CNNM2_log2norm = cnnm2_per_donor,
    level = label,
    stringsAsFactors = FALSE
  )

  list(plot_df = plot_df, res = cnnm2_res, n_donors = nrow(smd))
}

# =================================================================
# Step 1: All cells + Inhibitory broad (from 2023_NC.rds)
# =================================================================
cat("\n[", as.character(Sys.time()), "] Loading 2023_NC.rds ...\n")
obj <- readRDS(file.path(DATA_DIR, "2023_NC.rds"))
DefaultAssay(obj) <- "RNA"
if (inherits(obj[["RNA"]], "Assay5")) {
  obj <- JoinLayers(obj, assay = "RNA")
  obj[["RNA"]] <- as(obj[["RNA"]], "Assay")
}
cat("Total cells:", ncol(obj), "\n")

# All cells
all_res <- pseudobulk_cnnm2(obj, cell_filter = NULL, label = "All cells")

# Inhibitory broad
inn_res <- pseudobulk_cnnm2(obj,
  cell_filter = function(md) md$broad_cell_type == "InN",
  label = "Inhibitory neurons (InN)")

rm(obj); gc()

# =================================================================
# Step 2: In_PVALB from NC_inhibitory rds (recluster-based PVALB)
# =================================================================
cat("\n[", as.character(Sys.time()), "] Loading NC_inhibitory ...\n")
obj_in <- readRDS(file.path(DATA_DIR, "NC_inhibitory_neurons_reannotation.rds"))
DefaultAssay(obj_in) <- "RNA"
if (inherits(obj_in[["RNA"]], "Assay5")) {
  obj_in <- JoinLayers(obj_in, assay = "RNA")
  obj_in[["RNA"]] <- as(obj_in[["RNA"]], "Assay")
}
cat("Inhibitory cells:", ncol(obj_in), "\n")

pvalb_res <- pseudobulk_cnnm2(obj_in,
  cell_filter = function(md) md$RNA_snn_res.0.1_cell_type == "In_PVALB",
  label = "In_PVALB")

rm(obj_in); gc()

# =================================================================
# Combine and plot
# =================================================================
plot_df <- bind_rows(all_res$plot_df, inn_res$plot_df, pvalb_res$plot_df)
plot_df$level <- factor(plot_df$level,
  levels = c("All cells", "Inhibitory neurons (InN)", "In_PVALB"))
plot_df$phenotype <- factor(plot_df$phenotype, levels = c("Control","Case"))

# Annotation labels (DESeq2 result text per panel)
ann_df <- data.frame(
  level = c("All cells", "Inhibitory neurons (InN)", "In_PVALB"),
  log2FC = c(NA, inn_res$res$log2FoldChange, pvalb_res$res$log2FoldChange),
  P = c(NA, inn_res$res$pvalue, pvalb_res$res$pvalue),
  n_case = c(all_res$n_donors, inn_res$n_donors, pvalb_res$n_donors),
  stringsAsFactors = FALSE
)
ann_df$log2FC[1] <- all_res$res$log2FoldChange
ann_df$P[1]      <- all_res$res$pvalue
ann_df$label <- sprintf("DESeq2 P = %s\nlog2FC = %+.2f",
                         format(ann_df$P, digits = 2, scientific = TRUE),
                         ann_df$log2FC)
ann_df$level <- factor(ann_df$level, levels = levels(plot_df$level))

write.csv(ann_df, file.path(DAT, "donor_boxplot_CNNM2_stats.csv"),
          row.names = FALSE)
write.csv(plot_df, file.path(DAT, "donor_boxplot_CNNM2_data.csv"),
          row.names = FALSE)

# Plot: facet by level
p <- ggplot(plot_df, aes(x = phenotype, y = CNNM2_log2norm, fill = phenotype)) +
  geom_boxplot(width = 0.55, alpha = 0.75, outlier.shape = NA,
               color = "black", lwd = 0.4) +
  geom_jitter(aes(color = phenotype), shape = 21, size = 2.6,
              stroke = 0.4, width = 0.18, alpha = 0.9) +
  geom_text(data = ann_df,
            aes(x = 1.5, y = Inf, label = label),
            hjust = 0.5, vjust = 1.5, size = 3.5,
            inherit.aes = FALSE) +
  facet_wrap(~ level, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = c("Control" = "#3777ac", "Case" = "#cc3a36")) +
  scale_color_manual(values = c("Control" = "#1d4f7b", "Case" = "#8a1a18")) +
  labs(x = NULL, y = "Pseudobulk CNNM2 expression  log2(norm + 1)",
       title = "CNNM2 donor-level pseudobulk expression — MDD vs Control",
       subtitle = "Each point = 1 donor | DESeq2 + sex covariate | GSE213982 (n=37 MDD vs 35 Control)") +
  theme_bw(base_size = 12) +
  theme(strip.text = element_text(face = "bold", size = 11),
        strip.background = element_rect(fill = "gray95"),
        axis.text.x = element_text(size = 11, face = "bold"),
        legend.position = "none",
        panel.grid.minor = element_blank())

ggsave(file.path(FIG, "CNNM2_donor_boxplot_3levels.pdf"), p,
       width = 11, height = 5.5)
ggsave(file.path(FIG, "CNNM2_donor_boxplot_3levels.png"), p,
       width = 11, height = 5.5, dpi = 200, bg = "white")

options(width = 200)
cat("\n========== Annotation summary ==========\n")
print(ann_df)

cat("\nWrote:\n",
    file.path(FIG, "CNNM2_donor_boxplot_3levels.pdf"), "\n",
    file.path(FIG, "CNNM2_donor_boxplot_3levels.png"), "\n", sep = "")
cat("\n[", as.character(Sys.time()), "] DONE\n")
