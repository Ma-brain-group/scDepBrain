# =================================================================
# HW (Howard 2019) broad cell type dot charts — using OUR v2 results
# Adapted from the original ggdotchart_for_Pvalue_for_celltypes.R style
# =================================================================
suppressPackageStartupMessages({
  library(ggplot2); library(ggpubr)
})

# ---- output dir ----
OUT <- "/mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_HW_v2_figures_20260617"
FIG <- file.path(OUT, "figures")
DAT <- file.path(OUT, "data")

# ---- canonical 8 broad cell types (paper order, top→bottom on plot) ----
CT_ORDER <- c("Excitatory.neurons", "Inhibitory.neurons", "Purkinje.neurons",
              "Endothelial.cells", "Oligodendrocytes", "OPCs",
              "Microglia", "Astrocytes")
# rev() = bottom→top factor order; palette matches that order
PALETTE <- c("#d5231d","#e88f18","#b698c5","#e47faf","#3777ac","#a05528","#8e4c99","#4ea64a")

# ---- universal plotting fn (paper style) ----
plot_dotchart <- function(df, file_out, title = NULL) {
  df$cell_type <- factor(df$cell_type, levels = rev(CT_ORDER))
  df$`-log10(P)` <- -log10(df$P)
  df$Significance <- factor(ifelse(df$P < 0.05, 1, 0))

  pdf(file_out, width = 4.5, height = 3)
  p <- ggdotchart(df, x = "cell_type", y = "-log10(P)",
        color = "cell_type",
        palette = PALETTE,
        sorting = "descending",
        add = "segments",
        add.params = list(color = "lightgray", size = 2),
        group = "cell_type",
        dot.size = 7,
        label = round(df$`-log10(P)`, 3),
        font.label = list(color = "black", size = 8, vjust = 0.5),
        ggtheme = theme_pubr()) +
    coord_flip() +
    geom_hline(yintercept = -log10(0.05), linetype = 2, color = "lightgray") +
    theme(axis.text.x = element_text(size = 10, color = "black"),
          axis.text.y = element_text(size = 10, color = "black"),
          legend.position = "none") +
    scale_y_continuous(limits = c(0, 10), breaks = seq(0, 10, 2))
  if (!is.null(title)) p <- p + ggtitle(title)
  print(p)
  dev.off()
  cat("Wrote:", file_out, "\n")
}

# =================================================================
# 1. MAGMA  — our v2 set-annot HW result (matches paper's --set-annot top 10%)
# =================================================================
cat("\n=== 1. MAGMA v2 set-annot HW ===\n")
magma <- read.table(
  "/mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_MAGMA_PGC2025_20260617/enrichment/HW_anno_set_v2.gsa.out",
  header = TRUE, stringsAsFactors = FALSE)
# Columns: VARIABLE TYPE NGENES BETA BETA_STD SE P
magma_df <- data.frame(cell_type = magma$VARIABLE, P = magma$P)
# Rename to paper style (replace _ with .)
magma_df$cell_type <- gsub("_", ".", magma_df$cell_type)
print(magma_df)
write.csv(magma_df, file.path(DAT, "magma_v2_setannot_HW.csv"), row.names = FALSE)
plot_dotchart(magma_df, file.path(FIG, "ggdotchart_of_magma_v2.pdf"),
              title = "MAGMA (set-annot v2)")

# =================================================================
# 1b. MAGMA gene-covar (Bryois MAGMA_CellTyping standard) — alternative
# =================================================================
cat("\n=== 1b. MAGMA v2 gene-covar HW ===\n")
magma_cov <- read.table(
  "/mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_MAGMA_PGC2025_20260617/enrichment/HW_anno_covar_v2.gsa.out",
  header = TRUE, stringsAsFactors = FALSE)
magma_cov_df <- data.frame(cell_type = magma_cov$VARIABLE, P = magma_cov$P)
magma_cov_df$cell_type <- gsub("_", ".", magma_cov_df$cell_type)
print(magma_cov_df)
write.csv(magma_cov_df, file.path(DAT, "magma_v2_genecovar_HW.csv"), row.names = FALSE)
plot_dotchart(magma_cov_df, file.path(FIG, "ggdotchart_of_magma_v2_genecovar.pdf"),
              title = "MAGMA_CellTyping (gene-covar, Bryois)")

# =================================================================
# 2. scPagwas — HW paper-params (Merge_celltype_p)
# =================================================================
cat("\n=== 2. scPagwas HW paper-params ===\n")
scp <- read.csv(
  "/mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_scPagwas_HW_paperparams_20260615/cell_type_pvalues_paperparams.csv",
  stringsAsFactors = FALSE)
# Columns: celltype, pvalue
scp_df <- data.frame(cell_type = scp$celltype, P = scp$pvalue)
# Rename spaces → dots to match paper style
for (orig in c("Excitatory neurons","Inhibitory neurons","Purkinje neurons","Endothelial cells")) {
  scp_df$cell_type[scp_df$cell_type == orig] <- gsub(" ", ".", orig)
}
print(scp_df)
write.csv(scp_df, file.path(DAT, "scPagwas_HW.csv"), row.names = FALSE)
plot_dotchart(scp_df, file.path(FIG, "ggdotchart_of_scPagwas_v2.pdf"),
              title = "scPagwas (Merge_celltype_p)")

# =================================================================
# 3. LDSC-SEG  — will fill in once pipeline finishes
# =================================================================
LDSC_DIR <- "/mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_LDSC_PGC2025_20260617/enrichment"
ldsc_files <- list.files(LDSC_DIR, pattern = "^HW__major_.*\\.results$", full.names = TRUE)
if (length(ldsc_files) > 0) {
  cat("\n=== 3. LDSC-SEG HW ===\n")
  ldsc_rows <- lapply(ldsc_files, function(f) {
    ct <- sub("^HW__major_", "", sub("\\.results$", "", basename(f)))
    tab <- read.table(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    # The cell-type-specific row is the LAST annotation appended (after baseline)
    last_row <- tab[nrow(tab), ]
    data.frame(cell_type = ct,
               P = last_row$Enrichment_p,
               Coefficient_z = last_row$Coefficient_z_score)
  })
  ldsc_df <- do.call(rbind, ldsc_rows)
  ldsc_df$cell_type <- gsub("_", ".", ldsc_df$cell_type)
  print(ldsc_df)
  write.csv(ldsc_df, file.path(DAT, "ldsc_HW.csv"), row.names = FALSE)
  plot_dotchart(ldsc_df[, c("cell_type","P")],
                file.path(FIG, "ggdotchart_of_ldsc_v2.pdf"),
                title = "LDSC-SEG")
} else {
  cat("\n=== 3. LDSC-SEG: not ready yet. Rerun this script when done. ===\n")
}

# =================================================================
# 4. scDRS  — not yet run. Skip if file absent.
# =================================================================
scDRS_path <- "/mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_scDRS_HW_20260617/MDD.scdrs_ct.anno"
if (file.exists(scDRS_path)) {
  cat("\n=== 4. scDRS HW ===\n")
  scDRS <- read.table(scDRS_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  scDRS_df <- data.frame(cell_type = scDRS$X, P = scDRS$assoc_mcp)
  scDRS_df$cell_type <- gsub("_", ".", scDRS_df$cell_type)
  print(scDRS_df)
  write.csv(scDRS_df, file.path(DAT, "scDRS_HW.csv"), row.names = FALSE)
  plot_dotchart(scDRS_df, file.path(FIG, "ggdotchart_of_scDRS_v2.pdf"),
                title = "scDRS")
} else {
  cat("\n=== 4. scDRS: not run yet. Skip. ===\n")
}

cat("\n=== Outputs ===\n")
cat("Data:   ", DAT, "\n")
cat("Figures:", FIG, "\n")
system(paste("ls -lh", FIG))
