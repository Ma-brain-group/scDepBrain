# =================================================================
# TRS correlation: scPagwas vs scDRS — Inhibitory + Excitatory subtypes
# Two levels per class:
#   (a) Subtype-level — each dot = one subtype (median TRS), label points
#   (b) Single-cell level — each dot = one cell, regression line + Pearson r
# Style: mimics user's published Fig (panel d subtype-level + panel b single-cell)
# Based on HW (Howard 2019) GWAS — primary discovery GWAS
# =================================================================
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(dplyr); library(ggrepel)
})

DAT_F <- "/mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_HW_v2_figures_20260617/data"
SCD_F <- "/mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_scDRS_20260618/scores_1k"
FIG   <- "/mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_HW_v2_figures_20260617/figures"

# ---------- Load per-cell scPagwas TRS + final_anno + anno ----------
cat("[", as.character(Sys.time()), "] Loading data ...\n")
meta <- fread(cmd = paste("zcat", file.path(DAT_F, "scPagwas_HW_metadata.csv.gz")),
              select = c("cell_id","anno","final_anno","scPagwas.TRS.Score"))
setnames(meta, c("cell_id","anno","final_anno","scPagwas_TRS"))

# ---------- Load scDRS norm_score (1000 ctrl) ----------
scd <- fread(cmd = paste("zcat", file.path(SCD_F, "HW.full_score.gz")), select = c(1, 3))
setnames(scd, c("cell_id","scDRS_TRS"))

# Merge per cell
d <- merge(meta, scd, by = "cell_id")
cat("Merged cells:", nrow(d), "\n")
cat("Inhibitory cells:", nrow(d[anno == "Inhibitory neurons"]), "\n")
cat("Excitatory cells:", nrow(d[anno == "Excitatory neurons"]), "\n")

# ---------- Subtype name cleanup for labels ----------
# Inhibitory subtypes (use existing names)
IN_SUBTYPES <- c("In_PVALB","In_LAMP5","In_SST","In_VIP","In_SHANK2","In_CALM1")
EX_SUBTYPES <- c("Ex-L2/4","Ex-L2/3","Ex-L4/6","Ex-L5","Ex-L5/6","Ex-L6","Ex-NRGN","Ex_mix")

# ---------- Plot 1: Subtype-level correlation ----------
mk_subtype_plot <- function(d_sub, subtypes, title_str) {
  d_sub <- d_sub[final_anno %in% subtypes]
  agg <- d_sub[, .(
    median_scPagwas = median(scPagwas_TRS, na.rm = TRUE),
    median_scDRS   = median(scDRS_TRS,   na.rm = TRUE),
    n_cell = .N
  ), by = final_anno]

  # Pearson correlation on the 6 (or 8) points
  r <- cor(agg$median_scPagwas, agg$median_scDRS,
           method = "pearson", use = "complete.obs")
  pval <- cor.test(agg$median_scPagwas, agg$median_scDRS,
                    method = "pearson")$p.value

  cat("\n", title_str, "— Subtype-level:\n")
  cat("  Pearson r =", round(r, 3), "  P =", signif(pval, 3), "\n")
  print(agg)

  ann_text <- sprintf("r = %.2f, P = %s",
                       r,
                       ifelse(pval < 1e-3,
                              format(pval, scientific = TRUE, digits = 2),
                              sprintf("%.4f", pval)))

  ggplot(agg, aes(x = median_scPagwas, y = median_scDRS)) +
    geom_smooth(method = "lm", se = TRUE, color = "#bfa07a",
                fill = "gray85", alpha = 0.5, linewidth = 0.7,
                linetype = "dashed") +
    geom_point(size = 4, color = "#2d7d7e", alpha = 0.85) +
    geom_text_repel(aes(label = final_anno), size = 3.7,
                    color = "black", fontface = "bold",
                    box.padding = 0.5, max.overlaps = 30,
                    segment.size = 0.3, segment.color = "gray50") +
    annotate("text", x = -Inf, y = Inf, label = ann_text,
             hjust = -0.2, vjust = 1.5, size = 4.5, fontface = "bold") +
    labs(x = "scPagwas TRS score (median)",
         y = "scDRS score (median)",
         title = paste0("TRS correlation — ", title_str, " (subtype-level)")) +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 13),
          panel.grid.major = element_line(color = "gray92", linewidth = 0.3))
}

# ---------- Plot 2: Single-cell-level correlation ----------
mk_singlecell_plot <- function(d_sub, subtypes, title_str) {
  d_sub <- d_sub[final_anno %in% subtypes]
  # Subsample to ~30k cells if too many (for fast plotting)
  set.seed(42)
  if (nrow(d_sub) > 50000) {
    d_plot <- d_sub[sample(.N, 50000)]
  } else {
    d_plot <- d_sub
  }
  cat("\n", title_str, "— Single-cell-level:\n")
  cat("  Total cells:", nrow(d_sub), " | Plotted:", nrow(d_plot), "\n")

  # Full Pearson on all data
  r <- cor(d_sub$scPagwas_TRS, d_sub$scDRS_TRS,
           method = "pearson", use = "complete.obs")
  pval <- cor.test(d_sub$scPagwas_TRS, d_sub$scDRS_TRS,
                    method = "pearson")$p.value

  cat("  Pearson r =", round(r, 3), "  P =", signif(pval, 3), "\n")

  ann_text <- sprintf("r = %.2f, P %s",
                       r,
                       ifelse(pval < 2.2e-16, "< 2.2e-16",
                              paste("=", format(pval, scientific = TRUE, digits = 2))))

  ggplot(d_plot, aes(x = scPagwas_TRS, y = scDRS_TRS)) +
    geom_point(size = 0.4, color = "#2d7d7e", alpha = 0.45) +
    geom_smooth(method = "lm", se = FALSE, color = "#e7a05c",
                linewidth = 0.9) +
    annotate("text", x = -Inf, y = Inf, label = ann_text,
             hjust = -0.15, vjust = 1.8, size = 4.5, fontface = "bold") +
    labs(x = "scPagwas TRS score (per cell)",
         y = "scDRS score (per cell)",
         title = paste0("TRS correlation — ", title_str, " (single-cell level)")) +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 13),
          panel.grid.major = element_line(color = "gray92", linewidth = 0.3))
}

# ---------- Generate 4 figures ----------
d_in <- d[anno == "Inhibitory neurons"]
d_ex <- d[anno == "Excitatory neurons"]

p1 <- mk_subtype_plot(d_in, IN_SUBTYPES, "Inhibitory neuron subtypes")
ggsave(file.path(FIG, "TRS_correlation_HW_inhibitory_subtype.pdf"), p1,
       width = 6.5, height = 5.5)
ggsave(file.path(FIG, "TRS_correlation_HW_inhibitory_subtype.png"), p1,
       width = 6.5, height = 5.5, dpi = 300, bg = "white")

p2 <- mk_singlecell_plot(d_in, IN_SUBTYPES, "Inhibitory neurons")
ggsave(file.path(FIG, "TRS_correlation_HW_inhibitory_singlecell.pdf"), p2,
       width = 6.5, height = 5.5)
ggsave(file.path(FIG, "TRS_correlation_HW_inhibitory_singlecell.png"), p2,
       width = 6.5, height = 5.5, dpi = 300, bg = "white")

p3 <- mk_subtype_plot(d_ex, EX_SUBTYPES, "Excitatory neuron subtypes")
ggsave(file.path(FIG, "TRS_correlation_HW_excitatory_subtype.pdf"), p3,
       width = 6.5, height = 5.5)
ggsave(file.path(FIG, "TRS_correlation_HW_excitatory_subtype.png"), p3,
       width = 6.5, height = 5.5, dpi = 300, bg = "white")

p4 <- mk_singlecell_plot(d_ex, EX_SUBTYPES, "Excitatory neurons")
ggsave(file.path(FIG, "TRS_correlation_HW_excitatory_singlecell.pdf"), p4,
       width = 6.5, height = 5.5)
ggsave(file.path(FIG, "TRS_correlation_HW_excitatory_singlecell.png"), p4,
       width = 6.5, height = 5.5, dpi = 300, bg = "white")

cat("\nWrote 4 figures to:", FIG, "\n")
cat("\n[", as.character(Sys.time()), "] DONE\n")
