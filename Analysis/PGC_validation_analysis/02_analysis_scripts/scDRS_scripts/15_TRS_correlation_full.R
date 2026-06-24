# =================================================================
# TRS correlation — FULL pipeline
# scPagwas TRS (two sources) × scDRS norm_score, three GWAS, two classes, two levels
#
# Sources:
#   HW_old : scPagwas.TRS.Score_Howard2019_OLD vs HW scDRS  (legacy/published)
#   HW_new : scPagwas.TRS.Score (latest v2 run) vs HW scDRS
#   EUR    : scPagwas.TRS.Score (EUR run) vs EUR scDRS
#   DIV    : scPagwas.TRS.Score (DIV run) vs DIV scDRS
#
# For each source: 4 figures
#   - Inhibitory subtype-level correlation
#   - Inhibitory single-cell-level correlation
#   - Excitatory subtype-level correlation
#   - Excitatory single-cell-level correlation
# Total: 16 figures
# Style: matches user's published Fig (panel b + d)
# =================================================================
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(dplyr); library(ggrepel)
})

DAT_F <- "/mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_HW_v2_figures_20260617/data"
SCD_F <- "/mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_scDRS_20260618/scores_1k"
FIG   <- "/mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_HW_v2_figures_20260617/figures"

# Subtype sets
IN_SUBTYPES <- c("In_PVALB","In_LAMP5","In_SST","In_VIP","In_SHANK2","In_CALM1")
EX_SUBTYPES <- c("Ex-L2/4","Ex-L2/3","Ex-L4/6","Ex-L5","Ex-L5/6","Ex-L6","Ex-NRGN","Ex_mix")

# ---------- Loaders ----------
load_meta <- function(meta_path, scp_col) {
  d <- fread(cmd = paste("zcat", meta_path),
             select = c("cell_id","anno","final_anno", scp_col))
  setnames(d, c("cell_id","anno","final_anno","scPagwas_TRS"))
  d
}
load_scd <- function(path) {
  d <- fread(cmd = paste("zcat", path), select = c(1, 3))
  setnames(d, c("cell_id","scDRS_TRS")); d
}

# ---------- Subtype-level scatter ----------
mk_subtype_plot <- function(d_sub, subtypes, title_str) {
  d_sub <- d_sub[final_anno %in% subtypes]
  agg <- d_sub[, .(
    median_scPagwas = median(scPagwas_TRS, na.rm = TRUE),
    median_scDRS   = median(scDRS_TRS,   na.rm = TRUE),
    n_cell = .N), by = final_anno]

  r <- cor(agg$median_scPagwas, agg$median_scDRS, method = "pearson", use = "complete.obs")
  pval <- cor.test(agg$median_scPagwas, agg$median_scDRS, method = "pearson")$p.value

  ann_text <- sprintf("r = %.2f, P = %s", r,
                      ifelse(pval < 1e-3,
                             format(pval, scientific = TRUE, digits = 2),
                             sprintf("%.4f", pval)))
  cat("  ", title_str, ": subtype r =", round(r,3), " P =", signif(pval,3), "\n")

  ggplot(agg, aes(x = median_scPagwas, y = median_scDRS)) +
    geom_smooth(method = "lm", se = TRUE, color = "#bfa07a",
                fill = "gray85", alpha = 0.5, linewidth = 0.7, linetype = "dashed") +
    geom_point(size = 4, color = "#2d7d7e", alpha = 0.85) +
    geom_text_repel(aes(label = final_anno), size = 3.7, color = "black",
                    fontface = "bold", box.padding = 0.5, max.overlaps = 30,
                    segment.size = 0.3, segment.color = "gray50") +
    annotate("text", x = -Inf, y = Inf, label = ann_text,
             hjust = -0.18, vjust = 1.5, size = 4.5, fontface = "bold") +
    labs(x = "scPagwas TRS score (median)",
         y = "scDRS score (median)",
         title = paste0("TRS correlation — ", title_str, " (subtype-level)")) +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 13),
          panel.grid.major = element_line(color = "gray92", linewidth = 0.3))
}

# ---------- Single-cell-level scatter ----------
mk_singlecell_plot <- function(d_sub, subtypes, title_str) {
  d_sub <- d_sub[final_anno %in% subtypes]
  set.seed(42)
  d_plot <- if (nrow(d_sub) > 50000) d_sub[sample(.N, 50000)] else d_sub

  r <- cor(d_sub$scPagwas_TRS, d_sub$scDRS_TRS,
           method = "pearson", use = "complete.obs")
  pval <- cor.test(d_sub$scPagwas_TRS, d_sub$scDRS_TRS, method = "pearson")$p.value
  ann_text <- sprintf("r = %.2f, P %s", r,
                      ifelse(pval < 2.2e-16, "< 2.2e-16",
                             paste("=", format(pval, scientific = TRUE, digits = 2))))
  cat("  ", title_str, ": single-cell r =", round(r,3), " P =", signif(pval,3),
      " (n =", nrow(d_sub), ")\n")

  ggplot(d_plot, aes(x = scPagwas_TRS, y = scDRS_TRS)) +
    geom_point(size = 0.4, color = "#2d7d7e", alpha = 0.45) +
    geom_smooth(method = "lm", se = FALSE, color = "#e7a05c", linewidth = 0.9) +
    annotate("text", x = -Inf, y = Inf, label = ann_text,
             hjust = -0.15, vjust = 1.8, size = 4.5, fontface = "bold") +
    labs(x = "scPagwas TRS score (per cell)",
         y = "scDRS score (per cell)",
         title = paste0("TRS correlation — ", title_str, " (single-cell)")) +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 13),
          panel.grid.major = element_line(color = "gray92", linewidth = 0.3))
}

# ---------- Master function: one source combination -> 4 figures ----------
run_source <- function(meta_path, scd_path, scp_col, source_tag, source_label) {
  cat("\n==========", source_tag, "==========\n")
  m <- load_meta(meta_path, scp_col)
  s <- load_scd(scd_path)
  d <- merge(m, s, by = "cell_id")
  cat("Merged:", nrow(d), "cells\n")
  d_in <- d[anno == "Inhibitory neurons"]
  d_ex <- d[anno == "Excitatory neurons"]
  cat("  Inhibitory:", nrow(d_in), " Excitatory:", nrow(d_ex), "\n")

  p1 <- mk_subtype_plot(d_in, IN_SUBTYPES,
                         paste0("Inhibitory neuron subtypes (", source_label, ")"))
  ggsave(file.path(FIG, paste0("TRS_corr_", source_tag, "_inhibitory_subtype.pdf")), p1,
         width = 6.5, height = 5.5)
  ggsave(file.path(FIG, paste0("TRS_corr_", source_tag, "_inhibitory_subtype.png")), p1,
         width = 6.5, height = 5.5, dpi = 300, bg = "white")

  p2 <- mk_singlecell_plot(d_in, IN_SUBTYPES,
                            paste0("Inhibitory neurons (", source_label, ")"))
  ggsave(file.path(FIG, paste0("TRS_corr_", source_tag, "_inhibitory_singlecell.pdf")), p2,
         width = 6.5, height = 5.5)
  ggsave(file.path(FIG, paste0("TRS_corr_", source_tag, "_inhibitory_singlecell.png")), p2,
         width = 6.5, height = 5.5, dpi = 300, bg = "white")

  p3 <- mk_subtype_plot(d_ex, EX_SUBTYPES,
                         paste0("Excitatory neuron subtypes (", source_label, ")"))
  ggsave(file.path(FIG, paste0("TRS_corr_", source_tag, "_excitatory_subtype.pdf")), p3,
         width = 6.5, height = 5.5)
  ggsave(file.path(FIG, paste0("TRS_corr_", source_tag, "_excitatory_subtype.png")), p3,
         width = 6.5, height = 5.5, dpi = 300, bg = "white")

  p4 <- mk_singlecell_plot(d_ex, EX_SUBTYPES,
                            paste0("Excitatory neurons (", source_label, ")"))
  ggsave(file.path(FIG, paste0("TRS_corr_", source_tag, "_excitatory_singlecell.pdf")), p4,
         width = 6.5, height = 5.5)
  ggsave(file.path(FIG, paste0("TRS_corr_", source_tag, "_excitatory_singlecell.png")), p4,
         width = 6.5, height = 5.5, dpi = 300, bg = "white")
}

# ---------- Run all 4 source combinations ----------
SOURCES <- list(
  list(meta = file.path(DAT_F, "scPagwas_HW_metadata.csv.gz"),
       scd  = file.path(SCD_F, "HW.full_score.gz"),
       scp_col = "scPagwas.TRS.Score_Howard2019_OLD",
       tag = "HW_oldTRS",
       label = "HW, legacy scPagwas TRS"),
  list(meta = file.path(DAT_F, "scPagwas_HW_metadata.csv.gz"),
       scd  = file.path(SCD_F, "HW.full_score.gz"),
       scp_col = "scPagwas.TRS.Score",
       tag = "HW_newTRS",
       label = "HW, current scPagwas v2 TRS"),
  list(meta = file.path(DAT_F, "scPagwas_EUR_metadata.csv.gz"),
       scd  = file.path(SCD_F, "EUR.full_score.gz"),
       scp_col = "scPagwas.TRS.Score",
       tag = "EUR",
       label = "PGC MDD 2025 EUR"),
  list(meta = file.path(DAT_F, "scPagwas_DIV_metadata.csv.gz"),
       scd  = file.path(SCD_F, "DIV.full_score.gz"),
       scp_col = "scPagwas.TRS.Score",
       tag = "DIV",
       label = "PGC MDD 2025 trans-ancestry")
)

for (s in SOURCES) {
  run_source(s$meta, s$scd, s$scp_col, s$tag, s$label)
}

cat("\nWrote 16 TRS correlation figures to:", FIG, "\n")
system(paste("ls -1", FIG, "| grep TRS_corr_ | sort"))
