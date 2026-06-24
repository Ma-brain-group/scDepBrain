# =================================================================
# CNNM2 expression × scDRS TRS — per-cell correlation + quintile-bin boxplot
# 3 GWAS (HW, EUR, DIV) × 3 cell subsets (All / Inhibitory / In_PVALB)
# Two plot types per combo: scatter + quintile-bin boxplot
# Layout: 3 figures (one per GWAS), each 2 rows × 3 cols (6 panels)
# =================================================================
suppressPackageStartupMessages({
  library(Seurat); library(data.table); library(ggplot2); library(dplyr)
  library(patchwork); library(scales)
})
utils::assignInNamespace(".Deprecate", function(...) invisible(NULL),
                         ns = "SeuratObject")

OUT_F <- "/mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_HW_v2_figures_20260617"
FIG <- file.path(OUT_F, "figures"); DAT <- file.path(OUT_F, "data")
SCD_F <- "/mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_scDRS_20260618/scores_1k"

GENE <- "CNNM2"
BIN_COLORS <- c("1" = "#dfb6d3", "2" = "#5a9b56", "3" = "#e8a049",
                "4" = "#e8736f", "5" = "#a8d886")

# ---------- Step 1: Load CNNM2 expression + metadata ----------
cat("[", as.character(Sys.time()), "] Loading scPagwas rds for CNNM2 + metadata ...\n")
obj <- readRDS("/mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_scPagwas_HW_paperparams_20260615/MDD_scPagwas_HW_paperparams.rds")
DefaultAssay(obj) <- "RNA"
if (inherits(obj[["RNA"]], "Assay5")) {
  obj <- JoinLayers(obj, assay = "RNA")
  obj[["RNA"]] <- as(obj[["RNA"]], "Assay")
}
expr <- FetchData(obj, vars = GENE, slot = "data")
md <- obj@meta.data
base <- data.table(
  cell_id   = rownames(md),
  CNNM2     = expr[[GENE]],
  anno      = md$anno,
  final_anno = md$final_anno
)
rm(obj); gc()
cat("Total cells:", nrow(base), "\n")

# ---------- Step 2: Load scDRS scores for 3 GWAS ----------
read_scd <- function(path) {
  d <- fread(cmd = paste("zcat", path), select = c(1, 3))
  setnames(d, c("cell_id","scDRS")); d
}
scd_HW  <- read_scd(file.path(SCD_F, "HW.full_score.gz"))[, GWAS := "HW"]
scd_EUR <- read_scd(file.path(SCD_F, "EUR.full_score.gz"))[, GWAS := "EUR"]
scd_DIV <- read_scd(file.path(SCD_F, "DIV.full_score.gz"))[, GWAS := "DIV"]

# ---------- Step 3: Helper plot functions ----------
mk_scatter <- function(d_sub, gwas_label, subset_label) {
  d_sub <- d_sub[!is.na(CNNM2) & !is.na(scDRS)]
  # Subsample if too many cells
  set.seed(42)
  d_plot <- if (nrow(d_sub) > 50000) d_sub[sample(.N, 50000)] else d_sub
  r <- cor(d_sub$CNNM2, d_sub$scDRS, method = "pearson")
  pval <- cor.test(d_sub$CNNM2, d_sub$scDRS, method = "pearson")$p.value
  ann <- sprintf("r = %.2f, P %s", r,
                 ifelse(pval < 2.2e-16, "< 2.2e-16",
                        paste("=", format(pval, scientific = TRUE, digits = 2))))
  ggplot(d_plot, aes(x = CNNM2, y = scDRS)) +
    geom_point(size = 0.35, color = "#2d7d7e", alpha = 0.45) +
    geom_smooth(method = "lm", se = FALSE, color = "#e7a05c", linewidth = 0.9) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = Inf, ymax = Inf,
             fill = NA) +
    annotate("text", x = -Inf, y = Inf, label = ann,
             hjust = -0.1, vjust = 1.8, size = 3.6, fontface = "bold") +
    labs(x = "CNNM2 expression (log-normalized)", y = "scDRS TRS Score",
         title = paste0("CNNM2 × scDRS TRS in ", subset_label, " (", gwas_label, ")")) +
    theme_bw(base_size = 10) +
    theme(plot.title = element_text(face = "bold", size = 10),
          panel.grid.minor = element_blank())
}

mk_binbox <- function(d_sub, gwas_label, subset_label) {
  d_sub <- d_sub[!is.na(CNNM2) & !is.na(scDRS)]
  # Quintile bins by rank (handles many zeros via random tie-break)
  set.seed(42)
  d_sub[, bin := as.integer(cut(rank(CNNM2, ties.method = "random"),
                                  breaks = 5, labels = FALSE))]
  d_sub[, bin := factor(bin, levels = 1:5)]
  # Kruskal–Wallis test across bins
  kw <- kruskal.test(scDRS ~ bin, data = d_sub)
  ann <- sprintf("Kruskal–Wallis P %s",
                 ifelse(kw$p.value < 2.2e-16, "< 2.2e-16",
                        paste("=", format(kw$p.value, scientific = TRUE, digits = 2))))
  ymax <- quantile(d_sub$scDRS, 0.99, na.rm = TRUE)
  ymin <- quantile(d_sub$scDRS, 0.01, na.rm = TRUE)
  ggplot(d_sub, aes(x = bin, y = scDRS, fill = bin)) +
    geom_boxplot(outlier.size = 0.05, outlier.alpha = 0.2, lwd = 0.3) +
    scale_fill_manual(values = BIN_COLORS) +
    coord_cartesian(ylim = c(ymin, ymax * 1.05)) +
    annotate("text", x = 0.5, y = Inf, label = ann,
             hjust = 0, vjust = 1.8, size = 3.6, fontface = "bold") +
    labs(x = "CNNM2 expression gradient quintile bins", y = "scDRS TRS Score",
         title = paste0("CNNM2 quintile × scDRS TRS — ", subset_label, " (", gwas_label, ")")) +
    theme_bw(base_size = 10) +
    theme(legend.position = "none",
          plot.title = element_text(face = "bold", size = 10),
          panel.grid.minor = element_blank())
}

# ---------- Step 4: Run for each GWAS, make 6-panel figure ----------
SUBSETS <- list(
  "All cells"               = function(d) d,
  "Inhibitory neurons"      = function(d) d[anno == "Inhibitory neurons"],
  "In_PVALB"                = function(d) d[final_anno == "In_PVALB"]
)

stat_rows <- list()
for (g in c("HW","EUR","DIV")) {
  cat("\n==========", g, "==========\n")
  scd <- get(paste0("scd_", g))
  d <- merge(base, scd, by = "cell_id")
  cat("Merged:", nrow(d), "cells\n")

  panels <- list()
  i <- 1
  for (lvl in names(SUBSETS)) {
    d_lvl <- SUBSETS[[lvl]](d)
    panels[[i]] <- mk_scatter(d_lvl, g, lvl);          i <- i + 1
  }
  for (lvl in names(SUBSETS)) {
    d_lvl <- SUBSETS[[lvl]](d)
    panels[[i]] <- mk_binbox(d_lvl, g, lvl);           i <- i + 1
  }
  combined <- (panels[[1]] | panels[[2]] | panels[[3]]) /
              (panels[[4]] | panels[[5]] | panels[[6]]) +
    plot_annotation(
      title = paste0("CNNM2 expression × scDRS TRS — ", g),
      subtitle = "Top: per-cell scatter | Bottom: CNNM2-quintile bin boxplot | 3 cell-type levels (All, Inhibitory, In_PVALB)",
      tag_levels = "a",
      theme = theme(plot.title = element_text(face = "bold", size = 13),
                    plot.subtitle = element_text(size = 9, color = "gray35"))
    )
  ggsave(file.path(FIG, paste0("CNNM2_scDRS_correlation_", g, ".pdf")), combined,
         width = 14, height = 8)
  ggsave(file.path(FIG, paste0("CNNM2_scDRS_correlation_", g, ".png")), combined,
         width = 14, height = 8, dpi = 300, bg = "white")
  cat("Wrote: CNNM2_scDRS_correlation_", g, ".pdf\n", sep = "")

  # Stats summary
  for (lvl in names(SUBSETS)) {
    d_lvl <- SUBSETS[[lvl]](d)
    ct <- cor.test(d_lvl$CNNM2, d_lvl$scDRS, method = "pearson")
    stat_rows[[length(stat_rows)+1]] <- data.table(
      GWAS = g, level = lvl, n_cell = nrow(d_lvl),
      r = round(ct$estimate, 3), P = signif(ct$p.value, 3))
  }
}

stats_df <- rbindlist(stat_rows)
fwrite(stats_df, file.path(DAT, "CNNM2_scDRS_correlation_stats.csv"))
options(width = 200)
cat("\n========== CNNM2 × scDRS correlation summary ==========\n")
print(stats_df)

cat("\n[", as.character(Sys.time()), "] DONE\n")
