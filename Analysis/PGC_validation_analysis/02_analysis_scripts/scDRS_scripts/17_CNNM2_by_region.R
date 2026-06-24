# =================================================================
# CNNM2 expression by brain region — three cell-type levels
# Levels:
#   (1) All cells (n = ~310,833)
#   (2) Inhibitory neurons (n = ~57,453)
#   (3) PVALB+ inhibitory subset (n = ~16,363)
# X-axis: 14 brain regions ordered by median CNNM2 in panel 1
# Output: 3-panel figure + Wilcoxon test M1C vs other regions
# =================================================================
suppressPackageStartupMessages({
  library(Seurat); library(data.table); library(ggplot2); library(dplyr); library(patchwork)
})
utils::assignInNamespace(".Deprecate", function(...) invisible(NULL),
                         ns = "SeuratObject")

OUT_F <- "/mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_HW_v2_figures_20260617"
FIG <- file.path(OUT_F, "figures"); DAT <- file.path(OUT_F, "data")

# ---------- Step 1: Load any scPagwas rds (has full Seurat + tissue.region) ----------
# Use HW since smallest path
RDS_PATH <- "/mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_scPagwas_HW_paperparams_20260615/MDD_scPagwas_HW_paperparams.rds"
cat("[", as.character(Sys.time()), "] Loading scPagwas rds (6.5 GB) ...\n")
obj <- readRDS(RDS_PATH)
DefaultAssay(obj) <- "RNA"
if (inherits(obj[["RNA"]], "Assay5")) {
  obj <- JoinLayers(obj, assay = "RNA")
  obj[["RNA"]] <- as(obj[["RNA"]], "Assay")
}
cat("Cells:", ncol(obj), "\n")

# ---------- Step 2: Get CNNM2 expression + metadata ----------
GENE <- "CNNM2"
if (!GENE %in% rownames(obj)) stop("CNNM2 not in rds")
expr <- FetchData(obj, vars = GENE, slot = "data")
md <- obj@meta.data

d <- data.table(
  cell_id   = rownames(md),
  CNNM2     = expr[[GENE]],
  region    = md$tissue.region,
  anno      = md$anno,
  final_anno = md$final_anno
)
rm(obj); gc()
cat("Region counts:\n"); print(table(d$region))

# ---------- Step 3: Define 14 region order + colors ----------
REGION_ORDER <- c("M1C","S1C","MTG","A1C","PFC","V1C","CTX","DFC","FC",
                  "LA","ACC","SN","CN","CER")
REGION_COLORS <- c("M1C" = "#C5A4FE", "S1C" = "#D58A6C", "MTG" = "#A69333",
                   "A1C" = "#E54B45", "PFC" = "#B89F7E", "V1C" = "#5E2C82",
                   "CTX" = "#9DD2EA", "DFC" = "#7E303C", "FC"  = "#E48887",
                   "LA"  = "#E5A5B5", "ACC" = "#E68C45", "SN"  = "#E66B58",
                   "CN"  = "#3A7FCB", "CER" = "#3F8C3D")

d <- d[region %in% REGION_ORDER]
d[, region := factor(region, levels = REGION_ORDER)]

# ---------- Step 4: Make per-cell-class subsets ----------
SUBSETS <- list(
  "All cells"               = d,
  "Inhibitory neurons"      = d[anno == "Inhibitory neurons"],
  "In_PVALB inhibitory neurons" = d[final_anno == "In_PVALB"]
)

# ---------- Step 5: Stats + boxplot per subset ----------
summary_list <- list()
plot_list <- list()
for (lvl in names(SUBSETS)) {
  d_lvl <- SUBSETS[[lvl]]
  if (nrow(d_lvl) == 0) next
  cat("\n========== ", lvl, " (n =", nrow(d_lvl), ") ==========\n")

  # Median + Wilcoxon M1C vs other
  summ <- d_lvl[, .(n_cell = .N,
                    median_CNNM2 = round(median(CNNM2, na.rm=TRUE), 3),
                    mean_CNNM2 = round(mean(CNNM2, na.rm=TRUE), 3)),
                by = region]
  m1c_expr <- d_lvl[region == "M1C", CNNM2]
  other_expr <- d_lvl[region != "M1C", CNNM2]
  w <- if (length(m1c_expr) > 0 && length(other_expr) > 0) {
    wilcox.test(m1c_expr, other_expr, alternative = "greater")
  } else NULL

  cat("M1C median =", median(m1c_expr, na.rm=TRUE),
      "  vs other median =", median(other_expr, na.rm=TRUE), "\n")
  cat("Wilcoxon M1C vs other (one-sided, M1C > other) P =",
      ifelse(!is.null(w), signif(w$p.value, 3), NA), "\n")
  print(summ[order(-median_CNNM2)])
  summary_list[[lvl]] <- cbind(level = lvl, summ)

  # Boxplot
  ymax <- quantile(d_lvl$CNNM2, 0.99, na.rm = TRUE)
  p <- ggplot(d_lvl, aes(x = region, y = CNNM2, fill = region)) +
    geom_boxplot(outlier.size = 0.08, outlier.alpha = 0.2, lwd = 0.3) +
    scale_fill_manual(values = REGION_COLORS) +
    coord_cartesian(ylim = c(0, ymax * 1.05)) +
    labs(x = NULL, y = "CNNM2 expression (log-normalized)",
         title = paste0(lvl, " (n = ", scales::comma(nrow(d_lvl)), " cells)"),
         subtitle = sprintf("Wilcoxon M1C vs other regions (one-sided): P = %s",
                            ifelse(!is.null(w),
                                   ifelse(w$p.value < 1e-3,
                                          format(w$p.value, scientific = TRUE, digits = 2),
                                          sprintf("%.3f", w$p.value)),
                                   "N/A"))) +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(
            angle = 0, hjust = 0.5, size = 8.5,
            face = ifelse(REGION_ORDER == "M1C", "bold", "plain"),
            color = ifelse(REGION_ORDER == "M1C", "red", "black")),
          legend.position = "none",
          plot.title = element_text(face = "bold", size = 11),
          plot.subtitle = element_text(size = 9, color = "gray30"),
          panel.grid.minor = element_blank())
  plot_list[[lvl]] <- p
}
all_summ <- rbindlist(summary_list)
fwrite(all_summ, file.path(DAT, "CNNM2_by_region_summary.csv"))

# ---------- Step 6: Combine 3 panels vertically ----------
combo <- (plot_list[[1]] / plot_list[[2]] / plot_list[[3]]) +
  plot_annotation(
    title = "CNNM2 expression by brain region — three cell-type levels",
    subtitle = "M1C label highlighted in red | regions ordered as in published TRS figure",
    tag_levels = "a",
    theme = theme(plot.title = element_text(face = "bold", size = 13),
                  plot.subtitle = element_text(size = 10, color = "gray35"))
  )

ggsave(file.path(FIG, "CNNM2_by_region_3levels.pdf"), combo,
       width = 10, height = 10)
ggsave(file.path(FIG, "CNNM2_by_region_3levels.png"), combo,
       width = 10, height = 10, dpi = 300, bg = "white")

# Save also individual figures
for (lvl in names(SUBSETS)) {
  tag <- gsub("[ +]", "_", lvl)
  ggsave(file.path(FIG, paste0("CNNM2_by_region_", tag, ".pdf")), plot_list[[lvl]],
         width = 10, height = 4.5)
  ggsave(file.path(FIG, paste0("CNNM2_by_region_", tag, ".png")), plot_list[[lvl]],
         width = 10, height = 4.5, dpi = 300, bg = "white")
}

cat("\nWrote 4 figures to:", FIG, "\n")
cat("\n[", as.character(Sys.time()), "] DONE\n")
