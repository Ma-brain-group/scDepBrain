# =================================================================
# scDRS TRS by brain region, per cell-type subset (paper Figure style)
# 4 panels: Inhibitory, In_PVALB, Excitatory, Ex-L2/4
# 3 GWAS (HW, EUR, DIV)
# Region order matches user's published Figure (M1C first, highlighted)
# =================================================================
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(dplyr); library(patchwork)
})

DAT_F <- "/mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_HW_v2_figures_20260617/data"
SCD_F <- "/mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_scDRS_20260618/scores_1k"
FIG   <- "/mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_HW_v2_figures_20260617/figures"

# Region order from published paper Figure (M1C first)
REGION_ORDER <- c("M1C","S1C","MTG","A1C","PFC","V1C","CTX","DFC","FC",
                  "LA","ACC","SN","CN","CER")
# Color scheme matching paper (light purple M1C → green CER)
REGION_COLORS <- c("M1C" = "#C5A4FE", "S1C" = "#D58A6C", "MTG" = "#A69333",
                   "A1C" = "#E54B45", "PFC" = "#B89F7E", "V1C" = "#5E2C82",
                   "CTX" = "#9DD2EA", "DFC" = "#7E303C", "FC"  = "#E48887",
                   "LA"  = "#E5A5B5", "ACC" = "#E68C45", "SN"  = "#E66B58",
                   "CN"  = "#3A7FCB", "CER" = "#3F8C3D")

# ---------- Load metadata ----------
cat("[", as.character(Sys.time()), "] Loading scPagwas metadata (region + anno + final_anno) ...\n")
read_meta <- function(path) {
  d <- fread(cmd = paste("zcat", path),
             select = c("cell_id","tissue.region","anno","final_anno"))
  setnames(d, c("cell_id","region","anno","final_anno"))
  d
}
meta_HW  <- read_meta(file.path(DAT_F, "scPagwas_HW_metadata.csv.gz"))
meta_EUR <- read_meta(file.path(DAT_F, "scPagwas_EUR_metadata.csv.gz"))
meta_DIV <- read_meta(file.path(DAT_F, "scPagwas_DIV_metadata.csv.gz"))

# ---------- Load scDRS per-cell norm_score ----------
read_scd <- function(path) {
  d <- fread(cmd = paste("zcat", path), select = c(1, 3))
  setnames(d, c("cell_id","scDRS_TRS"))
  d
}
cat("[", as.character(Sys.time()), "] Loading 3 scDRS scores ...\n")
scd_HW  <- read_scd(file.path(SCD_F, "HW.full_score.gz"))
scd_EUR <- read_scd(file.path(SCD_F, "EUR.full_score.gz"))
scd_DIV <- read_scd(file.path(SCD_F, "DIV.full_score.gz"))

# Merge
merge_one <- function(meta, scd, gwas) {
  d <- merge(meta, scd, by = "cell_id")
  d[, GWAS := gwas]
  d
}
all <- rbindlist(list(
  merge_one(meta_HW,  scd_HW,  "HW"),
  merge_one(meta_EUR, scd_EUR, "EUR"),
  merge_one(meta_DIV, scd_DIV, "DIV")
))
all <- all[region %in% REGION_ORDER]
all[, region := factor(region, levels = REGION_ORDER)]
all[, GWAS := factor(GWAS, levels = c("HW","EUR","DIV"))]
cat("Merged rows:", nrow(all), "\n")

# ---------- Define 4 cell-type subsets ----------
CELL_SUBSETS <- list(
  "Inhibitory neurons" = all[anno == "Inhibitory neurons"],
  "In_PVALB"           = all[final_anno == "In_PVALB"],
  "Excitatory neurons" = all[anno == "Excitatory neurons"],
  "Ex-L2/4"            = all[final_anno == "Ex-L2/4"]
)

cat("Cells per subset:\n")
for (n in names(CELL_SUBSETS)) {
  cat("  ", n, ":", nrow(CELL_SUBSETS[[n]]) / 3, "cells per GWAS\n")
}

# ---------- Plot one panel ----------
mkpanel <- function(d, title_str) {
  ymax <- quantile(d$scDRS_TRS, 0.995, na.rm = TRUE)
  ymin <- quantile(d$scDRS_TRS, 0.005, na.rm = TRUE)
  ggplot(d, aes(x = region, y = scDRS_TRS, fill = region)) +
    geom_boxplot(outlier.size = 0.05, outlier.alpha = 0.2, lwd = 0.3) +
    scale_fill_manual(values = REGION_COLORS) +
    coord_cartesian(ylim = c(ymin, ymax)) +
    labs(x = "Brain region", y = "scDRS TRS", title = title_str) +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(
              angle = 0, hjust = 0.5, size = 8,
              face = ifelse(REGION_ORDER == "M1C", "bold", "plain"),
              color = ifelse(REGION_ORDER == "M1C", "red", "black")),
          legend.position = "none",
          strip.text = element_text(face = "bold", size = 11),
          plot.title = element_text(face = "bold", size = 12,
                                     hjust = 0.5),
          panel.grid.minor = element_blank())
}

# ---------- ONE FIGURE PER GWAS: 2x2 grid (a / b / c / d) ----------
for (g in c("HW","EUR","DIV")) {
  panels <- lapply(names(CELL_SUBSETS), function(n) {
    sub_d <- CELL_SUBSETS[[n]][GWAS == g]
    if (nrow(sub_d) == 0) return(NULL)
    p <- mkpanel(sub_d, n) + facet_wrap(~ GWAS)
    # Drop facet strip (since we have title)
    p + theme(strip.background = element_blank(),
              strip.text = element_blank())
  })
  combo <- (panels[[1]] | panels[[2]]) / (panels[[3]] | panels[[4]]) +
    plot_annotation(
      title = paste0("scDRS TRS by brain region — ", g),
      subtitle = "M1C highlighted in red | regions ordered as in published Figure",
      tag_levels = "a"
    )
  ggsave(file.path(FIG, paste0("TRS_scDRS_by_region_celltype_", g, ".pdf")),
         combo, width = 12, height = 9)
  ggsave(file.path(FIG, paste0("TRS_scDRS_by_region_celltype_", g, ".png")),
         combo, width = 12, height = 9, dpi = 200, bg = "white")
  cat("Wrote panel for", g, "\n")
}

# ---------- COMBINED FIGURE: 3 GWAS x 4 cell types grid ----------
plot_df <- rbindlist(lapply(names(CELL_SUBSETS), function(n) {
  CELL_SUBSETS[[n]][, cell_subset := n]
}))
plot_df[, cell_subset := factor(cell_subset, levels = names(CELL_SUBSETS))]

p_combined <- ggplot(plot_df, aes(x = region, y = scDRS_TRS, fill = region)) +
  geom_boxplot(outlier.size = 0.05, outlier.alpha = 0.2, lwd = 0.3) +
  scale_fill_manual(values = REGION_COLORS) +
  facet_grid(GWAS ~ cell_subset, scales = "free_y") +
  coord_cartesian(ylim = quantile(plot_df$scDRS_TRS, c(0.005, 0.995),
                                  na.rm = TRUE)) +
  labs(x = "Brain region", y = "scDRS TRS",
       title = "scDRS TRS by brain region — cell-type subsets across 3 GWAS",
       subtitle = "Cell-type panels: Inhibitory / In_PVALB / Excitatory / Ex-L2/4 | M1C in red") +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(
            angle = 0, hjust = 0.5, size = 7,
            face = ifelse(REGION_ORDER == "M1C", "bold", "plain"),
            color = ifelse(REGION_ORDER == "M1C", "red", "black")),
        legend.position = "none",
        strip.text = element_text(face = "bold", size = 11),
        strip.background = element_rect(fill = "gray95"),
        panel.grid.minor = element_blank())

ggsave(file.path(FIG, "TRS_scDRS_by_region_celltype_3GWAS_grid.pdf"),
       p_combined, width = 16, height = 10)
ggsave(file.path(FIG, "TRS_scDRS_by_region_celltype_3GWAS_grid.png"),
       p_combined, width = 16, height = 10, dpi = 200, bg = "white")
cat("\nWrote combined grid figure\n")

# ---------- Summary table per (GWAS, cell_subset, region) ----------
summ <- plot_df[, .(n_cell = .N,
                    median_TRS = round(median(scDRS_TRS, na.rm = TRUE), 3),
                    mean_TRS   = round(mean(scDRS_TRS, na.rm = TRUE), 3)),
                by = .(GWAS, cell_subset, region)]
summ <- summ[order(GWAS, cell_subset, -median_TRS)]
fwrite(summ, file.path(DAT_F, "TRS_scDRS_by_region_celltype_summary.csv"))

options(width = 200)
cat("\n========== Top 3 regions per (GWAS, cell subset) ==========\n")
top3 <- summ[, .SD[1:3], by = .(GWAS, cell_subset)]
print(top3)

cat("\n=== DONE ===\n")
