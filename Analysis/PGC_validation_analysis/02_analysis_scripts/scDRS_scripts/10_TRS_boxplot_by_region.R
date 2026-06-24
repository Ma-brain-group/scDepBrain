# =================================================================
# Per-cell TRS by brain region — replicate published paper Figure
# 14 regions: M1C, MTG, S1C, A1C, CTX, V1C, DFC, FC, PFC, LA, ACC, SN, CN, CER
# 3 GWAS × 2 methods (scPagwas TRS + scDRS norm_score)
# =================================================================
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(dplyr); library(patchwork)
})

DAT_F <- "/mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_HW_v2_figures_20260617/data"
SCD_F <- "/mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_scDRS_20260618/scores_1k"
FIG   <- "/mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_HW_v2_figures_20260617/figures"

# Region order from published paper Figure (highest TRS → lowest)
REGION_ORDER <- c("M1C","MTG","S1C","A1C","CTX","V1C","DFC","FC","PFC",
                  "LA","ACC","SN","CN","CER")
REGION_COLORS <- c("M1C"="#C5A4FE","MTG"="#A69333","S1C"="#D58A6C","A1C"="#E54B45",
                   "CTX"="#9DD2EA","V1C"="#5E2C82","DFC"="#7E303C","FC"="#E48887",
                   "PFC"="#B89F7E","LA"="#E5A5B5","ACC"="#E68C45","SN"="#E66B58",
                   "CN"="#3A7FCB","CER"="#3F8C3D")

# ---------- Load scPagwas per-cell TRS + region ----------
read_scp_meta <- function(path) {
  d <- fread(cmd = paste("zcat", path),
             select = c("cell_id","tissue.region",
                        "scPagwas.TRS.Score","anno","final_anno"))
  setnames(d, c("cell_id","region","scPagwas_TRS","anno","final_anno"))
  d
}

cat("[", as.character(Sys.time()), "] Loading 3 scPagwas csv.gz ...\n")
meta_HW  <- read_scp_meta(file.path(DAT_F, "scPagwas_HW_metadata.csv.gz"))
meta_EUR <- read_scp_meta(file.path(DAT_F, "scPagwas_EUR_metadata.csv.gz"))
meta_DIV <- read_scp_meta(file.path(DAT_F, "scPagwas_DIV_metadata.csv.gz"))
cat("Rows: HW=", nrow(meta_HW), " EUR=", nrow(meta_EUR), " DIV=", nrow(meta_DIV), "\n")

# ---------- Load scDRS per-cell norm_score ----------
# First column has empty header (cell_id), so use V1 and norm_score
read_scd <- function(path) {
  # Read just headers to find norm_score column position
  hd <- fread(cmd = paste("zcat", path, "| head -1"), header = TRUE)
  # First col empty header named V1 by fread; norm_score is col 3
  d <- fread(cmd = paste("zcat", path), select = c(1, 3))
  setnames(d, c("cell_id","scDRS_TRS"))
  d
}
cat("[", as.character(Sys.time()), "] Loading 3 scDRS scores ...\n")
scd_HW  <- read_scd(file.path(SCD_F, "HW.full_score.gz"))
scd_EUR <- read_scd(file.path(SCD_F, "EUR.full_score.gz"))
scd_DIV <- read_scd(file.path(SCD_F, "DIV.full_score.gz"))

# ---------- Merge ----------
merge_one <- function(meta, scd, gwas) {
  d <- merge(meta, scd, by = "cell_id")
  d[, GWAS := gwas]
  d
}
hw  <- merge_one(meta_HW,  scd_HW,  "HW")
eur <- merge_one(meta_EUR, scd_EUR, "EUR")
div <- merge_one(meta_DIV, scd_DIV, "DIV")
all <- rbindlist(list(hw, eur, div))
all <- all[region %in% REGION_ORDER]
all[, region := factor(region, levels = REGION_ORDER)]
all[, GWAS := factor(GWAS, levels = c("HW","EUR","DIV"))]
cat("Merged rows:", nrow(all), "\n")

# ---------- Summary ----------
summ <- all[, .(n_cell = .N,
                scPagwas_median = round(median(scPagwas_TRS, na.rm=TRUE), 3),
                scPagwas_mean   = round(mean(scPagwas_TRS, na.rm=TRUE), 3),
                scDRS_median = round(median(scDRS_TRS, na.rm=TRUE), 3),
                scDRS_mean   = round(mean(scDRS_TRS, na.rm=TRUE), 3)),
            by = .(GWAS, region)]
summ <- summ[order(GWAS, -scPagwas_median)]
fwrite(summ, file.path(DAT_F, "TRS_by_region_summary.csv"))
options(width = 200)
cat("\n========== TRS by region summary ==========\n")
print(summ)

# ---------- Plot helper ----------
mkplot <- function(d, y, title_str, ylim = NULL) {
  p <- ggplot(d, aes(x = region, y = .data[[y]], fill = region)) +
    geom_boxplot(outlier.size = 0.1, outlier.alpha = 0.2, lwd = 0.3) +
    scale_fill_manual(values = REGION_COLORS) +
    labs(x = "Brain region", y = title_str, title = NULL) +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 9),
          legend.position = "none",
          strip.text = element_text(face = "bold", size = 12),
          panel.grid.minor = element_blank())
  if (!is.null(ylim)) p <- p + coord_cartesian(ylim = ylim)
  p
}

# ---------- Plot 1: scPagwas TRS × 3 GWAS ----------
p_scp <- ggplot(all, aes(x = region, y = scPagwas_TRS, fill = region)) +
  geom_boxplot(outlier.size = 0.1, outlier.alpha = 0.2, lwd = 0.3) +
  scale_fill_manual(values = REGION_COLORS) +
  facet_wrap(~ GWAS, ncol = 1, scales = "free_y") +
  labs(x = "Brain region", y = "scPagwas TRS Score",
       title = "scPagwas per-cell TRS by brain region — All cells (3 GWAS)",
       subtitle = "Replicates published Figure | regions ordered by HW median TRS (M1C → CER)") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 9),
        legend.position = "none",
        strip.text = element_text(face = "bold", size = 12))

ggsave(file.path(FIG, "TRS_scPagwas_by_region.pdf"), p_scp,
       width = 12, height = 10)
ggsave(file.path(FIG, "TRS_scPagwas_by_region.png"), p_scp,
       width = 12, height = 10, dpi = 200, bg = "white")

# ---------- Plot 2: scDRS TRS × 3 GWAS ----------
p_scd <- ggplot(all, aes(x = region, y = scDRS_TRS, fill = region)) +
  geom_boxplot(outlier.size = 0.1, outlier.alpha = 0.2, lwd = 0.3) +
  scale_fill_manual(values = REGION_COLORS) +
  facet_wrap(~ GWAS, ncol = 1, scales = "free_y") +
  labs(x = "Brain region", y = "scDRS TRS (norm_score)",
       title = "scDRS per-cell TRS by brain region — All cells (3 GWAS)",
       subtitle = "1000 control gene sets | regions ordered by HW scPagwas median TRS") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 9),
        legend.position = "none",
        strip.text = element_text(face = "bold", size = 12))

ggsave(file.path(FIG, "TRS_scDRS_by_region.pdf"), p_scd,
       width = 12, height = 10)
ggsave(file.path(FIG, "TRS_scDRS_by_region.png"), p_scd,
       width = 12, height = 10, dpi = 200, bg = "white")

# ---------- Plot 3: HW scPagwas + HW scDRS — side by side (validate paper) ----------
hw_only <- all[GWAS == "HW"]
p_hw_scp <- mkplot(hw_only, "scPagwas_TRS", "scPagwas TRS Score") +
  ggtitle("HW (Howard 2019) — scPagwas TRS")
p_hw_scd <- mkplot(hw_only, "scDRS_TRS", "scDRS TRS (norm_score)") +
  ggtitle("HW (Howard 2019) — scDRS TRS")
p_validation <- p_hw_scp / p_hw_scd
ggsave(file.path(FIG, "TRS_HW_validation_paper.pdf"), p_validation,
       width = 12, height = 8)
ggsave(file.path(FIG, "TRS_HW_validation_paper.png"), p_validation,
       width = 12, height = 8, dpi = 200, bg = "white")

cat("\nWrote:\n",
    file.path(FIG, "TRS_scPagwas_by_region.pdf"), "\n",
    file.path(FIG, "TRS_scDRS_by_region.pdf"), "\n",
    file.path(FIG, "TRS_HW_validation_paper.pdf"), "\n", sep = "")
cat("\n[", as.character(Sys.time()), "] DONE\n")
