# =================================================================
# Subtype-resolved scDRS TRS distribution boxplots
# Style: mimics published scPagwas TRS boxplots (Excitatory + Inhibitory)
# 3 GWAS (HW + EUR + DIV) as facets
# Subtype order: sorted by HW scDRS median (descending), held constant
# =================================================================
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(dplyr); library(patchwork)
})

DAT_F <- "/mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_HW_v2_figures_20260617/data"
SCD_F <- "/mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_scDRS_20260618/scores_1k"
FIG   <- "/mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_HW_v2_figures_20260617/figures"

# ---------- Load metadata ----------
cat("[", as.character(Sys.time()), "] Loading metadata ...\n")
read_meta <- function(path) {
  d <- fread(cmd = paste("zcat", path), select = c("cell_id","final_anno"))
  setnames(d, c("cell_id","subtype"))
  d
}
meta_HW  <- read_meta(file.path(DAT_F, "scPagwas_HW_metadata.csv.gz"))
meta_EUR <- read_meta(file.path(DAT_F, "scPagwas_EUR_metadata.csv.gz"))
meta_DIV <- read_meta(file.path(DAT_F, "scPagwas_DIV_metadata.csv.gz"))

# ---------- Load scDRS norm_score (1000 ctrl) ----------
read_scd <- function(path) {
  d <- fread(cmd = paste("zcat", path), select = c(1, 3))
  setnames(d, c("cell_id","TRS"))
  d
}
scd_HW  <- read_scd(file.path(SCD_F, "HW.full_score.gz"))
scd_EUR <- read_scd(file.path(SCD_F, "EUR.full_score.gz"))
scd_DIV <- read_scd(file.path(SCD_F, "DIV.full_score.gz"))

merge_one <- function(meta, scd, gwas) {
  d <- merge(meta, scd, by = "cell_id")
  d[, GWAS := gwas]; d
}
all <- rbindlist(list(
  merge_one(meta_HW,  scd_HW,  "HW"),
  merge_one(meta_EUR, scd_EUR, "EUR"),
  merge_one(meta_DIV, scd_DIV, "DIV")
))
all[, GWAS := factor(GWAS, levels = c("HW","EUR","DIV"))]
cat("Total rows:", nrow(all), "\n")

# ---------- Define subtypes ----------
EX_SUBTYPES <- c("Ex-L2/4","Ex-L2/3","Ex-L4/6","Ex-L5","Ex-L5/6","Ex-L6",
                 "Ex-NRGN","Ex_mix")
IN_SUBTYPES <- c("In_PVALB","In_LAMP5","In_SST","In_VIP","In_SHANK2","In_CALM1")

# Color palette (matches user's published figure)
EX_COLORS <- c(
  "Ex-L2/4"  = "#A69333",   # olive
  "Ex-L2/3"  = "#D58A6C",   # salmon
  "Ex-L5"    = "#4F9D4A",   # green
  "Ex-L4/6"  = "#7DBE6F",   # light green
  "Ex-L5/6"  = "#3E8E59",   # darker green
  "Ex-L6"    = "#E54B45",   # red
  "Ex-NRGN"  = "#E66B58",   # coral
  "Ex_mix"   = "#9DA9B1"    # gray-blue
)
IN_COLORS <- c(
  "In_PVALB"  = "#E5C100",   # yellow
  "In_LAMP5"  = "#E68C45",   # orange
  "In_SST"    = "#5BAEAE",   # teal
  "In_VIP"    = "#A69333",   # olive
  "In_SHANK2" = "#E54B45",   # red
  "In_CALM1"  = "#3F8C8C"    # blue-teal
)

# ---------- Excitatory ----------
ex_d <- all[subtype %in% EX_SUBTYPES]
# order by HW median descending
ex_order <- ex_d[GWAS == "HW", .(med = median(TRS, na.rm = TRUE)), by = subtype][order(-med)]$subtype
ex_d[, subtype := factor(subtype, levels = ex_order)]
cat("\nExcitatory order (by HW scDRS median):\n"); print(ex_order)

p_ex <- ggplot(ex_d, aes(x = subtype, y = TRS, fill = subtype)) +
  geom_boxplot(outlier.size = 0.08, outlier.alpha = 0.2, lwd = 0.3) +
  scale_fill_manual(values = EX_COLORS) +
  facet_wrap(~ GWAS, ncol = 1, scales = "free_y",
             labeller = labeller(GWAS = c(
               "HW" = "MDD 2019 EUR (Howard et al.)",
               "EUR" = "MDD 2025 EUR",
               "DIV" = "MDD 2025 trans-ancestry"))) +
  labs(x = NULL, y = "scDRS TRS Score (norm_score)",
       title = "Excitatory neuronal subtype-resolved TRS distribution",
       subtitle = "Per-cell scDRS norm_score, 1,000 control gene sets | subtypes ordered by HW median") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 10),
        legend.position = "none",
        strip.text = element_text(face = "bold", size = 11),
        strip.background = element_rect(fill = "gray95"),
        panel.grid.minor = element_blank())

ggsave(file.path(FIG, "subtype_TRS_excitatory_scDRS.pdf"), p_ex,
       width = 9, height = 9)
ggsave(file.path(FIG, "subtype_TRS_excitatory_scDRS.png"), p_ex,
       width = 9, height = 9, dpi = 200, bg = "white")

# ---------- Inhibitory ----------
in_d <- all[subtype %in% IN_SUBTYPES]
in_order <- in_d[GWAS == "HW", .(med = median(TRS, na.rm = TRUE)), by = subtype][order(-med)]$subtype
in_d[, subtype := factor(subtype, levels = in_order)]
cat("\nInhibitory order (by HW scDRS median):\n"); print(in_order)

p_in <- ggplot(in_d, aes(x = subtype, y = TRS, fill = subtype)) +
  geom_boxplot(outlier.size = 0.08, outlier.alpha = 0.2, lwd = 0.3) +
  scale_fill_manual(values = IN_COLORS) +
  facet_wrap(~ GWAS, ncol = 1, scales = "free_y",
             labeller = labeller(GWAS = c(
               "HW" = "MDD 2019 EUR (Howard et al.)",
               "EUR" = "MDD 2025 EUR",
               "DIV" = "MDD 2025 trans-ancestry"))) +
  labs(x = NULL, y = "scDRS TRS Score (norm_score)",
       title = "Inhibitory neuronal subtype-resolved TRS distribution",
       subtitle = "Per-cell scDRS norm_score, 1,000 control gene sets | subtypes ordered by HW median") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 10),
        legend.position = "none",
        strip.text = element_text(face = "bold", size = 11),
        strip.background = element_rect(fill = "gray95"),
        panel.grid.minor = element_blank())

ggsave(file.path(FIG, "subtype_TRS_inhibitory_scDRS.pdf"), p_in,
       width = 8, height = 9)
ggsave(file.path(FIG, "subtype_TRS_inhibitory_scDRS.png"), p_in,
       width = 8, height = 9, dpi = 200, bg = "white")

# ---------- Summary table ----------
summ_ex <- ex_d[, .(n_cell = .N,
                    median_TRS = round(median(TRS, na.rm = TRUE), 3),
                    mean_TRS   = round(mean(TRS, na.rm = TRUE), 3)),
                by = .(GWAS, subtype)]
summ_in <- in_d[, .(n_cell = .N,
                    median_TRS = round(median(TRS, na.rm = TRUE), 3),
                    mean_TRS   = round(mean(TRS, na.rm = TRUE), 3)),
                by = .(GWAS, subtype)]
fwrite(rbind(cbind(class = "Ex", summ_ex), cbind(class = "In", summ_in)),
       file.path(DAT_F, "subtype_TRS_scDRS_summary.csv"))

options(width = 200)
cat("\n========== Excitatory subtype TRS ==========\n")
print(summ_ex[order(GWAS, -median_TRS)])
cat("\n========== Inhibitory subtype TRS ==========\n")
print(summ_in[order(GWAS, -median_TRS)])

cat("\nWrote:\n",
    file.path(FIG, "subtype_TRS_excitatory_scDRS.pdf"), "\n",
    file.path(FIG, "subtype_TRS_inhibitory_scDRS.pdf"), "\n", sep = "")
cat("\n[", as.character(Sys.time()), "] DONE\n")
