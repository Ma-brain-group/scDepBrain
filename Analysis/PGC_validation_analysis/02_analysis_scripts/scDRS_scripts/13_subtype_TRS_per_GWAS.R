# =================================================================
# Single-panel scDRS TRS boxplots — per GWAS × per neuronal class
# 6 figures: {HW, EUR, DIV} × {Excitatory, Inhibitory}
# Style: single panel like user's published Fig (matches scPagwas screenshots)
# Each panel sorted by its OWN median TRS (descending)
# =================================================================
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(dplyr)
})

DAT_F <- "/mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_HW_v2_figures_20260617/data"
SCD_F <- "/mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_scDRS_20260618/scores_1k"
FIG   <- "/mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_HW_v2_figures_20260617/figures"

# ---------- Load metadata + scDRS ----------
cat("[", as.character(Sys.time()), "] Loading data ...\n")
read_meta <- function(path) {
  d <- fread(cmd = paste("zcat", path), select = c("cell_id","final_anno"))
  setnames(d, c("cell_id","subtype")); d
}
read_scd <- function(path) {
  d <- fread(cmd = paste("zcat", path), select = c(1, 3))
  setnames(d, c("cell_id","TRS")); d
}

GWAS_META <- list(
  HW  = list(meta = file.path(DAT_F, "scPagwas_HW_metadata.csv.gz"),
             scd  = file.path(SCD_F, "HW.full_score.gz"),
             label = "MDD 2019 EUR (Howard et al.)"),
  EUR = list(meta = file.path(DAT_F, "scPagwas_EUR_metadata.csv.gz"),
             scd  = file.path(SCD_F, "EUR.full_score.gz"),
             label = "MDD 2025 EUR"),
  DIV = list(meta = file.path(DAT_F, "scPagwas_DIV_metadata.csv.gz"),
             scd  = file.path(SCD_F, "DIV.full_score.gz"),
             label = "MDD 2025 trans-ancestry")
)

# ---------- Subtypes + Color palette (matches user's published Fig) ----------
EX_SUBTYPES <- c("Ex-L2/4","Ex-L2/3","Ex-L4/6","Ex-L5","Ex-L5/6","Ex-L6",
                 "Ex-NRGN","Ex_mix")
IN_SUBTYPES <- c("In_PVALB","In_LAMP5","In_SST","In_VIP","In_SHANK2","In_CALM1")

EX_COLORS <- c(
  "Ex-L2/4"  = "#A69333",
  "Ex-L2/3"  = "#D58A6C",
  "Ex-L5"    = "#4F9D4A",
  "Ex-L4/6"  = "#7DBE6F",
  "Ex-L5/6"  = "#3E8E59",
  "Ex-L6"    = "#E54B45",
  "Ex-NRGN"  = "#E66B58",
  "Ex_mix"   = "#9DA9B1"
)
IN_COLORS <- c(
  "In_PVALB"  = "#E5C100",
  "In_LAMP5"  = "#E68C45",
  "In_SST"    = "#5BAEAE",
  "In_VIP"    = "#A69333",
  "In_SHANK2" = "#E54B45",
  "In_CALM1"  = "#3F8C8C"
)

# ---------- Plot one panel ----------
mkpanel <- function(d, subtypes, colors, title_str) {
  d <- d[subtype %in% subtypes]
  # Order by THIS panel's median (descending)
  ord <- d[, .(med = median(TRS, na.rm = TRUE)), by = subtype][order(-med)]$subtype
  d[, subtype := factor(subtype, levels = ord)]

  ggplot(d, aes(x = subtype, y = TRS, fill = subtype)) +
    geom_boxplot(outlier.size = 0.08, outlier.alpha = 0.2, lwd = 0.3) +
    scale_fill_manual(values = colors) +
    labs(x = NULL, y = "scDRS TRS Score", title = title_str) +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 10,
                                      face = "bold"),
          legend.position = "none",
          plot.title = element_text(face = "bold", size = 12),
          panel.grid.minor = element_blank())
}

# ---------- Loop over 3 GWAS × 2 cell classes ----------
for (g in names(GWAS_META)) {
  cfg <- GWAS_META[[g]]
  meta <- read_meta(cfg$meta)
  scd  <- read_scd(cfg$scd)
  d <- merge(meta, scd, by = "cell_id")
  cat(g, ":", nrow(d), "cells\n")

  # Excitatory panel
  p_ex <- mkpanel(d, EX_SUBTYPES, EX_COLORS,
                  paste0("Excitatory neuronal subtype-resolved TRS distribution — ", cfg$label))
  ggsave(file.path(FIG, paste0("subtype_TRS_excitatory_scDRS_", g, ".pdf")),
         p_ex, width = 7, height = 4.5)
  ggsave(file.path(FIG, paste0("subtype_TRS_excitatory_scDRS_", g, ".png")),
         p_ex, width = 7, height = 4.5, dpi = 200, bg = "white")

  # Inhibitory panel
  p_in <- mkpanel(d, IN_SUBTYPES, IN_COLORS,
                  paste0("Inhibitory neuronal subtype-resolved TRS distribution — ", cfg$label))
  ggsave(file.path(FIG, paste0("subtype_TRS_inhibitory_scDRS_", g, ".pdf")),
         p_in, width = 6, height = 4.5)
  ggsave(file.path(FIG, paste0("subtype_TRS_inhibitory_scDRS_", g, ".png")),
         p_in, width = 6, height = 4.5, dpi = 200, bg = "white")

  cat("  Wrote 2 figures for", g, "\n")
}

cat("\n[", as.character(Sys.time()), "] DONE — 6 single-panel figures generated\n")
system(paste("ls -lh", FIG, "| grep subtype_TRS_.*scDRS_(HW|EUR|DIV) | head"))
