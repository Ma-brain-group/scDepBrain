# =================================================================
# Cross-GWAS scDRS TRS correlation — subtype level
# Compares per-subtype median scDRS scores between GWAS pairs:
#   (1) HW vs EUR
#   (2) HW vs DIV
# Separately for Inhibitory (6) + Excitatory (8) subtypes
# Style: matches user's published Fig (panel e — labeled scatter + dashed reg line)
# =================================================================
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(dplyr); library(ggrepel); library(patchwork)
})

DAT_F <- "/mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_HW_v2_figures_20260617/data"
SCD_F <- "/mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_scDRS_20260618/scores_1k"
FIG   <- "/mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_HW_v2_figures_20260617/figures"

IN_SUBTYPES <- c("In_PVALB","In_LAMP5","In_SST","In_VIP","In_SHANK2","In_CALM1")
EX_SUBTYPES <- c("Ex-L2/4","Ex-L2/3","Ex-L4/6","Ex-L5","Ex-L5/6","Ex-L6","Ex-NRGN","Ex_mix")

# ---------- Load per-cell scDRS + final_anno for all 3 GWAS ----------
cat("[", as.character(Sys.time()), "] Loading data ...\n")
meta <- fread(cmd = paste("zcat", file.path(DAT_F, "scPagwas_HW_metadata.csv.gz")),
              select = c("cell_id","anno","final_anno"))

read_scd <- function(path) {
  d <- fread(cmd = paste("zcat", path), select = c(1, 3))
  setnames(d, c("cell_id","scDRS")); d
}
scd_HW  <- read_scd(file.path(SCD_F, "HW.full_score.gz"))[, GWAS:="HW"]
scd_EUR <- read_scd(file.path(SCD_F, "EUR.full_score.gz"))[, GWAS:="EUR"]
scd_DIV <- read_scd(file.path(SCD_F, "DIV.full_score.gz"))[, GWAS:="DIV"]

all_scd <- rbindlist(list(scd_HW, scd_EUR, scd_DIV))
all_scd <- merge(all_scd, meta, by = "cell_id")
cat("Total cell-GWAS rows:", nrow(all_scd), "\n")

# ---------- Aggregate to subtype median ----------
agg <- all_scd[, .(median_scDRS = median(scDRS, na.rm = TRUE),
                    n_cell = .N), by = .(GWAS, final_anno)]
# Wide format: subtype × GWAS
wide <- dcast(agg, final_anno ~ GWAS, value.var = "median_scDRS")
wide_n <- dcast(agg, final_anno ~ GWAS, value.var = "n_cell")
setnames(wide_n, c("final_anno","DIV_n","EUR_n","HW_n"))
wide <- merge(wide, wide_n, by = "final_anno")
fwrite(wide, file.path(DAT_F, "TRS_scDRS_subtype_medians.csv"))
cat("\nSubtype median scDRS (Wide):\n")
print(wide)

# ---------- Plot function ----------
mk_scatter <- function(d, x_gwas, y_gwas, class_label) {
  d_p <- d[, .(final_anno, x = get(x_gwas), y = get(y_gwas))]
  r <- cor(d_p$x, d_p$y, method = "pearson", use = "complete.obs")
  pval <- cor.test(d_p$x, d_p$y, method = "pearson")$p.value
  ann_text <- sprintf("r = %.2f, P = %s", r,
                      ifelse(pval < 1e-3,
                             format(pval, scientific = TRUE, digits = 2),
                             sprintf("%.4f", pval)))
  cat("  ", class_label, x_gwas, "vs", y_gwas, ":  r =", round(r,3),
      "  P =", signif(pval,3), "\n")

  ggplot(d_p, aes(x = x, y = y)) +
    geom_smooth(method = "lm", se = TRUE, color = "#bfa07a",
                fill = "gray85", alpha = 0.5, linewidth = 0.7, linetype = "dashed") +
    geom_point(size = 4, color = "#2d7d7e", alpha = 0.85) +
    geom_text_repel(aes(label = final_anno), size = 3.7, color = "black",
                    fontface = "bold", box.padding = 0.5, max.overlaps = 30,
                    segment.size = 0.3, segment.color = "gray50") +
    annotate("text", x = -Inf, y = Inf, label = ann_text,
             hjust = -0.18, vjust = 1.5, size = 4.5, fontface = "bold") +
    labs(x = paste("scDRS TRS score (", x_gwas, ", median)", sep=""),
         y = paste("scDRS TRS score (", y_gwas, ", median)", sep=""),
         title = paste0("TRS correlation — ", class_label,
                        " (scDRS ", x_gwas, " vs ", y_gwas, ")")) +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
          panel.grid.major = element_line(color = "gray92", linewidth = 0.3))
}

# ---------- Generate figures ----------
cat("\n========== Cross-GWAS scDRS subtype correlations ==========\n")
for (cls in c("Inhibitory","Excitatory")) {
  subs <- if (cls=="Inhibitory") IN_SUBTYPES else EX_SUBTYPES
  d_sub <- wide[final_anno %in% subs]

  for (pair in list(c("HW","EUR"), c("HW","DIV"))) {
    p <- mk_scatter(d_sub, pair[1], pair[2], paste0(cls, " neuron subtypes"))
    fname_base <- sprintf("TRS_scDRS_acrossGWAS_%s_%s_vs_%s",
                          tolower(substr(cls,1,3)), pair[1], pair[2])
    ggsave(file.path(FIG, paste0(fname_base, ".pdf")), p, width = 6.5, height = 5.5)
    ggsave(file.path(FIG, paste0(fname_base, ".png")), p, width = 6.5, height = 5.5,
           dpi = 300, bg = "white")
  }
}

# ---------- Combined 2x2 figure ----------
panels <- list()
labels <- c("a","b","c","d")
i <- 1
for (cls in c("Inhibitory","Excitatory")) {
  subs <- if (cls=="Inhibitory") IN_SUBTYPES else EX_SUBTYPES
  d_sub <- wide[final_anno %in% subs]
  for (pair in list(c("HW","EUR"), c("HW","DIV"))) {
    panels[[i]] <- mk_scatter(d_sub, pair[1], pair[2], paste0(cls, " neurons"))
    i <- i + 1
  }
}
combined <- (panels[[1]] | panels[[2]]) /
            (panels[[3]] | panels[[4]]) +
  plot_annotation(
    title = "Cross-GWAS scDRS TRS correlation — subtype level",
    subtitle = "Top: Inhibitory subtypes | Bottom: Excitatory subtypes | HW = Howard 2019, EUR = PGC MDD 2025 European, DIV = PGC MDD 2025 trans-ancestry",
    tag_levels = "a",
    theme = theme(plot.title = element_text(face = "bold", size = 14),
                  plot.subtitle = element_text(size = 10, color = "gray35"))
  )

ggsave(file.path(FIG, "TRS_scDRS_acrossGWAS_combined.pdf"), combined,
       width = 13, height = 11)
ggsave(file.path(FIG, "TRS_scDRS_acrossGWAS_combined.png"), combined,
       width = 13, height = 11, dpi = 300, bg = "white")

cat("\nWrote 5 figures to:", FIG, "\n")
cat("\n[", as.character(Sys.time()), "] DONE\n")
