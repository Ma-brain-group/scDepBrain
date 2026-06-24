# =================================================================
# Per-cell TRS (scDRS norm_score) by brain region (batch proxy)
# 3 GWAS (HW + EUR + DIV, 1000 ctrl)
# Check: is M1C (Allen-M1) highest?
# =================================================================
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(dplyr); library(ggpubr)
})

PP   <- "/mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_scDRS_20260618"
FIG  <- "/mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_HW_v2_figures_20260617/figures"
DAT  <- "/mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_HW_v2_figures_20260617/data"
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)

# ---------- Batch → brain region label ----------
BATCH2REGION <- c(
  "Allen-M1"    = "M1C\n(motor cortex)",
  "ACC"         = "ACC",
  "CER"         = "Cerebellum",
  "GSE97942"    = "Cerebellum\n(GSE97942)",
  "CN"          = "Caudate N.",
  "PRJNA434002" = "ACC+DLPFC\n(BA24/46)",
  "Psychencode" = "DLPFC\n(Psychencode)",
  "GSE126836"   = "GSE126836",
  "GSE140231"   = "GSE140231",
  "Allen-Multi" = "Allen-Multi"
)

# ---------- per-cell metadata (anno, final_anno, batch) ----------
md <- fread("/mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_PGC295_geneset_20260618/data/PGC_MDD_AddModuleScore_per_cell.csv.gz")
setnames(md, 1, "cell_id")
cat("Total cells in metadata:", nrow(md), "\n")
# Need batch too — re-read from h5ad via python (separate step below).
# For now, load batch from a side file we'll create

# Load batch from h5ad via system call
BATCH_CSV <- file.path(DAT, "cell_batch.csv")
if (!file.exists(BATCH_CSV)) {
  cat("[", as.character(Sys.time()), "] Extract batch from h5ad ...\n")
  system(paste0("source /home/may2/miniconda3/etc/profile.d/conda.sh && conda activate pyscenic && ",
                "python -c \"",
                "import anndata as ad, pandas as pd; ",
                "a = ad.read_h5ad('", PP, "/data/MDD_scDRS.h5ad', backed='r'); ",
                "df = a.obs[['batch']].reset_index().rename(columns={'index':'cell_id'}); ",
                "df.to_csv('", BATCH_CSV, "', index=False)",
                "\""))
}
batch <- fread(BATCH_CSV)
cat("Batch rows:", nrow(batch), "\n")
md <- merge(md, batch, by = "cell_id")
md[, region := BATCH2REGION[batch]]
cat("Cells per region:\n")
print(md[, .N, by = region][order(-N)])

# ---------- scDRS per-cell norm_score (TRS) for HW + EUR + DIV ----------
read_scdrs_score <- function(path, gwas) {
  d <- fread(cmd = paste("zcat", path))
  d[, .(cell_id = index, norm_score, mc_pval, GWAS = gwas)]
}
scd <- rbindlist(list(
  read_scdrs_score(paste0(PP, "/scores_1k/HW.full_score.gz"),  "HW"),
  read_scdrs_score(paste0(PP, "/scores_1k/EUR.full_score.gz"), "EUR"),
  read_scdrs_score(paste0(PP, "/scores_1k/DIV.full_score.gz"), "DIV")
))
setnames(scd, "norm_score", "TRS")
cat("scDRS rows:", nrow(scd), "\n")

# ---------- Merge ----------
d <- merge(md[, .(cell_id, anno, final_anno, region)], scd,
           by = "cell_id", allow.cartesian = TRUE)
d[, GWAS := factor(GWAS, levels = c("HW","EUR","DIV"))]
cat("Final merged rows:", nrow(d), "\n")

# Save per-cell merged data (could be large — gzip)
fwrite(d, file.path(DAT, "TRS_by_region_long.csv.gz"))

# ---------- Mean / median TRS per region per GWAS ----------
summ <- d[, .(n_cell = .N,
              mean_TRS = round(mean(TRS, na.rm=TRUE), 3),
              median_TRS = round(median(TRS, na.rm=TRUE), 3),
              pct_sig_5pct = round(mean(mc_pval < 0.05, na.rm=TRUE) * 100, 1)),
          by = .(GWAS, region)]
summ <- summ[order(GWAS, -mean_TRS)]
fwrite(summ, file.path(DAT, "TRS_by_region_summary.csv"))

cat("\n========== TRS by region summary ==========\n")
options(width = 200)
print(summ)

# ---------- Region order (by HW mean TRS) ----------
region_order <- summ[GWAS == "HW"][order(-mean_TRS)][["region"]]
d[, region := factor(region, levels = region_order)]

# ---------- Plot: ALL cells boxplot by region (3 GWAS facet) ----------
p_all <- ggplot(d[!is.na(region)], aes(x = region, y = TRS, fill = region)) +
  geom_boxplot(outlier.size = 0.1, outlier.alpha = 0.2) +
  geom_hline(yintercept = 0, linetype = 2, color = "gray50") +
  facet_wrap(~ GWAS, nrow = 1) +
  labs(x = NULL, y = "scDRS TRS (norm_score)",
       title = "Per-cell MDD trait relevance score by brain region — all cells",
       subtitle = "M1C = Allen-M1 (primary motor cortex)") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 8),
        legend.position = "none",
        strip.text = element_text(face = "bold", size = 12))
ggsave(file.path(FIG, "TRS_by_region_all_cells.pdf"), p_all,
       width = 14, height = 5)
ggsave(file.path(FIG, "TRS_by_region_all_cells.png"), p_all,
       width = 14, height = 5, dpi = 200, bg = "white")

# ---------- Plot: Excitatory neurons only ----------
ex_d <- d[grepl("^Excitatory", anno) & !is.na(region)]
p_ex <- ggplot(ex_d, aes(x = region, y = TRS, fill = region)) +
  geom_boxplot(outlier.size = 0.1, outlier.alpha = 0.2) +
  geom_hline(yintercept = 0, linetype = 2, color = "gray50") +
  facet_wrap(~ GWAS, nrow = 1) +
  labs(x = NULL, y = "scDRS TRS",
       title = "TRS in Excitatory neurons only — by brain region",
       subtitle = sprintf("n_Ex_cells = %d", uniqueN(ex_d$cell_id))) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 8),
        legend.position = "none",
        strip.text = element_text(face = "bold", size = 12))
ggsave(file.path(FIG, "TRS_by_region_excitatory.pdf"), p_ex,
       width = 14, height = 5)
ggsave(file.path(FIG, "TRS_by_region_excitatory.png"), p_ex,
       width = 14, height = 5, dpi = 200, bg = "white")

# ---------- Plot: Inhibitory neurons only ----------
in_d <- d[grepl("^Inhibitory", anno) & !is.na(region)]
p_in <- ggplot(in_d, aes(x = region, y = TRS, fill = region)) +
  geom_boxplot(outlier.size = 0.1, outlier.alpha = 0.2) +
  geom_hline(yintercept = 0, linetype = 2, color = "gray50") +
  facet_wrap(~ GWAS, nrow = 1) +
  labs(x = NULL, y = "scDRS TRS",
       title = "TRS in Inhibitory neurons only — by brain region",
       subtitle = sprintf("n_In_cells = %d", uniqueN(in_d$cell_id))) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 8),
        legend.position = "none",
        strip.text = element_text(face = "bold", size = 12))
ggsave(file.path(FIG, "TRS_by_region_inhibitory.pdf"), p_in,
       width = 14, height = 5)
ggsave(file.path(FIG, "TRS_by_region_inhibitory.png"), p_in,
       width = 14, height = 5, dpi = 200, bg = "white")

cat("\nWrote 3 figures (all / Ex / In) to:", FIG, "\n")
cat("=== DONE ===\n")
