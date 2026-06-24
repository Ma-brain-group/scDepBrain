# =================================================================
# Per-cell TRS (scPagwas + scDRS) by brain region — full pipeline
# 3 GWAS (HW + EUR + DIV)
# Goal: validate published Figure with M1C boxplot
# =================================================================
suppressPackageStartupMessages({
  library(Seurat); library(data.table); library(ggplot2); library(dplyr)
})
utils::assignInNamespace(".Deprecate", function(...) invisible(NULL),
                         ns = "SeuratObject")

FIG <- "/mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_HW_v2_figures_20260617/figures"
DAT <- "/mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_HW_v2_figures_20260617/data"

# ---------- 1) Load scPagwas EUR (smaller than 24G original rds) ----------
# scPagwas output rds inherits all meta cols (including region if present)
cat("[", as.character(Sys.time()), "] Loading scPagwas EUR rds (6.5 GB) ...\n")
obj <- readRDS("/mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_scPagwas_MDD2025_EUR_20260610/MDD_scPagwas_PGC2025_EUR_v2026_v2.rds")
md <- obj@meta.data
cat("Cells:", nrow(md), "\n")
cat("All meta columns:\n"); print(colnames(md))

# Detect brain-region column heuristically
region_cols <- grep("region|brain|location|tissue|area", colnames(md), ignore.case = TRUE, value = TRUE)
cat("\nCandidate region columns:", region_cols, "\n")
for (c in region_cols) {
  cat("  ", c, ":", paste(head(unique(md[[c]]), 25), collapse=", "), "\n")
}

# scPagwas per-cell TRS column
trs_cols <- grep("trs|score|pagwas|relevance", colnames(md), ignore.case = TRUE, value = TRUE)
cat("\nTRS-related columns:", trs_cols, "\n")
saveRDS(md, file.path(DAT, "scPagwas_EUR_metadata.rds"))
fwrite(as.data.table(cbind(cell_id = rownames(md), md)),
       file.path(DAT, "scPagwas_EUR_metadata.csv.gz"))
cat("Saved EUR scPagwas metadata.\n\n")

# ---------- 2) Same for DIV ----------
cat("[", as.character(Sys.time()), "] Loading scPagwas DIV rds ...\n")
obj_div <- readRDS("/mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_scPagwas_MDD2025_TransAnc_20260610/MDD_scPagwas_PGC2025_TransAnc_v2026_v2.rds")
md_div <- obj_div@meta.data
fwrite(as.data.table(cbind(cell_id = rownames(md_div), md_div)),
       file.path(DAT, "scPagwas_DIV_metadata.csv.gz"))
cat("Saved DIV scPagwas metadata.\n\n")

# ---------- 3) Same for HW ----------
cat("[", as.character(Sys.time()), "] Loading scPagwas HW rds ...\n")
obj_hw <- readRDS("/mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_scPagwas_HW_paperparams_20260615/MDD_scPagwas_HW_paperparams.rds")
md_hw <- obj_hw@meta.data
fwrite(as.data.table(cbind(cell_id = rownames(md_hw), md_hw)),
       file.path(DAT, "scPagwas_HW_metadata.csv.gz"))
cat("Saved HW scPagwas metadata.\n")

cat("\n=== DONE — check the 3 metadata files for region column ===\n")
