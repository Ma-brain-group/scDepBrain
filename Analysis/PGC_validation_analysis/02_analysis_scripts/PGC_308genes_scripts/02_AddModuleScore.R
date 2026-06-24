# =================================================================
# AddModuleScore for PGC 295 high-confidence MDD genes
# Per-cell score + boxplot by anno (8 major) and final_anno (20 subtype)
# =================================================================
suppressPackageStartupMessages({
  library(Seurat); library(ggplot2); library(ggpubr); library(dplyr)
})
utils::assignInNamespace(".Deprecate", function(...) invisible(NULL),
                         ns = "SeuratObject")

RDS <- "/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/integrated_heathy_brain_cells_short_read/MDD_singlecell_data_reannotation.rds"
PGC_FILE <- "/mnt/isilon/gandal_lab/mayl/07_scDepBrain/PGC_308_high_confidence_genes.txt"
OUT <- "/mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_PGC295_geneset_20260618"
FIG <- file.path(OUT, "figures"); DAT <- file.path(OUT, "data")

pgc <- readLines(PGC_FILE); pgc <- pgc[nchar(trimws(pgc)) > 0]; pgc <- unique(pgc)
cat("PGC gene list:", length(pgc), "\n")

cat("[", as.character(Sys.time()), "] Loading 24GB rds ...\n")
t0 <- Sys.time()
obj <- readRDS(RDS)
DefaultAssay(obj) <- "RNA"
if (inherits(obj[["RNA"]], "Assay5")) {
  obj <- JoinLayers(obj, assay = "RNA")
  obj[["RNA"]] <- as(obj[["RNA"]], "Assay")
}
if ("SCT" %in% Assays(obj)) { obj[["SCT"]] <- NULL; gc() }
cat(sprintf("Loaded in %.1f min\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))

# Filter PGC genes to those in our scRNA
pgc_in_data <- intersect(pgc, rownames(obj))
cat("PGC genes in scRNA data:", length(pgc_in_data), "/", length(pgc), "\n")

# AddModuleScore (default: 100 ctrl genes per bin, 24 bins)
cat("[", as.character(Sys.time()), "] AddModuleScore ...\n")
set.seed(42)
obj <- AddModuleScore(obj, features = list(PGC_MDD = pgc_in_data),
                      ctrl = 100, nbin = 24, name = "PGC_MDD_score")
score_col <- grep("^PGC_MDD_score", colnames(obj@meta.data), value = TRUE)[1]
cat("Score column:", score_col, "\n")

# Save per-cell scores + cell type info
md <- obj@meta.data
md_save <- md[, c("anno","final_anno", score_col)]
colnames(md_save)[3] <- "PGC_MDD_score"
write.csv(md_save, file.path(DAT, "PGC_MDD_AddModuleScore_per_cell.csv.gz"),
          row.names = TRUE)

# Summary stats per cell type
for (level in c("anno","final_anno")) {
  agg <- md_save %>%
    group_by(.data[[level]]) %>%
    summarise(n = n(),
              mean_score = mean(PGC_MDD_score, na.rm = TRUE),
              median_score = median(PGC_MDD_score, na.rm = TRUE),
              sd_score = sd(PGC_MDD_score, na.rm = TRUE)) %>%
    arrange(desc(mean_score))
  write.csv(agg, file.path(DAT, paste0("PGC_MDD_AddModuleScore_summary_", level, ".csv")),
            row.names = FALSE)
  cat("\n=== ", level, " — top 5 by mean score ===\n")
  print(head(agg, 5))
}

# ---------- Boxplot (anno = 8 major) ----------
p1 <- ggplot(md_save, aes(x = reorder(anno, PGC_MDD_score, FUN = median),
                          y = PGC_MDD_score, fill = anno)) +
  geom_boxplot(outlier.size = 0.2, outlier.alpha = 0.3) +
  geom_hline(yintercept = 0, linetype = 2, color = "gray50") +
  scale_fill_brewer(palette = "Set2") +
  labs(x = NULL,
       y = "PGC-MDD-295 module score",
       title = "Per-cell expression of PGC 295 MDD genes by major cell type",
       subtitle = sprintf("n_cells = %d | n_genes used = %d",
                          nrow(md_save), length(pgc_in_data))) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        legend.position = "none")

ggsave(file.path(FIG, "PGC_MDD_AddModuleScore_major.pdf"), p1,
       width = 8, height = 5)
ggsave(file.path(FIG, "PGC_MDD_AddModuleScore_major.png"), p1,
       width = 8, height = 5, dpi = 200, bg = "white")

# ---------- Boxplot (final_anno = 20 subtype) ----------
p2 <- ggplot(md_save, aes(x = reorder(final_anno, PGC_MDD_score, FUN = median),
                          y = PGC_MDD_score, fill = final_anno)) +
  geom_boxplot(outlier.size = 0.2, outlier.alpha = 0.3) +
  geom_hline(yintercept = 0, linetype = 2, color = "gray50") +
  labs(x = NULL,
       y = "PGC-MDD-295 module score",
       title = "Per-cell expression of PGC 295 MDD genes by subtype",
       subtitle = sprintf("n_cells = %d | n_genes used = %d",
                          nrow(md_save), length(pgc_in_data))) +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        legend.position = "none")

ggsave(file.path(FIG, "PGC_MDD_AddModuleScore_subtype.pdf"), p2,
       width = 12, height = 5)
ggsave(file.path(FIG, "PGC_MDD_AddModuleScore_subtype.png"), p2,
       width = 12, height = 5, dpi = 200, bg = "white")

cat("\nWrote figures to:", FIG, "\n")
cat("=== DONE ===\n")
