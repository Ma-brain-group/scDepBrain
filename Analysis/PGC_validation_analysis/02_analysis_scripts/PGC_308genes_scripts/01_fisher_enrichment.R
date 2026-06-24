# =================================================================
# Fisher enrichment: 295 PGC genes ∩ cell-type top 10% specific genes
# 8 major + 20 subtype (using existing MAGMA-pipeline gene sets)
# =================================================================
suppressPackageStartupMessages({ library(dplyr) })

PGC_FILE <- "/mnt/isilon/gandal_lab/mayl/07_scDepBrain/PGC_308_high_confidence_genes.txt"
GS_DIR   <- "/mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_LDSC_PGC2025_20260617/celltype_genesets"
OUT_DIR  <- "/mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_PGC295_geneset_20260618/data"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

pgc <- readLines(PGC_FILE)
pgc <- pgc[nchar(trimws(pgc)) > 0]
pgc <- unique(pgc)
cat("PGC genes:", length(pgc), "\n")

# All genes universe = union of all cell-type top10 lists
all_files <- list.files(GS_DIR, pattern = "\\.geneset\\.txt$", full.names = TRUE)
all_genes <- unique(unlist(lapply(all_files, readLines)))
N <- length(all_genes)
cat("Background universe:", N, "genes (union of all top-10% sets)\n")

# Restrict PGC to genes in universe
pgc_universe <- intersect(pgc, all_genes)
K <- length(pgc_universe)
cat("PGC ∩ universe:", K, "\n")

results <- list()
for (f in all_files) {
  ct <- sub("\\.geneset\\.txt$", "", basename(f))
  level <- ifelse(grepl("^major_", ct), "major", "subtype")
  ct_clean <- sub("^(major|subtype)_", "", ct)

  ct_genes <- readLines(f)
  n  <- length(ct_genes)              # cell-type set size
  k  <- length(intersect(pgc_universe, ct_genes))  # overlap
  # Hypergeometric one-sided P (enrichment)
  p_val <- phyper(k - 1, K, N - K, n, lower.tail = FALSE)
  OR <- (k * (N - K - n + k)) / ((K - k) * (n - k))
  results[[ct]] <- data.frame(
    level = level,
    cell_type = ct_clean,
    n_overlap = k,
    n_celltype = n,
    n_PGC_in_universe = K,
    n_universe = N,
    OR = round(OR, 2),
    P = signif(p_val, 3)
  )
}
df <- do.call(rbind, results)
df$FDR <- signif(p.adjust(df$P, "BH"), 3)
rownames(df) <- NULL

# Save + display
write.csv(df, file.path(OUT_DIR, "fisher_enrichment.csv"), row.names = FALSE)

cat("\n========== MAJOR (8) — Fisher P sorted ==========\n")
print(df %>% filter(level == "major") %>% arrange(P), row.names = FALSE)
cat("\n========== SUBTYPE (20) — Fisher P sorted ==========\n")
print(df %>% filter(level == "subtype") %>% arrange(P), row.names = FALSE)
