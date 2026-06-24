# ============================================================
# HW paper-params: Step 1.7 finalize, matching paper code style:
#   - 2 features (top + down) NOT 3 (no upTRS)
#   - n_topgenes = 1000 in Get_CorrectBg_p (paper used 1000)
#   - PCC via Corr_Random (scGet_PCC would OOM on 310k cells)
# Includes the same patches as v2 (fast row variance + fast empirical p)
# ============================================================
suppressPackageStartupMessages({
  library(Seurat); library(scPagwas); library(Matrix)
})

# Compat patches
if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
  RhpcBLASctl::blas_set_num_threads(1)
}
options(default.nproc.blas = 1)
options(bigstatsr.check.parallel.blas = FALSE)
suppressMessages({
  utils::assignInNamespace(".Deprecate",
                           function(...) invisible(NULL),
                           ns = "SeuratObject")
})

# Speed patches
patched_Get_CorrectBg_p <- function(Single_data, scPagwas.TRS.Score,
                                    iters_singlecell, n_topgenes,
                                    scPagwas_topgenes, assay = "RNA") {
  gene_matrix <- GetAssayData(Single_data, slot = "data", assay = assay)
  mat_ctrl_raw_score <- matrix(0, nrow = ncol(gene_matrix),
                               ncol = iters_singlecell)
  dic_ctrl_list <- list()
  pb <- txtProgressBar(style = 3)
  for (i in 1:iters_singlecell) {
    set.seed(i)
    dic_ctrl_list[[i]] <- sample(rownames(Single_data), n_topgenes)
    Single_data <- Seurat::AddModuleScore(Single_data, assay = assay,
                                          list(dic_ctrl_list[[i]]),
                                          name = c("contr_genes"))
    mat_ctrl_raw_score[, i] <- Single_data$contr_genes1
    Single_data$contr_genes1 <- NULL
    setTxtProgressBar(pb, i/iters_singlecell)
  }
  close(pb)
  genes <- intersect(rownames(Single_data), rownames(gene_matrix))
  scPagwas_topgenes <- intersect(scPagwas_topgenes, genes)
  gene_matrix <- gene_matrix[genes, ]
  cat("\n[PATCH 1] fast row variance ...\n")
  t0 <- Sys.time()
  n <- ncol(gene_matrix)
  mu <- Matrix::rowMeans(gene_matrix)
  ss <- Matrix::rowSums(gene_matrix^2)
  rv <- as.numeric((ss - n * mu^2) / (n - 1))
  cat(sprintf("[PATCH 1] done in %.1f sec\n",
              as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  df_gene <- data.frame(gene = genes, var = rv)
  rownames(df_gene) <- df_gene$gene
  v_var_ratio_c2t <- vapply(seq_len(iters_singlecell),
                            function(j) sum(df_gene[dic_ctrl_list[[j]], "var"]),
                            numeric(1))
  v_var_ratio_c2t <- v_var_ratio_c2t / sum(df_gene[scPagwas_topgenes, "var"])
  correct_pdf <- scPagwas:::correct_background(scPagwas.TRS.Score,
                                               mat_ctrl_raw_score,
                                               v_var_ratio_c2t)
  rownames(correct_pdf) <- colnames(Single_data)
  return(correct_pdf)
}
utils::assignInNamespace("Get_CorrectBg_p", patched_Get_CorrectBg_p, ns = "scPagwas")

patched_get_p_from_empi_null <- function(v_t, v_t_null) {
  cat("[PATCH 2] empirical p via findInterval ...\n")
  t0 <- Sys.time()
  sorted_null <- sort(v_t_null)
  v_pos <- findInterval(v_t, sorted_null)
  v_p <- (length(v_t_null) - v_pos + 1) / (length(v_t_null) + 1)
  cat(sprintf("[PATCH 2] done in %.1f sec\n",
              as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  return(v_p)
}
utils::assignInNamespace("get_p_from_empi_null", patched_get_p_from_empi_null,
                         ns = "scPagwas")
cat("[PATCH] both speed patches loaded.\n")

# === Paths ===
SHARED <- "/mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_scPagwas_shared_20260612"
OUT    <- "/mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_scPagwas_DIV_paperparams_20260615"
CHUNKS_DIR <- file.path(OUT, "sol1_chunks")
N_CHUNKS   <- 5
RDS_OUT    <- file.path(OUT, "MDD_scPagwas_DIV_paperparams.rds")
PMAT_FILE  <- file.path(OUT, "pmat_merge.RData")

# === Merge chunks ===
if (file.exists(PMAT_FILE)) {
  load(PMAT_FILE)
} else {
  cat("[", as.character(Sys.time()), "] Merging chunks ...\n")
  pmat_merge <- NULL
  for (i in seq_len(N_CHUNKS)) {
    f <- file.path(CHUNKS_DIR, sprintf("chunk_%d_results.RData", i))
    load(f)
    pmat_merge <- if (is.null(pmat_merge)) pmat else rbind(pmat_merge, pmat)
    rm(pmat)
  }
  save(pmat_merge, file = PMAT_FILE)
}
cat("pmat_merge dim:", paste(dim(pmat_merge), collapse=" x "), "\n")

# === Load shared Pagwas + Single_data ===
cat("[", as.character(Sys.time()), "] Loading shared pagwas_preprocessed.RData ...\n")
load(file.path(SHARED, "pagwas_preprocessed.RData"))

# === Merge_gPas ===
cat("[", as.character(Sys.time()), "] Merge_gPas ...\n")
scPagwas.gPAS.score <- scPagwas::Merge_gPas(Pagwas, pmat_merge)

# === PCC via Corr_Random (scGet_PCC OOMs on 310k cells) ===
cat("[", as.character(Sys.time()), "] Corr_Random PCC ...\n")
PCC <- scPagwas::Corr_Random(Pagwas$data_mat, scPagwas.gPAS.score,
                             seed = 1234, random = TRUE,
                             Nrandom = 5, Nselect = 200)

# === PAPER STYLE: top + down only (2 sets, NO up) with n=500 ===
n_topgenes <- 500
scPagwas_topgenes  <- names(sort(PCC, decreasing = TRUE))[1:n_topgenes]
scPagwas_downgenes <- names(sort(PCC, decreasing = FALSE))[1:n_topgenes]

cat("[", as.character(Sys.time()), "] AddModuleScore (TRS, downTRS) — 2 features ...\n")
Single_data <- Seurat::AddModuleScore(
  Single_data,
  assay = "RNA",
  features = list(scPagwas_topgenes, scPagwas_downgenes),
  name = c("scPagwas.TRS.Score", "scPagwas.downTRS.Score")
)
Single_data$scPagwas.TRS.Score     <- Single_data$scPagwas.TRS.Score1
Single_data$scPagwas.downTRS.Score <- Single_data$scPagwas.downTRS.Score2
Single_data$scPagwas.gPAS.score    <- scPagwas.gPAS.score[colnames(Single_data)]

# === Get_CorrectBg_p with PAPER PARAM n_topgenes=1000 ===
cat("[", as.character(Sys.time()), "] Get_CorrectBg_p (PATCHED, n_topgenes=1000 PAPER) ...\n")
correct_pdf <- scPagwas::Get_CorrectBg_p(
  Single_data = Single_data,
  scPagwas.TRS.Score = Single_data$scPagwas.TRS.Score,
  iters_singlecell = 100,
  n_topgenes = 1000,            # PAPER PARAM (was 500 in v2)
  scPagwas_topgenes = scPagwas_topgenes
)
# === BUG FIX: properly attribute pooled_p by cell name ===
pooled_p <- correct_pdf$pooled_p
names(pooled_p) <- rownames(correct_pdf)
Single_data$Random_Correct_BG_adjp <- pooled_p[colnames(Single_data)]
saveRDS(pooled_p, file = file.path(OUT, "pooled_p_paperparams.rds"))

# === Merge_celltype_p: BOTH major (anno) and subtype (final_anno) ===
cat("[", as.character(Sys.time()), "] Merge_celltype_p — major (anno, 8 categories) ...\n")
ct_major <- scPagwas::Merge_celltype_p(
  single_p = pooled_p,
  celltype = Single_data$anno[names(pooled_p)]
)
cat("\n=== Major cell-type p (paper-params) ===\n"); print(ct_major)
write.csv(ct_major,
          file = file.path(OUT, "cell_type_pvalues_paperparams.csv"),
          row.names = FALSE)

cat("[", as.character(Sys.time()), "] Merge_celltype_p — subtype (final_anno, 20 categories) ...\n")
ct_sub <- scPagwas::Merge_celltype_p(
  single_p = pooled_p,
  celltype = Single_data$final_anno[names(pooled_p)]
)
cat("\n=== Subtype p (paper-params) ===\n"); print(ct_sub)
write.csv(ct_sub,
          file = file.path(OUT, "cell_type_pvalues_paperparams_subtype.csv"),
          row.names = FALSE)
cell_type_p <- ct_major

cat("Saving final Seurat to", RDS_OUT, "\n")
saveRDS(Single_data, file = RDS_OUT)

cat("\n=== DONE (paper-params) ===\n")
