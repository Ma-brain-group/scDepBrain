# ============================================================
# Solution 1 — Step 1.4: Link_pathway_blocks_gwas per chunk
#
# Reads chunk_id from SLURM_ARRAY_TASK_ID, processes ONE chunk
# of single-cell data for THIS GWAS, saves Pathway_sclm_results.
#
# Memory: ~ 200-400 GB per chunk (LD matrices dominate)
# Runtime: ~2-12 h per chunk depending on pathway sizes
# ============================================================
suppressPackageStartupMessages({
  library(scPagwas)
})

# BLAS / parallel patches
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

# --- Determine chunk_id from SLURM array task ---
chunk_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))
cat("Processing chunk_id =", chunk_id, "\n")

# --- Paths ---
SHARED <- "/mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_scPagwas_shared_20260612"
OUT    <- "/mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_scPagwas_MDD2025_EUR_20260610"
CHUNKS_DIR <- file.path(OUT, "sol1_chunks")
dir.create(CHUNKS_DIR, showWarnings = FALSE, recursive = TRUE)

CHUNK_RDATA <- file.path(SHARED, sprintf("chunk_%d_pagwas.RData", chunk_id))
GWAS_RDATA  <- file.path(OUT, "gwas_pagwas.RData")
OUT_FILE    <- file.path(CHUNKS_DIR, sprintf("chunk_%d_results.RData", chunk_id))
BACKINGPATH <- file.path(CHUNKS_DIR, sprintf("temp_%d", chunk_id))
dir.create(BACKINGPATH, showWarnings = FALSE, recursive = TRUE)

# --- Load GWAS Pagwas first ---
cat("[", as.character(Sys.time()), "] Loading GWAS Pagwas ...\n")
load(GWAS_RDATA)        # brings 'Pagwas' (gwas_data + snp_gene_df + rawPathway_list)
Pagwas_gwas <- Pagwas
rm(Pagwas)

# --- Load chunk Pagwas ---
cat("[", as.character(Sys.time()), "] Loading chunk", chunk_id, "...\n")
load(CHUNK_RDATA)       # brings 'chunk' (Single_data_input output for these cells)

# --- Merge into one Pagwas list ---
Pagwas <- c(Pagwas_gwas, chunk)   # list concat preserves all named elements
rm(Pagwas_gwas, chunk); gc()
cat("Merged Pagwas. data_mat:", paste(dim(Pagwas$data_mat), collapse=" x "),
    " pca_scCell_mat:", paste(dim(Pagwas$pca_scCell_mat), collapse=" x "), "\n")

# --- Pathway_annotation_input (matches scPagwas internal usage) ---
cat("[", as.character(Sys.time()), "] Pathway_annotation_input ...\n")
Pagwas <- scPagwas::Pathway_annotation_input(
  Pagwas = Pagwas,
  block_annotation = block_annotation
)

# --- The heavy step: Link_pathway_blocks_gwas (singlecell=TRUE, celltype=FALSE) ---
cat("[", as.character(Sys.time()), "] Link_pathway_blocks_gwas (chunk", chunk_id, ") ...\n")
Pagwas <- scPagwas::Link_pathway_blocks_gwas(
  Pagwas      = Pagwas,
  chrom_ld    = chrom_ld,
  singlecell  = TRUE,
  celltype    = FALSE,
  backingpath = BACKINGPATH,
  n.cores     = 1
)

# --- Save Pathway_sclm_results (per-cell regression results) ---
pmat <- Pagwas$Pathway_sclm_results
cat("Saving Pathway_sclm_results dim:", paste(dim(pmat), collapse=" x "),
    "to", OUT_FILE, "\n")
save(pmat, file = OUT_FILE, compress = TRUE)

# Cleanup temp .bk files
unlink(BACKINGPATH, recursive = TRUE)

cat("[", as.character(Sys.time()), "] Chunk", chunk_id, "DONE.\n")
