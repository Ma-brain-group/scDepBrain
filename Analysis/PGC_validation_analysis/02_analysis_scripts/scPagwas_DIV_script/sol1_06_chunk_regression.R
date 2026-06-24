# ============================================================
# HW paper-params: Link_pathway_blocks_gwas per chunk
# Uses chunks from SHARED folder + paper-params gwas_pagwas.RData
# ============================================================
suppressPackageStartupMessages({
  library(scPagwas)
})

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

chunk_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))
cat("Processing chunk_id =", chunk_id, "\n")

SHARED <- "/mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_scPagwas_shared_20260612"
OUT    <- "/mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_scPagwas_DIV_paperparams_20260615"
CHUNKS_DIR <- file.path(OUT, "sol1_chunks")
dir.create(CHUNKS_DIR, showWarnings = FALSE, recursive = TRUE)

CHUNK_RDATA <- file.path(SHARED, sprintf("chunk_%d_pagwas.RData", chunk_id))
GWAS_RDATA  <- file.path(OUT, "gwas_pagwas.RData")
OUT_FILE    <- file.path(CHUNKS_DIR, sprintf("chunk_%d_results.RData", chunk_id))
BACKINGPATH <- file.path(CHUNKS_DIR, sprintf("temp_%d", chunk_id))
dir.create(BACKINGPATH, showWarnings = FALSE, recursive = TRUE)

cat("[", as.character(Sys.time()), "] Loading paper-params GWAS Pagwas ...\n")
load(GWAS_RDATA)
Pagwas_gwas <- Pagwas
rm(Pagwas)

cat("[", as.character(Sys.time()), "] Loading chunk", chunk_id, "from SHARED ...\n")
load(CHUNK_RDATA)

Pagwas <- c(Pagwas_gwas, chunk)
rm(Pagwas_gwas, chunk); gc()

cat("[", as.character(Sys.time()), "] Pathway_annotation_input ...\n")
Pagwas <- scPagwas::Pathway_annotation_input(
  Pagwas = Pagwas,
  block_annotation = block_annotation
)

cat("[", as.character(Sys.time()), "] Link_pathway_blocks_gwas (chunk", chunk_id, ") ...\n")
Pagwas <- scPagwas::Link_pathway_blocks_gwas(
  Pagwas      = Pagwas,
  chrom_ld    = chrom_ld,
  singlecell  = TRUE,
  celltype    = FALSE,
  backingpath = BACKINGPATH,
  n.cores     = 1
)

pmat <- Pagwas$Pathway_sclm_results
cat("Saving Pathway_sclm_results dim:", paste(dim(pmat), collapse=" x "), "\n")
save(pmat, file = OUT_FILE, compress = TRUE)

unlink(BACKINGPATH, recursive = TRUE)
cat("[", as.character(Sys.time()), "] Chunk", chunk_id, "DONE.\n")
