# =================================================================
# Convert MDD single-cell rds → h5ad for scDRS
# Keep: raw counts + metadata (anno, final_anno, batch, sex, n_genes)
# =================================================================
suppressPackageStartupMessages({
  library(Seurat); library(SeuratDisk); library(Matrix)
})
utils::assignInNamespace(".Deprecate", function(...) invisible(NULL),
                         ns = "SeuratObject")

RDS <- "/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/integrated_heathy_brain_cells_short_read/MDD_singlecell_data_reannotation.rds"
OUT_DIR <- "/mnt/isilon/gandal_lab/mayl/07_scDepBrain/results_scDRS_20260618/data"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

cat("[", as.character(Sys.time()), "] Loading 24GB rds ...\n")
t0 <- Sys.time()
obj <- readRDS(RDS)
cat(sprintf("Loaded in %.1f min\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))

DefaultAssay(obj) <- "RNA"
if (inherits(obj[["RNA"]], "Assay5")) {
  obj <- JoinLayers(obj, assay = "RNA")
  obj[["RNA"]] <- as(obj[["RNA"]], "Assay")
}
# Drop SCT assay if exists (saves memory + we don't need it)
if ("SCT" %in% Assays(obj)) { obj[["SCT"]] <- NULL; gc() }
# Drop integrated assay if exists
for (a in setdiff(Assays(obj), "RNA")) { obj[[a]] <- NULL }
gc()

# Slim metadata to essentials
md <- obj@meta.data
keep_cols <- intersect(colnames(md),
  c("anno","final_anno","Disease","sample","donor","batch","sex",
    "Sex","n_genes","nFeature_RNA","nCount_RNA","percent.mt"))
cat("Keeping metadata cols:", paste(keep_cols, collapse=", "), "\n")
obj@meta.data <- md[, keep_cols, drop = FALSE]
gc()

cat("Cells:", ncol(obj), "  Genes:", nrow(obj), "\n")
cat("Metadata cols:\n"); print(colnames(obj@meta.data))

# SaveH5Seurat → convert to h5ad
H5S <- file.path(OUT_DIR, "MDD_scDRS.h5Seurat")
H5A <- file.path(OUT_DIR, "MDD_scDRS.h5ad")
file.remove(c(H5S, H5A))  # clean up if any

cat("[", as.character(Sys.time()), "] SaveH5Seurat -> ", H5S, "\n")
SaveH5Seurat(obj, filename = H5S, overwrite = TRUE)
cat("[", as.character(Sys.time()), "] Convert -> ", H5A, "\n")
Convert(H5S, dest = "h5ad", overwrite = TRUE)

file.remove(H5S)  # cleanup intermediate
cat("[", as.character(Sys.time()), "] DONE.\n")
system(paste("ls -lh", OUT_DIR))
