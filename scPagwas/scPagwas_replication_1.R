library(Seurat)
library(Matrix)
library(scPagwas)
library(arrow)
load("/public/ojsys/eye/sujianzhong/chencheng/Singlecell/MDD/scPagwas/block_annotation_hg37.RData")
load("protein_gene.Rdata")


##>> Because the cell number in this dataset was very large, we randomly sampled 500,000 cells for scPagwas analysis. A random seed was set in the Python script to ensure the reproducibility of the sampling results.

#### 设置路径
setwd("/public/ojsys/eye/sujianzhong/chencheng/Singlecell/MDD/300M_singlecell_data/")
i="merged_3M_singlecell_data_random_1M"

##### 读入count matrix
data_df <- arrow::read_feather(paste0(i,'_matrix.feather'))
my_r_matrix<-data.matrix(data_df)

##### 读入基因信息和metadata
obs<-read.csv(paste0(i,"_obs.csv"),header=T)
var<-read.csv(paste0(i,"_var_names.csv"))
my_r_matrix <- my_r_matrix[, -ncol(my_r_matrix)]
rownames(obs)<-colnames(my_r_matrix)

##### get the gene in block_annotation_hg37
rownames(my_r_matrix)<-var[,"Gene"]
my_r_matrix<-my_r_matrix[row.names(my_r_matrix)%in%protein_gene$gene_name,]
dim(my_r_matrix)
#my_r_matrix <- Matrix(my_r_matrix, sparse = TRUE)

##### 构建seurat对象
Single_data<-CreateSeuratObject(counts = my_r_matrix,assay = "RNA",meta.data = obs,min.cells=1)
Single_data
Idents(Single_data)<-Single_data$supercluster_term
#saveRDS(Single_data,file="merged_3M_singlecell_data_random_1M.rds")

##### remove the matrix
rm(my_r_matrix)
rm(data_df)
gc()

##### load gene pathway, LD matrix, gene annotation file
data(Genes_by_pathway_kegg)
data(chrom_ld)

##### Single data input
Pagwas <- list()
Pagwas <- Single_data_input(
  Pagwas = Pagwas,
  assay = "RNA",
  Single_data = Single_data,
  Pathway_list = Genes_by_pathway_kegg
)
Pagwas$VariableFeatures<-rownames(Pagwas$data_mat)
Single_data <- Single_data[, colnames(Pagwas$data_mat)]


##### Run pathway pca score
Pagwas <- Pathway_pcascore_run(
  Pagwas = Pagwas,
  Pathway_list = Genes_by_pathway_kegg
)


##### GWAS summary data input(LD filter)
gwas_data <- bigreadr::fread2("/public/ojsys/eye/sujianzhong/chencheng/Singlecell/MDD/scPagwas/scPagwas_input_LD_prune/ieu_b_102.pagwas.txt")
Pagwas <- GWAS_summary_input(
  Pagwas = Pagwas,
  gwas_data = gwas_data,
  maf_filter = 0.1
)


##### Mapping Snps to Genes
Pagwas$snp_gene_df <- SnpToGene(
  gwas_data = Pagwas$gwas_data,
  block_annotation = block_annotation_hg37,
  marg = 10000
)


##### Pathway-SNP annotation
Pagwas <- Pathway_annotation_input(
  Pagwas = Pagwas,
  block_annotation = block_annotation_hg37
)


##### Link the pathway blocks to pca score
Pagwas <- Link_pathway_blocks_gwas(
  Pagwas = Pagwas,
  chrom_ld = chrom_ld,
  singlecell = T,
  celltype = T,
  backingpath="./temp")
##save(Pagwas,file=paste0(sub_dir,"_Pagwas.Rdata"))


###### Perform regression for celltypes
Pagwas$lm_results <- Pagwas_perform_regression(Pathway_ld_gwas_data = Pagwas$Pathway_ld_gwas_data)
Pagwas <- Boot_evaluate(Pagwas, bootstrap_iters = 200, part = 0.5)
##remove the Pathway_ld_gwas_data, it takes a lot of memory.
Pagwas$Pathway_ld_gwas_data <- NULL


##### Construct the scPagwas score
Pagwas <- scPagwas_perform_score(
  Pagwas = Pagwas,
  remove_outlier = TRUE
)

##### Get the pearson correlation coefficients for gene(PCC)
Pagwas$PCC <- scPagwas::scGet_PCC(scPagwas.gPAS.score=Pagwas$scPagwas.gPAS.score,
                                  data_mat=Pagwas$data_mat)


##### Calculate the TRS score for top genes
##Calculate the TRS score for top genes
n_topgenes=500
scPagwas_topgenes <- rownames(Pagwas$PCC)[order(Pagwas$PCC, decreasing = T)[1:n_topgenes]]
scPagwas_downgenes <- rownames(Pagwas$PCC)[order(Pagwas$PCC, decreasing =F)[1:n_topgenes]]
Single_data <- Seurat::AddModuleScore(Single_data, assay = "RNA", list(scPagwas_topgenes,scPagwas_downgenes), name = c("scPagwas.TRS.Score","scPagwas.downTRS.Score"))

##################### Calculate the p-values for scPagwas.TRS.Score of single cells after background correction.
correct_pdf<-Get_CorrectBg_p(Single_data=Single_data,
                             scPagwas.TRS.Score=Single_data$scPagwas.TRS.Score1,
                             iters_singlecell=100,    ## iters_singlecell=100
                             n_topgenes=1000,
                             scPagwas_topgenes=scPagwas_topgenes,
                             assay="RNA")
Pagwas$Random_Correct_BG_pdf <- correct_pdf

##################### Merge the p-values of cells belonging to the same cell type into a single p-value for each cell type.
Pagwas$Merged_celltype_pvalue<-Merge_celltype_p(single_p=correct_pdf$pooled_p,celltype=Pagwas$Celltype_anno$annotation)
#### save the celltype result
write.csv(Pagwas$Merged_celltype_pvalue,file = "celltype_result.csv")

#### Integrate and output the results of single-cell analysis.
a <- data.frame(
  scPagwas.TRS.Score = Single_data$scPagwas.TRS.Score1,
  scPagwas.downTRS.Score = Single_data$scPagwas.downTRS.Score2,
  scPagwas.gPAS.score = Pagwas$scPagwas.gPAS.score,
  Random_Correct_BG_p = correct_pdf$pooled_p,
  Random_Correct_BG_adjp = correct_pdf$adj_p,
  Random_Correct_BG_z = correct_pdf$pooled_z)

#### save the singlecell result
write.csv(a,file="singlecell_scPagwas_score_pvalue_Result.csv",quote=F)
