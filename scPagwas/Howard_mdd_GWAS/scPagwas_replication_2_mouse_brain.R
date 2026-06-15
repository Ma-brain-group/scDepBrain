############ scPagwas运行标准代码 ##################
library(scPagwas)
library(Seurat)
library(dplyr)
library(data.table)
library(Matrix)
library(stringr)
library(parallel)
library(irlba)
library(glmnet)
library(GenomicRanges)
library(utils)
library(ggplot2)
library(ggthemes)
library(ggpubr)

############# 设置路径
setwd("/public/ojsys/eye/sujianzhong/chencheng/Singlecell/MDD/Mouse_data/")

############# 输入数据，通路数据，基因注释信息，LD文件都可以用自带的（这里要注意GWAS数据和单细胞的参考基因组版本）
data(Genes_by_pathway_kegg)
data(chrom_ld)
load("/public/ojsys/eye/sujianzhong/chencheng/Singlecell/MDD/scPagwas/block_annotation_hg37.RData")

############# 这里采用分步计算的方法
############# 单细胞数据读入
Pagwas <- list()
Single_data <- readRDS("/public/ojsys/eye/sujianzhong/chencheng/Singlecell/MDD/Mouse_data/mouse_merge1.rds")

############# 指定细胞类型计算SVD矩阵
Idents(Single_data)<-Single_data$brain_region
Pagwas <- Single_data_input(
  Pagwas = Pagwas,
  assay = "RNA",
  Single_data = Single_data,
  Pathway_list = Genes_by_pathway_kegg
)
Single_data <- Single_data[, colnames(Pagwas$data_mat)]


############## 计算通路的PCA score
Pagwas <- Pathway_pcascore_run(
  Pagwas = Pagwas,
  Pathway_list = Genes_by_pathway_kegg
)


############### GWAS summary数据输入(SNP数目大的数据，建议进行SNP过滤)
gwas_data <- bigreadr::fread2("/public/ojsys/eye/sujianzhong/chencheng/Singlecell/MDD/scPagwas/scPagwas_input_LD_prune/ieu_b_102.pagwas.txt")
Pagwas <- GWAS_summary_input(
  Pagwas = Pagwas,
  gwas_data = gwas_data,
  maf_filter = 0.1
)

############### 将SNP注释到基因上
Pagwas$snp_gene_df <- SnpToGene(
  gwas_data = Pagwas$gwas_data,
  block_annotation = block_annotation_hg37,
  marg = 10000
)


################ Pathway-SNP的注释
Pagwas <- Pathway_annotation_input(
  Pagwas = Pagwas,
  block_annotation = block_annotation_hg37
)


################ Link the pathway blocks to pca score
Pagwas <- Link_pathway_blocks_gwas(
  Pagwas = Pagwas,
  chrom_ld = chrom_ld,
  singlecell = T,
  celltype = T,
  backingpath="./temp")
#save(Pagwas,file=paste0(sub_dir,"_Pagwas.Rdata")) #这一步运行时间很长，建议保存一下

################ Perform regression for celltypes
Pagwas$lm_results <- Pagwas_perform_regression(Pathway_ld_gwas_data = Pagwas$Pathway_ld_gwas_data)
Pagwas <- Boot_evaluate(Pagwas, bootstrap_iters = 200, part = 0.5)
Pagwas$Pathway_ld_gwas_data <- NULL ##remove the Pathway_ld_gwas_data, it takes a lot of memory.

################ Construct the scPagwas score
Pagwas <- scPagwas_perform_score(
  Pagwas = Pagwas,
  remove_outlier = TRUE
)

############### Get the pearson correlation coefficients for gene(PCC)
Pagwas$PCC <- scPagwas::scGet_PCC(scPagwas.gPAS.score=Pagwas$scPagwas.gPAS.score,
                                  data_mat=Pagwas$data_mat)

############## Calculate the TRS score for top genes
############## Calculate the TRS score for top genes
############## 先取出正负top500的基因
n_topgenes=500
scPagwas_topgenes <- rownames(Pagwas$PCC)[order(Pagwas$PCC, decreasing = T)[1:n_topgenes]]
scPagwas_downgenes <- rownames(Pagwas$PCC)[order(Pagwas$PCC, decreasing =F)[1:n_topgenes]]
Single_data <- Seurat::AddModuleScore(Single_data, assay = "RNA", list(scPagwas_topgenes,scPagwas_downgenes), name = c("scPagwas.TRS.Score","scPagwas.downTRS.Score"))

############## Calculate the p-values for scPagwas.TRS.Score of single cells after background correction.
correct_pdf<-Get_CorrectBg_p(Single_data=Single_data,
                             scPagwas.TRS.Score=Single_data$scPagwas.TRS.Score1,
                             iters_singlecell=100,    ## iters_singlecell=100计算出的结果比较准确
                             n_topgenes=1000,
                             scPagwas_topgenes=scPagwas_topgenes,
                             assay="RNA")
Pagwas$Random_Correct_BG_pdf <- correct_pdf

############## Merge the p-values of cells belonging to the same cell type into a single p-value for each cell type.
############## 这一步计算出细胞类型的结果
Pagwas$Merged_celltype_pvalue<-Merge_celltype_p(single_p=correct_pdf$pooled_p,celltype=Pagwas$Celltype_anno$annotation)
celltype_result<-Pagwas$Merged_celltype_pvalue
write.csv(celltype_result,file = "celltype_result.csv")

############## 保存Pagwas的结果（所有结果都保存在Pagwas里面）
#save(Pagwas,file="Pagwas_all_result.Rdata")

############## Integrate and output the results of single-cell analysis.
############## 整合单细胞的结果
a <- data.frame(
  scPagwas.TRS.Score = Single_data$scPagwas.TRS.Score1,
  scPagwas.downTRS.Score = Single_data$scPagwas.downTRS.Score2,
  scPagwas.gPAS.score = Pagwas$scPagwas.gPAS.score,
  Random_Correct_BG_p = correct_pdf$pooled_p,
  Random_Correct_BG_adjp = correct_pdf$adj_p,
  Random_Correct_BG_z = correct_pdf$pooled_z)

############## 保存TRS score结果
utils::write.csv(a,file="singlecell_scPagwas_score_pvalue_Result.csv",quote=F)
