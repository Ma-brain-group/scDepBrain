##################### According to the celltype percent in each brain region, select cells to run scPagwas ######################
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
set.seed(123)

##################### set the workspace
setwd("/public/ojsys/eye/sujianzhong/chencheng/Singlecell/MDD/Brain_region_scPagaws_random/In_all_cells/")

##################### import the 30w singlecell data
Single_data<-readRDS("MDD_singlecell_data_reannotation_simple.rds")
## min cells:2469
min(table(Single_data$tissue.region))
Idents(Single_data)<-Single_data$anno

##################### according to the celltype percent in each brain region, select 1000 cells randomly to run scPagwas
### get the brain region
regions <- unique(Single_data$tissue.region)
### construct an empty list 
cells_selected <- list()
for (region in regions) {
  ## get the cells in each brain region
  cells_in_region <- subset(Single_data, subset = tissue.region == region)
  ## calculate the celltype percent in each brain region
  cell_type_proportions <- prop.table(table(cells_in_region$anno))
  ## calculate the cell numbers
  n_per_type <- round(cell_type_proportions * 1000)
  ## according to the cell numbers, select cells randomly
  sampled_cells <- lapply(names(n_per_type), function(anno) {
    cells <- WhichCells(cells_in_region, idents = anno)
    sample(cells, n_per_type[anno])
  })
  cells_selected[[region]] <- unlist(sampled_cells)
}

#################### get the new singlecell data
selected_cells <- unlist(cells_selected)
Single_data <- subset(Single_data, cells = selected_cells)

#################### load gene pathway, LD matrix, gene annotation file
data(Genes_by_pathway_kegg)
data(chrom_ld)
load("/public/ojsys/eye/sujianzhong/chencheng/Singlecell/MDD/scPagwas/block_annotation_hg37.RData")


#################### Single data input
Pagwas <- list()
Pagwas <- Single_data_input(
  Pagwas = Pagwas,
  assay = "RNA",
  Single_data = Single_data,
  Pathway_list = Genes_by_pathway_kegg
)
Single_data <- Single_data[, colnames(Pagwas$data_mat)]


#################### Run pathway pca score
Pagwas <- Pathway_pcascore_run(
  Pagwas = Pagwas,
  Pathway_list = Genes_by_pathway_kegg
)


##################### GWAS summary data input(LD filter)
gwas_data <- bigreadr::fread2("/public/ojsys/eye/sujianzhong/chencheng/Singlecell/MDD/scPagwas/scPagwas_input_LD_prune/ieu_b_102.pagwas.txt")
Pagwas <- GWAS_summary_input(
  Pagwas = Pagwas,
  gwas_data = gwas_data,
  maf_filter = 0.1
)


##################### Mapping Snps to Genes
Pagwas$snp_gene_df <- SnpToGene(
  gwas_data = Pagwas$gwas_data,
  block_annotation = block_annotation_hg37,
  marg = 10000
)


##################### Pathway-SNP annotation
Pagwas <- Pathway_annotation_input(
  Pagwas = Pagwas,
  block_annotation = block_annotation_hg37
)


##################### Link the pathway blocks to pca score
Pagwas <- Link_pathway_blocks_gwas(
  Pagwas = Pagwas,
  chrom_ld = chrom_ld,
  singlecell = T,
  celltype = T,
  backingpath="./temp")
##save(Pagwas,file=paste0(sub_dir,"_Pagwas.Rdata"))


##################### Perform regression for celltypes
Pagwas$lm_results <- Pagwas_perform_regression(Pathway_ld_gwas_data = Pagwas$Pathway_ld_gwas_data)
Pagwas <- Boot_evaluate(Pagwas, bootstrap_iters = 200, part = 0.5)
##remove the Pathway_ld_gwas_data, it takes a lot of memory.
Pagwas$Pathway_ld_gwas_data <- NULL


##################### Construct the scPagwas score
Pagwas <- scPagwas_perform_score(
  Pagwas = Pagwas,
  remove_outlier = TRUE
)

##################### Get the pearson correlation coefficients for gene(PCC)
Pagwas$PCC <- scPagwas::scGet_PCC(scPagwas.gPAS.score=Pagwas$scPagwas.gPAS.score,
                                  data_mat=Pagwas$data_mat)


##################### Calculate the TRS score for top genes
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


###################### plot the TRS result of brain region
setwd("D:/Project/SingleCell_MDD/SingleCell_analysis/Brain_region_scPagaws_random/In_all_cells/")
single_data<-readRDS("D:/Project/SingleCell_MDD/SingleCell_analysis/annotation/MDD_singlecell_data_reannotation_simple.rds")
scPagwas.TRS.Score<-read.csv("singlecell_scPagwas_score_pvalue_Result.csv",row.names = 1)

###################### match the brain region
brain_region<-single_data$tissue.region
brain_region<-brain_region[match(row.names(scPagwas.TRS.Score),names(brain_region))]
scPagwas.TRS.Score$brain_region<-brain_region
scPagwas.TRS.Score1<-aggregate(scPagwas.TRS.Score$scPagwas.TRS.Score,by=list(scPagwas.TRS.Score$brain_region),median)
scPagwas.TRS.Score1<-scPagwas.TRS.Score1[order(scPagwas.TRS.Score1$x,decreasing = T),]
scPagwas.TRS.Score$brain_region<-factor(scPagwas.TRS.Score$brain_region,levels = scPagwas.TRS.Score1$Group.1)

pdf("Random_scPagwas.TRS.Score_in_all_cells_of_brain_region.pdf",height = 5,width = 6)
ggplot(scPagwas.TRS.Score, aes(x=brain_region, y=scPagwas.TRS.Score)) + 
  stat_boxplot(geom = "errorbar",width=0.05, size=0.5,position=position_dodge(0.6),color= "black")+
  theme_bw()+
  geom_boxplot(position = position_dodge(0.6),
               size = 0.5,
               width = 0.7,
               fill = c("#CAA2F4","#96873B","#B37557","#ED4437","#89C8E8","#6A3D9A","#B3446C","#EBD57C","#B49D99",
                        "#E68FAC","#E1884A","#FC9A9A","#1F78B4","#8ACC72"),
               color = "black",
               outlier.color = "black",
               outlier.fill = "black",
               outlier.shape = 19,
               outlier.size = 1,
               outlier.stroke = 0.5,
               outlier.alpha = 45,
               notch = F,
               notchwidth = 0.5)+
  xlab("Brain region")+
  ylab("Random scPagwas.TRS.Score")+
  theme(axis.text.x = element_text(size = 10,color = 'black',vjust = 1,hjust = 1,angle = 35),
        axis.text.y = element_text(size = 10,color = 'black'), 
        legend.text = element_text(size = 10,color = 'black'), 
        legend.title = element_text(size = 10,color = 'black'), 
        panel.grid.major =  element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5))+
  ggtitle("Randomly select")
dev.off()
