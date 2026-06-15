# Step 1 TRS plot.R
```R
#######################scPagwas TRS plot#############################
#######################单细胞数据读入
library(Seurat)
library(ggplot2)
library(cowplot)
library(ggpubr)
library(Nebulosa)
library(BiocFileCache)
library(paletteer)
library(scCustomize)
setwd("D:/Project/SingleCell_MDD/GWAS_combine_singlecell_analysis/1_GWAS_analysis_for_broad_celltype/scPagwas/scPagwas_singlecell_result_block_annotation_hg37_1/ieu_b_102/")
single_data<-readRDS("D:/Project/SingleCell_MDD/SingleCell_analysis/1_Singlecell_data/simple.rds")

######################按照single_data的列名，match一下metadata的行名，细胞顺序一致
#colnames<-colnames(single_data)
#metadata<-single_data@meta.data
#metadata1<-metadata[match(colnames,row.names(metadata)),]
#single_data@meta.data<-metadata1

##保存数据
#saveRDS(single_data,file = "C:/Users/Lenovo/Desktop/simple.rds")

#######################导入TRS score数据
#load("C:/Users/Lenovo/Desktop/TRS_plot/scPagwas_singlecell_result_block_annotation_hg37_1/2013_NG_15_risk_loci_MDD/2013_NG_15_risk_loci_MDD_result.Rdata")
#######################将TRS数据加到单细胞数据中
#single_data$scPagwas.TRS.Score<-a$scPagwas.TRS.Score
metadata<-single_data@meta.data

####################### boxplot #############################
metadata_anno<-metadata[,c("anno","scPagwas.TRS.Score")]
metadata_anno1<-aggregate(metadata_anno$scPagwas.TRS.Score,by=list(metadata_anno$anno),FUN=median)
metadata_anno1$color<-c("#d5231d","#3777ac","#4ea64a","#8e4c99","#e88f18","#e47faf","#b698c5","#a05528")
metadata_anno1<-metadata_anno1[order(metadata_anno1$x,decreasing = T),]
metadata_anno$anno<-factor(metadata_anno$anno,levels = metadata_anno1$Group.1)

pdf("TRS_boxplot.pdf",height = 5,width = 6)
ggplot(metadata_anno, aes(x=anno, y=scPagwas.TRS.Score)) + 
  stat_boxplot(geom = "errorbar",width=0.05, size=0.5,position=position_dodge(0.6),color= "black")+
  theme_bw()+
  geom_boxplot(position = position_dodge(0.6),
               size = 0.5,
               width = 0.8,
               fill = metadata_anno1$color,
               color = "black",
               outlier.color = "black",
               outlier.fill = "black",
               outlier.shape = 19,
               outlier.size = 1,
               outlier.stroke = 0.5,
               outlier.alpha = 45,
               notch = F,
               notchwidth = 0.5)+
  xlab("Cell type")+
  ylab("scPagwas TRS Score")+
  theme(axis.text.x = element_text(size = 10,color = 'black',vjust = 1,hjust = 1,angle = 35), # 调整x轴坐标文字
        axis.text.y = element_text(size = 10,color = 'black'), # 调整y轴坐标文字
        legend.text = element_text(size = 10,color = 'black'), # 调整legend字体大小
        legend.title = element_text(size = 10,color = 'black'), # 调整legend title大小
        panel.grid.major =  element_blank(),
        panel.grid.minor = element_blank()
        )
dev.off()

##################### Featureplot ##############################
####取出细胞的坐标值
embedding<-single_data@reductions$umap@cell.embeddings
embedding<-cbind(embedding,metadata_anno$scPagwas.TRS.Score)
embedding<-as.data.frame(embedding)
colnames(embedding)[3]<-"scPagwas.TRS.Score"

pdf("TRS_featureplot.pdf",height = 6)
ggplot(embedding,aes(x=UMAP_1,y=UMAP_2,color=scPagwas.TRS.Score))+
  geom_point(size=0.2,alpha=0.2)+
  theme_bw()+
  scale_colour_gradient2(low="blue",mid="white",high="red",midpoint=0.2)+
  labs(color="scPagwas.TRS.Score")+
  theme_cowplot()
dev.off()


################## density plot
pdf("scPagwas.TRS.Score_density_plot.pdf",height = 4.5,width = 5)
Plot_Density_Custom(seurat_object =single_data, features = "scPagwas.TRS.Score",reduction = "umap",
                    custom_palette = c("#B0CFE4","#FACABC","#E77A77","#DC0000FF"))+
  theme(plot.title = element_text(hjust = 0.5))
dev.off()


################## 导入细胞类型的结果 ###########################
load("C:/Users/Lenovo/Desktop/TRS_plot/scPagwas_singlecell_result_block_annotation_hg37_1/2013_NG_15_risk_loci_MDD/2013_NG_15_risk_loci_MDD_celltype_result.Rdata")
write.csv(temp,file = "C:/Users/Lenovo/Desktop/TRS_plot/scPagwas_singlecell_result_block_annotation_hg37_1/2013_NG_15_risk_loci_MDD/2013_NG_15_risk_loci_MDD_celltype_result.csv",row.names = F,quote = F)



```


# Step 2 scPagwas cell type inference
```R
##scPagwas运行代码
##单细胞的结果
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

##设置路径
main_dir<-"/public/ojsys/eye/sujianzhong/chencheng/Singlecell/MDD/scPagwas/scPagwas_singlecell_result_block_annotation_hg37_1/"
sub_dir<-"uk_self_depression"

if (file.exists(sub_dir)){
  setwd(file.path(main_dir, sub_dir))
} else {
  dir.create(file.path(main_dir, sub_dir))
  setwd(file.path(main_dir, sub_dir))
}

##输入数据，通路数据，基因注释信息，LD文件都可以用自带的
data(Genes_by_pathway_kegg)
data(chrom_ld)
load("/public/ojsys/eye/sujianzhong/chencheng/Singlecell/MDD/scPagwas/scPagwas_singlecell_result_block_annotation_hg37/block_annotation_hg37.RData")

##这里采用分步计算的方法

##3.1Single data input
Pagwas <- list()
Single_data <- readRDS("/public/ojsys/eye/sujianzhong/chencheng/Singlecell/MDD/scPagwas/simple.rds")
##指定anno计算SVD矩阵
Idents(Single_data)<-Single_data$anno
Pagwas <- Single_data_input(
      Pagwas = Pagwas,
      assay = "RNA",
      Single_data = Single_data,
      Pathway_list = Genes_by_pathway_kegg
    )
Single_data <- Single_data[, colnames(Pagwas$data_mat)]


##3.2Run pathway pca score
Pagwas <- Pathway_pcascore_run(
        Pagwas = Pagwas,
        Pathway_list = Genes_by_pathway_kegg
      )
      

##3.3GWAS summary data input(这里的GWAS数据，小数据我没有LD过滤，大数据我进行LD过滤了)
gwas_data <- bigreadr::fread2("/public/ojsys/eye/sujianzhong/chencheng/Singlecell/MDD/scPagwas/scPagwas_input_LD_prune/uk_self_depression.pagwas.txt")
Pagwas <- GWAS_summary_input(
    Pagwas = Pagwas,
    gwas_data = gwas_data,
    maf_filter = 0.1
  )


##3.4Mapping Snps to Genes
Pagwas$snp_gene_df <- SnpToGene(
        gwas_data = Pagwas$gwas_data,
        block_annotation = block_annotation_hg37,
        marg = 10000
      )


##3.5Pathway-SNP annotation
Pagwas <- Pathway_annotation_input(
      Pagwas = Pagwas,
      block_annotation = block_annotation_hg37
    )
    
    
##3.6Link the pathway blocks to pca score
Pagwas <- Link_pathway_blocks_gwas(
      Pagwas = Pagwas,
      chrom_ld = chrom_ld,
      singlecell = T,
      celltype = T,
      backingpath="./temp")
##这步运行很长时间，保存一下结果
###save(Pagwas,file=paste0(sub_dir,"_Pagwas.Rdata"))
      

##3.7Perform regression for celltypes
Pagwas$lm_results <- Pagwas_perform_regression(Pathway_ld_gwas_data = Pagwas$Pathway_ld_gwas_data)
Pagwas <- Boot_evaluate(Pagwas, bootstrap_iters = 200, part = 0.5)
##remove the Pathway_ld_gwas_data, it takes a lot of memory.
Pagwas$Pathway_ld_gwas_data <- NULL


##3.8Construct the scPagwas score
Pagwas <- scPagwas_perform_score(
      Pagwas = Pagwas,
      remove_outlier = TRUE
    )


##3.9Get the pearson correlation coefficients for gene(PCC)
Pagwas$PCC <- scPagwas::scGet_PCC(scPagwas.gPAS.score=Pagwas$scPagwas.gPAS.score,
                                    data_mat=Pagwas$data_mat)
                                    

##3.10 Calculate the TRS score for top genes
##Calculate the TRS score for top genes
##先取出正负top500的基因
n_topgenes=500
scPagwas_topgenes <- rownames(Pagwas$PCC)[order(Pagwas$PCC, decreasing = T)[1:n_topgenes]]
scPagwas_downgenes <- rownames(Pagwas$PCC)[order(Pagwas$PCC, decreasing =F)[1:n_topgenes]]

Single_data <- Seurat::AddModuleScore(Single_data, assay = "RNA", list(scPagwas_topgenes,scPagwas_downgenes), name = c("scPagwas.TRS.Score","scPagwas.downTRS.Score"))

#Calculate the p-values for scPagwas.TRS.Score of single cells after background correction.
correct_pdf<-Get_CorrectBg_p(Single_data=Single_data,
                             scPagwas.TRS.Score=Single_data$scPagwas.TRS.Score1,
                             iters_singlecell=100,    ## iters_singlecell=100计算出的结果比较准确
                             n_topgenes=1000,
                             scPagwas_topgenes=scPagwas_topgenes,
                             assay="RNA")
Pagwas$Random_Correct_BG_pdf <- correct_pdf

#Merge the p-values of cells belonging to the same cell type into a single p-value for each cell type.
Pagwas$Merged_celltype_pvalue<-Merge_celltype_p(single_p=correct_pdf$pooled_p,celltype=Pagwas$Celltype_anno$annotation)
temp<-Pagwas$Merged_celltype_pvalue
save(temp,file=paste0(sub_dir,"_celltype_result.Rdata"))


##保存Pagwas的结果（所有结果都保存在Pagwas里面）
save(Pagwas,file=paste0(sub_dir,"_Pagwas_all_result.Rdata"))

#Integrate and output the results of single-cell analysis.
a <- data.frame(
     scPagwas.TRS.Score = Single_data$scPagwas.TRS.Score1,
     scPagwas.downTRS.Score = Single_data$scPagwas.downTRS.Score2,
     scPagwas.gPAS.score = Pagwas$scPagwas.gPAS.score,
     Random_Correct_BG_p = correct_pdf$pooled_p,
     Random_Correct_BG_adjp = correct_pdf$adj_p,
     Random_Correct_BG_z = correct_pdf$pooled_z)

##保存TRS score结果
save(a,file=paste0(sub_dir,"_result.Rdata"))
utils::write.csv(a,file=paste0(sub_dir,"_singlecell_scPagwas_score_pvalue_Result.csv"),quote=F)

```
