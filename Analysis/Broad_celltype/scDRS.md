
# Step 1 scDRS.score_plot.R
```r
############## scDRS score plot #####################
library(data.table)
library(ggplot2)
library(cowplot)
library(Nebulosa)
library(BiocFileCache)
library(paletteer)
library(scCustomize)

setwd("D:/Project/SingleCell_MDD/GWAS_combine_singlecell_analysis/1_GWAS_analysis_for_broad_celltype/scDRS/")
single_data<-readRDS("D:/Project/SingleCell_MDD/SingleCell_analysis/1_Singlecell_data/simple.rds")

############# 读入scRDS score数据 ###################
scRDS_score = read.table( gzfile("MDD.full_score.gz"), header = T)

############# 绘图 ##########################
single_data$scRDS_score<-scRDS_score$norm_score
metadata<-single_data@meta.data
metadata_anno<-metadata[,c("anno","scRDS_score")]
metadata_anno1<-aggregate(metadata_anno$scRDS_score,by=list(metadata_anno$anno),FUN=median)
metadata_anno1$color<-c("#d5231d","#3777ac","#4ea64a","#8e4c99","#e88f18","#e47faf","#b698c5","#a05528")
metadata_anno1<-metadata_anno1[order(metadata_anno1$x,decreasing = T),]
metadata_anno$anno<-factor(metadata_anno$anno,levels = metadata_anno1$Group.1)

#############boxplot###############
pdf("scRDS_score_boxplot.pdf",height = 5,width = 6)
ggplot(metadata_anno, aes(x=anno, y=scRDS_score)) + 
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
  ylab("scRDS score")+
  theme(axis.text.x = element_text(size = 10,color = 'black',vjust = 1,hjust = 1,angle = 35), # 调整x轴坐标文字
        axis.text.y = element_text(size = 10,color = 'black'), # 调整y轴坐标文字
        legend.text = element_text(size = 10,color = 'black'), # 调整legend字体大小
        legend.title = element_text(size = 10,color = 'black'), # 调整legend title大小
        panel.grid.major =  element_blank(),
        panel.grid.minor = element_blank()
  )
dev.off()


################# Featureplot ######################
embedding<-single_data@reductions$umap@cell.embeddings
embedding<-cbind(embedding,scRDS_score$norm_score)
embedding<-as.data.frame(embedding)
colnames(embedding)[3]<-"scRDS_score"

pdf("scDRS_score_featureplot.pdf",height = 6)
ggplot(embedding,aes(x=UMAP_1,y=UMAP_2,color=scRDS_score))+
  geom_point(size=0.2,alpha=0.2)+
  theme_bw()+
  scale_colour_gradient2(low="#3A71AA",mid="white",high="#B22028",midpoint=0)+
  labs(color="scRDS score")+
  theme_cowplot()
dev.off()

################# density plot
pdf("scRDS_score_density_plot.pdf",height = 4.5,width = 5)
Plot_Density_Custom(seurat_object =single_data, features = "scRDS_score",reduction = "umap",
                    custom_palette = c("#B0CFE4","#FACABC","#E77A77","#DC0000FF"))+
  theme(plot.title = element_text(hjust = 0.5))
dev.off()






```

# Step 2 correlation between scPagwas TRS score and scDRS score.R
```R
############ scPagwas score和scDRS score之间的相关性 ###################
setwd("D:/Project/SingleCell_MDD/GWAS_combine_singlecell_analysis/1_GWAS_analysis_for_broad_celltype/scDRS/")
library(Seurat)
library(ggplot2)
library(cowplot)
library(ggrepel)

single_data<-readRDS("D:/Project/SingleCell_MDD/SingleCell_analysis/1_Singlecell_data/simple.rds")

############ 按照细胞类型，对TRS求均值
TRS<-FetchData(single_data,vars = c("anno","scPagwas.TRS.Score"))
TRS1<-aggregate(TRS$scPagwas.TRS.Score,by=list(TRS$anno),FUN=median)
colnames(TRS1)<-c("cell_type","TRS")

############# 读入scRDS score数据 ###################
scRDS_score = read.table(gzfile("MDD.full_score.gz"), header = T)
write.csv(scRDS_score[,1:6],file = "scDRS_singlecell_result.csv",quote = F)
single_data$scRDS_score<-scRDS_score$norm_score

scDRS<-FetchData(single_data,vars = c("anno","scRDS_score"))
scDRS1<-aggregate(scDRS$scRDS_score,by=list(scDRS$anno),FUN=median)
colnames(scDRS1)<-c("cell_type","scDRS")

############  合并数据 #############
temp<-cbind(TRS1,scDRS1)
temp[,3]<-NULL
cor.test(temp$TRS,temp$scDRS) 

############# 绘图 ##########################
library(ggplot2)
library(cowplot)
library(ggpubr)

pdf("correlation_between_scPagwas_TRS_score_and_scDRS_score.pdf",height = 5,width = 6)
ggplot(temp,aes(TRS,scDRS))+
  geom_point(color="#2BA5B9",size=3.5)+
  geom_smooth(method = 'lm', formula = y ~ x, se = T,color="#D8AD63")+
  stat_cor(data=temp, method = "pearson")+   #相关性检验的R包
  theme_classic()+
  theme(panel.grid = element_blank(),
        plot.title = element_text (hjust = 0.5))+
  xlab("scPagwas TRS score")+
  ylab("scDRS score")+
  geom_text_repel(aes(label = cell_type),size = 3.5)+
  theme(axis.text.x = element_text(size = 10,color = 'black',vjust = 0.5), # 调整x轴坐标文字
    axis.text.y = element_text(size = 10,color = 'black'), # 调整y轴坐标文字
    legend.text = element_text(size = 10,color = 'black'), # 调整legend字体大小
    legend.title = element_text(size = 10,color = 'black'), # 调整legend title大小
    panel.grid.major =  element_blank(),
    plot.margin=unit(c(1,1,1,1),'lines'),#设置边距
    legend.position = 'right')
dev.off()



############### plot celltypes result
############### import the celltypes result
celltype_result<-read.table("MDD.scdrs_ct.anno",header = T,sep = "\t")
celltype_result<-celltype_result[,c("X","assoc_mcp")]

############### add another column for Significance
celltype_result[celltype_result$assoc_mcp<0.05,3]<-1
celltype_result[celltype_result$assoc_mcp>0.05,3]<-0
colnames(celltype_result)[3]<-"Significance"
celltype_result$Significance<-factor(celltype_result$Significance)

celltype_result$X<-factor(celltype_result$X,levels=rev(c("Excitatory neurons","Inhibitory neurons","Purkinje neurons","Endothelial cells",
                                                                 "Oligodendrocytes","OPCs","Microglia","Astrocytes")))


pdf("D:/Project/SingleCell_MDD/GWAS_combine_singlecell_analysis/1_GWAS_analysis_for_broad_celltype/ieu_b_102_scDRS.pdf",height = 4,width = 6)
ggplot(celltype_result,aes(X,-log10(assoc_mcp)))+
  geom_col(aes(fill=Significance),width = 0.8)+
  scale_fill_manual(values = c("#346EA5","#F69311"))+  ##手动指定颜色
  theme_bw()+
  theme(panel.grid = element_blank(),
        plot.title = element_text (hjust = 0.5),
        axis.text.x = element_text(angle = 90,vjust = 0.85,hjust = 0.75),
        legend.position = "none")+##调整字体角度
  xlab("Cell type")+
  ylab("-log10(P_value)")+
  theme_cowplot()+
  coord_flip()+  ##x，y轴互换
  scale_y_continuous(expand=c(0,0),limits = c (0, 10), breaks = seq (0, 10, 2))+
  geom_hline(aes(yintercept=(-log10(0.05))), colour="#F1F2F1", linetype="dashed") ##y轴上加刻度线
dev.off()




```












