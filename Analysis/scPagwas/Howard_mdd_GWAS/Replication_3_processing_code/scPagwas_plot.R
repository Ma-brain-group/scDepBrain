################# Brain organoid scPagwas plot ####################
library(ggplot2)
library(Seurat)
library(ggplot2)
library(ggvenn) 
library(cowplot)
library(ggpubr)
library(Nebulosa)
library(BiocFileCache)
library(paletteer)
library(scCustomize)
setwd("D:/Project/SingleCell_MDD/organoid/scPagwas/")

################# import the singlecell data
single_data<-readRDS("D:/Project/SingleCell_MDD/organoid/Sample_merge/Brain_organoid.rds")
single_data
################# import the scPagwas result
scPagwas<-read.csv("singlecell_scPagwas_score_pvalue_Result.csv",row.names = 1)

################# import the metadata
identical(row.names(scPagwas),colnames(single_data))
single_data$scPagwas.TRS.Score<-scPagwas$scPagwas.TRS.Score
metadata<-FetchData(single_data,vars = c("cell_type","scPagwas.TRS.Score"))

################# boxplot for celltypes
metadata1<-aggregate(metadata$scPagwas.TRS.Score,by=list(metadata$cell_type),median)
metadata1$color<-c("#4ea64a","#8e4c99","#d5231d","#e88f18","#3777ac","#e47faf","#b698c5","#a05528","#B3446C")
metadata1<-metadata1[order(metadata1$x,decreasing = T),]
metadata$cell_type<-factor(metadata$cell_type,levels = metadata1$Group.1)

pdf("boxplot_scPagwas_score_for_cell_type.pdf",height = 5,width = 6)
ggplot(metadata, aes(x=cell_type, y=scPagwas.TRS.Score)) + 
  stat_boxplot(geom = "errorbar",width=0.05, size=0.5,position=position_dodge(0.6),color= "black")+
  theme_bw()+
  geom_boxplot(position = position_dodge(0.6),
               size = 0.5,
               width = 0.85,
               fill = metadata1$color,
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
  theme(axis.text.x = element_text(size = 10,color = 'black',vjust = 1,hjust = 1,angle = 35), 
        axis.text.y = element_text(size = 10,color = 'black'), 
        legend.text = element_text(size = 10,color = 'black'), 
        legend.title = element_text(size = 10,color = 'black'),
        panel.grid.major =  element_blank(),
        panel.grid.minor = element_blank())
dev.off()


############# density plot for scPagwas.TRS.Score
pdf("Density_plot_for_scPagwas_TRS_Score.pdf",height = 5.5,width = 5.5)
Plot_Density_Custom(seurat_object =single_data, features = "scPagwas.TRS.Score",reduction = "umap",
                    custom_palette = c("#B0CFE4","#FACABC","#E77A77","#DC0000FF"))+
  theme(plot.title = element_text(hjust = 0.5))
dev.off()


############# scPagwas celltypes result
scPagwas_celltype<-read.csv("celltype_result.csv",row.names = 1)
scPagwas_celltype<-scPagwas_celltype[order(scPagwas_celltype$pvalue),]
colnames(metadata1)[1]<-"celltype"
scPagwas_celltype<-merge(scPagwas_celltype,metadata1,by="celltype")
scPagwas_celltype$trans_P<--log10(scPagwas_celltype$pvalue)
scPagwas_celltype<-scPagwas_celltype[order(scPagwas_celltype$trans_P,decreasing = T),]
scPagwas_celltype$celltype<-factor(scPagwas_celltype$celltype,levels = rev(scPagwas_celltype$celltype))

pdf("ggdotchart_for_celltype_result.pdf",height = 4,width = 5)
ggdotchart(scPagwas_celltype, x = "celltype", y = "trans_P",
           color = "celltype",                              
           palette = rev(scPagwas_celltype$color),
           sorting = "descending",                       # Sort value in descending order
           add = "segments",                             # Add segments from y = 0 to dots
           add.params = list(color = "lightgray", size = 2), # Change segment color and size
           group = "celltype",                                # Order by groups
           dot.size = 7,                                 # Large dot size
           label = round(scPagwas_celltype$trans_P,3), # Add mpg values as dot labels
           font.label = list(color = "black", size = 8, vjust = 0.5),               # Adjust label parameters
           ggtheme = theme_pubr())+
  coord_flip()+
  geom_hline(yintercept = -log10(0.05), linetype = 2, color = "lightgray")+
  theme(axis.text.x = element_text (size = 10,color = 'black',angle = 0,hjust = 0.5),
        axis.text.y = element_text(size = 10,color = 'black'),
        legend.text = element_text(size = 10,color = 'black'),
        legend.title = element_text(size = 10,color = 'black'),
        legend.position = 'none'
  )+
  ylab("-log10(P)")+
  scale_y_continuous(limits = c (0, 10), breaks = seq (0, 10, 2))
dev.off()
