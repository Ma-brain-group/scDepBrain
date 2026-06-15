################# scPagwas plot for Development singlecell data ################
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
setwd("D:/Project/SingleCell_MDD/Human_million_Brain_singlecell/Development_data/scPagwas/")

################# import the development singlecell data
#### 60w development brain singlecell
single_data<-readRDS("D:/Project/SingleCell_MDD/Human_million_Brain_singlecell/Development_data/annotation/clean.Dev_metaatlas_integrated.rds")
DefaultAssay(single_data)<-"RNA"
single_data[["integrated"]]<-NULL
Idents(single_data)<-single_data$Type.v2
table(single_data$Type.v2)

################# import the scPagwas result
scPagwas<-read.csv("singlecell_scPagwas_score_pvalue_Result.csv",row.names = 1)
identical(row.names(scPagwas),colnames(single_data))

################# get the metadata
metadata<-single_data@meta.data
metadata1<-metadata[!is.na(metadata$Type.v2),]
identical(row.names(scPagwas),row.names(metadata1))
metadata1$scPagwas.TRS.Score<-scPagwas$scPagwas.TRS.Score
metadata1<-metadata1[,c("Type.v2","scPagwas.TRS.Score")]

################# boxplot
metadata2<-aggregate(metadata1$scPagwas.TRS.Score,by=list(metadata1$Type.v2),median)
metadata2$color<-c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
                   '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
                   '#9FA3A8')
metadata2<-metadata2[order(metadata2$x,decreasing = T),]
metadata1$Type.v2<-factor(metadata1$Type.v2,levels = metadata2$Group.1)

pdf("boxplot_scPagwas_score_for_Type.v2.pdf",height = 5,width = 7)
ggplot(metadata1, aes(x=Type.v2, y=scPagwas.TRS.Score)) + 
  stat_boxplot(geom = "errorbar",width=0.05, size=0.5,position=position_dodge(0.6),color= "black")+
  theme_bw()+
  geom_boxplot(position = position_dodge(0.6),
               size = 0.5,
               width = 0.85,
               fill = metadata2$color,
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


################ scPagwas celltypes result
scPagwas_celltype<-read.csv("celltype_result.csv",row.names = 1)
scPagwas_celltype<-scPagwas_celltype[order(scPagwas_celltype$pvalue),]
colnames(metadata2)[1]<-"celltype"
scPagwas_celltype<-merge(scPagwas_celltype,metadata2,by="celltype")
scPagwas_celltype$trans_P<--log10(scPagwas_celltype$pvalue)
scPagwas_celltype<-scPagwas_celltype[order(scPagwas_celltype$trans_P,decreasing = T),]
scPagwas_celltype$celltype<-factor(scPagwas_celltype$celltype,levels = rev(scPagwas_celltype$celltype))

pdf("ggdotchart_for_celltype_result.pdf",height = 4.6,width = 5.5)
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
  scale_y_continuous(limits = c (0, 6), breaks = seq (0, 6, 2))
dev.off()


################## density plot for scPagwas score
###### remove the cells without annotation
single_data<-readRDS("D:/Project/SingleCell_MDD/Human_million_Brain_singlecell/Development_data/annotation/clean.Dev_metaatlas_integrated.rds")
DefaultAssay(single_data)<-"RNA"
single_data[["integrated"]]<-NULL
Idents(single_data)<-single_data$Type.v2
single_data<-subset(single_data,idents = c("Astrocyte","Div","Endothelial","Excitatory Neuron","Inhibitory Neuron","IPC","Microglia",
                                           "Newborn Neuron","Oligodendrocyte","OPC","OPC.div","Pericyte","Red blood cells","RG","RG.div"))

###### import the scPagwas score
scPagwas<-read.csv("singlecell_scPagwas_score_pvalue_Result.csv",row.names = 1)
identical(row.names(scPagwas),colnames(single_data))
single_data$scPagwas.TRS.Score<-scPagwas$scPagwas.TRS.Score

###### density plot 
pdf("Density_plot_for_scPagwas_TRS_Score.pdf",height = 6,width = 6)
Plot_Density_Custom(seurat_object =single_data, features = "scPagwas.TRS.Score",reduction = "umap",
                    custom_palette = c("blue","#57C3F3","#B0CFE4","#CCE0F5"))+
  theme(plot.title = element_text(hjust = 0.5))
dev.off()




