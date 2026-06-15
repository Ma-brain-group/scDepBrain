################# nature communication singlecell data reannotation ####################
library(Seurat)
library(ggplot2)
library(harmony)
library(dplyr)
library(ggpubr)
setwd("D:/Project/SingleCell_MDD/Figure/NC_In_PVALB/annotation/")

################# import the nc singlecell data
single_data_INs<-readRDS("D:/Project/SingleCell_MDD/SingleCell_analysis/16_validation_NC/2_INs_subset/NC_inhibitory_neurons_reannotation.rds")
single_data<-readRDS("D:/Project/SingleCell_MDD/SingleCell_analysis/16_validation_NC/1_NC_data_process/GSE213982/2023_NC.rds")

################# dimplot for celltypes 
color<-c("#e0bc58","#64abc0","#fab37f","#e98741","#f2d3ca","#eebd85")
pdf("DimPlot_NC_inhibitory_neurons_reannotation.pdf",height =4.3,width = 5)
DimPlot(single_data_INs,group.by = "RNA_snn_res.0.1_cell_type",reduction = "umap",cols = color,label = T,repel = T)                                                                                                       
dev.off()                                                                                                       


################# dotplot for celltypes                                                                                                       
pdf("DotPlot_NC_inhibitory_neurons_reannotation_marker_gene.pdf",height = 4.5,width = 4.5)
DotPlot(single_data_INs,features=rev(c("VIP", ## In_VIP
                                       "SST", ## In_SST
                                       "PVALB", ## In_PVALB
                                       "LAMP5", ## In_LAMP5
                                       "RELN" ## In_Reelin
                                       )),group.by="RNA_snn_res.0.1_cell_type",assay="RNA")+
  scale_colour_gradient2(low="#3A71AA",mid="white",high="#B22028",midpoint=0)+
  theme_bw()+
  coord_flip()+
  theme(axis.text.x = element_text (size = 10,color = 'black',angle = 35,vjust = 1, hjust = 1),
        axis.text.y = element_text(size = 10,color = 'black'),
        legend.text = element_text(size = 10,color = 'black'),
        legend.title = element_text(size = 10,color = 'black'))
dev.off()


############### dotplot for CNNM2 in Inhibitory neurons subtypes
pdf("DotPlot_NC_CNNM2_in_inhibitory_neurons_dotplot.pdf",height = 3.5,width =5.5)
DotPlot(single_data_INs,features="CNNM2",group.by="RNA_snn_res.0.1_cell_type",assay="RNA")+
  scale_colour_gradient2(low="#3A71AA",mid="white",high="#B22028",midpoint=-0.5)+
  theme_bw()+
  coord_flip()+
  theme(axis.text.x = element_text (size = 10,color = 'black',angle = 35,vjust = 1, hjust = 1),
        axis.text.y = element_text(size = 10,color = 'black'),
        legend.text = element_text(size = 10,color = 'black'),
        legend.title = element_text(size = 10,color = 'black'))
dev.off()



############### boxplot for CNNM2 between In_PVALB and other Inhibitory neurons subtypes
expression<-single_data_INs@assays$RNA@counts
expression<-expression[row.names(expression)%in%c("CNNM2"),]
expression<-as.data.frame(as.matrix(expression))
colnames(expression)<-"CNNM2_expression"
expression$cell_type<-single_data_INs$RNA_snn_res.0.1_cell_type
expression$cell_type<-as.character(expression$cell_type)
expression[expression$cell_type!="In_PVALB","cell_type"]<-"Other Inhibitory neurons"

### log10(expression+1)
expression$log_expression<-log10(expression$CNNM2_expression+1)

### calculate mean expression
expression1<-aggregate(expression$log_expression,by=list(expression$cell_type),mean)
colnames(expression1)<-c("cell_type","CNNM2_mean_expression")
write.csv(expression1,file = "CNNM2_mean_expression.csv",quote = F)

### plot boxplot
mycomparisons<-list(c("In_PVALB","Other Inhibitory neurons"))
pdf("boxplot_CNNM2_expression_in_Inhibitory_neurons_subtypes.pdf",height = 4,width = 6)
ggplot(expression, aes(x=cell_type,y=log_expression,fill=cell_type)) + 
  stat_boxplot(geom = "errorbar",width=0.05, size=0.5,position=position_dodge(0.6),color= "black")+
  theme_bw()+
  geom_boxplot(position = position_dodge(0.6),
               size = 0.5,
               width = 0.5,
               color = "black",
               outlier.color = "black",
               outlier.fill = "black",
               outlier.shape = 19,
               outlier.size = 1,
               outlier.stroke = 0.5,
               outlier.alpha = 45,
               notch = F,
               notchwidth = 0.5)+
  scale_fill_manual(values = c("#e97371","#5ac6e9"))+
  xlab("")+
  ylab("CNNM2 expression")+
  theme(axis.text.x = element_text(size = 10,color = 'black',vjust = 0.5,hjust=0.5,angle = 0), 
        axis.text.y = element_text(size = 10,color = 'black'), 
        legend.text = element_text(size = 10,color = 'black'), 
        legend.title = element_text(size = 10,color = 'black'), 
        plot.title = element_text(size=12,hjust=0.5), 
        panel.grid.major =  element_blank(),
        panel.grid.minor = element_blank())+
  stat_compare_means(comparisons = mycomparisons,
                     method = "wilcox.test",
                     label = "p.signif",
                     size=4.5)
dev.off()


############### dotplot for CNNM2 in broad celltypes
pdf("DotPlot_NC_CNNM2_in_broad_cell_type_dotplot.pdf",height = 3.5,width =5.5)
DotPlot(single_data,features="CNNM2",group.by="broad_cell_type",assay="RNA")+
  scale_colour_gradient2(low="#3A71AA",mid="white",high="#B22028",midpoint=0)+
  theme_bw()+
  coord_flip()+
  theme(axis.text.x = element_text (size = 10,color = 'black',angle = 35,vjust = 1, hjust = 1),
        axis.text.y = element_text(size = 10,color = 'black'),
        legend.text = element_text(size = 10,color = 'black'),
        legend.title = element_text(size = 10,color = 'black'))
dev.off()





