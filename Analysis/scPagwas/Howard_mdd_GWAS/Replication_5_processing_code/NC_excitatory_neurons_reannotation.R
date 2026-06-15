################# nature communication singlecell data reannotation ####################
library(Seurat)
library(ggplot2)
library(harmony)
library(dplyr)
library(ggpubr)
setwd("D:/Project/SingleCell_MDD/Figure/NC_Ex-L2-4/annotation/")

################# import the nc singlecell data
single_data_Ex<-readRDS("D:/Project/SingleCell_MDD/SingleCell_analysis/Exitatory_neurons_subtype_analysis/NC_Excitatory_neurons_annotation/NC_Excitatory_neuron_subtypes.rds")
single_data<-readRDS("D:/Project/SingleCell_MDD/SingleCell_analysis/16_validation_NC/1_NC_data_process/GSE213982/2023_NC.rds")


############### dotplot for NEGR1 in Excitatory neurons subtypes
pdf("DotPlot_NC_NEGR1_in_Excitatory_neurons_dotplot.pdf",height = 3.5,width =5.5)
DotPlot(single_data_Ex,features="NEGR1",group.by="Excitatory_neuron_subtypes",assay="RNA")+
  scale_colour_gradient2(low="#3A71AA",mid="white",high="#B22028",midpoint=0)+
  theme_bw()+
  coord_flip()+
  theme(axis.text.x = element_text (size = 10,color = 'black',angle = 35,vjust = 1, hjust = 1),
        axis.text.y = element_text(size = 10,color = 'black'),
        legend.text = element_text(size = 10,color = 'black'),
        legend.title = element_text(size = 10,color = 'black'))
dev.off()



############### boxplot for NEGR1 between Ex-L2/4 and other Exitatory neurons subtypes
expression<-single_data_Ex@assays$RNA@counts
expression<-expression[row.names(expression)%in%c("NEGR1"),]
expression<-as.data.frame(as.matrix(expression))
colnames(expression)<-"NEGR1_expression"
expression$cell_type<-single_data_Ex$Excitatory_neuron_subtypes
expression$cell_type<-as.character(expression$cell_type)
expression[expression$cell_type!="Ex-L2/4","cell_type"]<-"Other Excitatory neurons"

### log10(expression+1)
expression$log_expression<-log10(expression$NEGR1_expression+1)

### calculate mean expression
expression1<-aggregate(expression$log_expression,by=list(expression$cell_type),mean)
colnames(expression1)<-c("cell_type","NEGR1_mean_expression")
write.csv(expression1,file = "NEGR1_mean_expression.csv",quote = F)

### plot boxplot
mycomparisons<-list(c("Ex-L2/4","Other Excitatory neurons"))
pdf("boxplot_NEGR1_expression_in_Excitatory_neurons_subtypes.pdf",height = 4,width = 6)
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
  ylab("NEGR1 expression")+
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
pdf("DotPlot_NC_NEGR1_in_broad_cell_type_dotplot.pdf",height = 3.5,width =5.5)
DotPlot(single_data,features="NEGR1",group.by="broad_cell_type",assay="RNA")+
  scale_colour_gradient2(low="#3A71AA",mid="white",high="#B22028",midpoint=0)+
  theme_bw()+
  coord_flip()+
  theme(axis.text.x = element_text (size = 10,color = 'black',angle = 35,vjust = 1, hjust = 1),
        axis.text.y = element_text(size = 10,color = 'black'),
        legend.text = element_text(size = 10,color = 'black'),
        legend.title = element_text(size = 10,color = 'black'))
dev.off()
