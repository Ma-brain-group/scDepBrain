############# NC Ex-L2/4 pyscenic ##################
library(Seurat)
library(AUCell)
library(reshape2)
library(pheatmap)
library(ggplot2)
library(grid)
library(gtable)
library(ggpubr)
setwd("D:/Project/SingleCell_MDD/SingleCell_analysis/Exitatory_neurons_subtype_analysis/NC_pyscenic/")

################# import the NC Ex-L2/4 singlecell data
single_data<-readRDS("D:/Project/SingleCell_MDD/SingleCell_analysis/Exitatory_neurons_subtype_analysis/NC_Excitatory_neurons_annotation/NC_NEGR1_EX_L2_4.rds")
table(single_data$NEGR1_type)
Idents(single_data)<-single_data$NEGR1_type
single_data<-subset(single_data,idents = c("NEGR1+","NEGR1-"))

################# save the count matrix
write.csv(t(as.matrix(single_data@assays$RNA@counts)),file = "SingleCell_count.csv")

################# import the pyscenic result
auc<-read.csv("SingleCell.out.auc_mtx.csv",header = T,row.names = 1)
auc<-t(auc)
rownames(auc)<-gsub('\\...','', rownames(auc))

################# get the phenotype and NEGR1_type
single_data$NEGR1_type<-paste(single_data$NEGR1_type,single_data$phenotype,sep = "_")


################# heatmap Control group
cellTypes1<-subset(single_data,NEGR1_type=="NEGR1+_Control")
matrix1<-auc[,colnames(auc)%in%colnames(cellTypes1)]
cellTypes2<-subset(single_data,NEGR1_type=="NEGR1-_Control")
matrix2<-auc[,colnames(auc)%in%colnames(cellTypes2)]
matrix_AUC<-cbind(matrix1,matrix2)

################# scale the auc matrix
matrix_AUC<-log10(matrix_AUC+1)
matrix_AUC_Scaled <- t(scale(t(matrix_AUC),center = T, scale=T)) 


############## visualize the result
range(matrix_AUC_Scaled)
bk<-c(seq(-4,-0.1,by=0.01),seq(0,4,by=0.01))

############## construct the column annotation
cellTypes<-data.frame(celltype=c(rep("NEGR1+_Control",2065),rep("NEGR1-_Control",1715)))
row.names(cellTypes)<-colnames(matrix_AUC_Scaled)

p<-pheatmap(matrix_AUC_Scaled,cluster_cols = F,cluster_rows = T,
            show_rownames = T,show_colnames = F,
            annotation_col = cellTypes,
            cellheight = 10,
            cellwidth = 0.05,
            color = c(colorRampPalette(colors = c("#3A71AA","white"))(length(bk)/2),colorRampPalette(colors = c("white","red"))(length(bk)/2)),
            breaks = bk)
ggsave(filename = "NC_Ex-L-2-4_Control_NEGR1.pdf",p,height = 5,width = 6)



################# heatmap Case group
table(single_data$NEGR1_type)
cellTypes1<-subset(single_data,NEGR1_type=="NEGR1+_Case")
matrix1<-auc[,colnames(auc)%in%colnames(cellTypes1)]
cellTypes2<-subset(single_data,NEGR1_type=="NEGR1-_Case")
matrix2<-auc[,colnames(auc)%in%colnames(cellTypes2)]
matrix_AUC<-cbind(matrix1,matrix2)

################# scale the auc matrix
matrix_AUC<-log10(matrix_AUC+1)
matrix_AUC_Scaled <- t(scale(t(matrix_AUC),center = T, scale=T)) 


############## visualize the result
range(matrix_AUC_Scaled)
bk<-c(seq(-4,-0.1,by=0.01),seq(0,4,by=0.01))

############## construct the column annotation
cellTypes<-data.frame(celltype=c(rep("NEGR1+_Case",2859),rep("NEGR1-_Case",3207)))
row.names(cellTypes)<-colnames(matrix_AUC_Scaled)

p<-pheatmap(matrix_AUC_Scaled,cluster_cols = F,cluster_rows = T,
            show_rownames = T,show_colnames = F,
            annotation_col = cellTypes,
            cellheight = 10,
            cellwidth = 0.05,
            color = c(colorRampPalette(colors = c("#3A71AA","white"))(length(bk)/2),colorRampPalette(colors = c("white","red"))(length(bk)/2)),
            breaks = bk)
ggsave(filename = "NC_Ex-L-2-4_Case_NEGR1.pdf",p,height = 5,width = 7.5)








