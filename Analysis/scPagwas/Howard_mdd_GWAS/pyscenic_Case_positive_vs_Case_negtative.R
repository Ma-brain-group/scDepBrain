############## pyscenic input ##################
library(Seurat)
library(AUCell)
library(reshape2)
library(pheatmap)
library(ggplot2)
library(grid)
library(gtable)
library(ggpubr)

############## iuput the singlecell data
setwd("D:/Project/SingleCell_MDD/SingleCell_analysis/16_validation_NC/pyscenic_Case_positive_vs_Case_negtative/")
single_data<-readRDS("D:/Project/SingleCell_MDD/SingleCell_analysis/16_validation_NC/2_INs_subset/NC_In_PVALB.rds")
table(single_data$CNNM2_type)

############## subset the CNNM2+_Case, CNNM2-_Case cells
Idents(single_data)<-single_data$CNNM2_type
single_data<-subset(single_data,idents=c("CNNM2+_Case","CNNM2-_Case"))

############## save the count matrix
write.csv(t(as.matrix(single_data@assays$RNA@counts)),file = "SingleCell_count.csv")


############## import the pyscenic result
auc<-read.csv("SingleCell.out.auc_mtx.csv",header = T,row.names = 1)
auc<-t(auc)
rownames(auc)<-gsub('\\...','', rownames(auc))


############## adjust the cell order
cellTypes1<-subset(single_data,CNNM2_type=="CNNM2+_Case")
matrix1<-auc[,colnames(auc)%in%colnames(cellTypes1)]
cellTypes2<-subset(single_data,CNNM2_type=="CNNM2-_Case")
matrix2<-auc[,colnames(auc)%in%colnames(cellTypes2)]
matrix_AUC<-cbind(matrix1,matrix2)

############## reprocess the auc matrix
matrix_AUC<-log10(matrix_AUC+1)
matrix_AUC_Scaled <- t(scale(t(matrix_AUC),center = T, scale=T)) 


############## visualize the result
range(matrix_AUC_Scaled)
bk<-c(seq(-6,-0.1,by=0.01),seq(0,6,by=0.01))

############## construct the column annotation
cellTypes<-data.frame(celltype=c(rep("CNNM2+_Case",2365),rep("CNNM2-_Case",2519)))
row.names(cellTypes)<-colnames(matrix_AUC_Scaled)

p<-pheatmap(matrix_AUC_Scaled,cluster_cols = F,cluster_rows = T,
            show_rownames = T,show_colnames = F,
            annotation_col = cellTypes,
            cellheight = 20,
            cellwidth = 0.04,
            color = c(colorRampPalette(colors = c("#3A71AA","white"))(length(bk)/2),colorRampPalette(colors = c("white","red"))(length(bk)/2)),
            breaks = bk)
ggsave(filename = "pyscenic_Case_positive_vs_Case_negtative.pdf",p,height = 5,width = 6)



############### plot the boxplot of THRB in different celltypes 
matrix_AUC_Scaled1<-matrix_AUC_Scaled[row.names(matrix_AUC_Scaled)%in%c("THRB"),]
matrix_AUC_Scaled1<-as.data.frame(matrix_AUC_Scaled1)
matrix_AUC_Scaled1$cell_type<-c(rep("CNNM2+_Case",2365),rep("CNNM2-_Case",2519))
matrix_AUC_Scaled1$cell_type<-factor(matrix_AUC_Scaled1$cell_type,levels = c("CNNM2+_Case","CNNM2-_Case"))
groups <- list(c("CNNM2+_Case","CNNM2-_Case"))

pdf("The_AUC_scores_of_THRB.pdf",height = 4.2,width = 4)
ggplot(matrix_AUC_Scaled1, aes(x=cell_type, y=matrix_AUC_Scaled1)) + 
  stat_boxplot(geom = "errorbar",width=0.05, size=0.5,position=position_dodge(0.6),color= "black")+
  theme_bw()+
  geom_boxplot(position = position_dodge(0.6),
               size = 0.5,
               width = 0.5,
               fill = c("#E77A77","#67ADB7"),
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
  ylab("The AUC scores of THRB")+
  theme(axis.text.x = element_text(size = 10,color = 'black',vjust = 1,hjust = 0.5,angle = 0), # 调整x轴坐标文字
        axis.text.y = element_text(size = 10,color = 'black'), # 调整y轴坐标文字
        legend.text = element_text(size = 10,color = 'black'), # 调整legend字体大小
        legend.title = element_text(size = 10,color = 'black'), # 调整legend title大小
        panel.grid.major =  element_blank(),
        panel.grid.minor = element_blank())+
  stat_compare_means(comparisons=groups,method="t.test")
dev.off()





