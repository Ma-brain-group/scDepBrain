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
setwd("D:/Project/SingleCell_MDD/SingleCell_analysis/12_In6_pyscenic/")
single_data<-readRDS("D:/Project/SingleCell_MDD/SingleCell_analysis/11_CNNM2_positive_vs_CNNM2_negnative_IN6/CNNM2_IN6.rds")
table(single_data$CNNM2_type)

##############  remove the mitochondrial and ribosomal genes
#temp<-data.frame(as.matrix(single_data@assays$RNA@counts))
#temp<-temp[-grep("^MT-", rownames(temp)),]  
#temp<-temp[-grep("^RPS", rownames(temp)),]  
#temp<-temp[-grep("^RPL", rownames(temp)),] 

############## save the count matrix
write.csv(t(as.matrix(single_data@assays$RNA@counts)),file = "SingleCell_count.csv")


############## import the pyscenic result
auc<-read.csv("SingleCell.out.auc_mtx.csv",header = T,row.names = 1)
auc<-t(auc)
rownames(auc)<-gsub('\\...','', rownames(auc))
colnames(auc)<-gsub('\\.','-', colnames(auc))

############## adjust the cell order
Idents(single_data)<-single_data$CNNM2_type
cellTypes1<-subset(single_data,CNNM2_type=="CNNM2+")
matrix1<-auc[,colnames(auc)%in%colnames(cellTypes1)]
cellTypes2<-subset(single_data,CNNM2_type=="CNNM2-")
matrix2<-auc[,colnames(auc)%in%colnames(cellTypes2)]
matrix_AUC<-cbind(matrix1,matrix2)

############## reprocess the auc matrix
matrix_AUC<-log10(matrix_AUC+1)
matrix_AUC_Scaled <- t(scale(t(matrix_AUC),center = T, scale=T)) 

############## visualize the result
range(matrix_AUC_Scaled)
bk<-c(seq(-5,-0.1,by=0.01),seq(0,5,by=0.01))

############## construct the column annotation
cellTypes<-data.frame(celltype=c(rep("CNNM2+",12269),rep("CNNM2-",4197)))
row.names(cellTypes)<-colnames(matrix_AUC_Scaled)

p<-pheatmap(matrix_AUC_Scaled,cluster_cols = F,cluster_rows = T,
            show_rownames = T,show_colnames = F,
            annotation_col = cellTypes,
            #cellheight = 20,
            #cellwidth = 0.04,
            color = c(colorRampPalette(colors = c("#3A71AA","white"))(length(bk)/2),colorRampPalette(colors = c("white","red"))(length(bk)/2)),
            breaks = bk)
ggsave(filename = "In6_TF_activity.pdf",p,width = 6,height = 5)



############### plot the boxplot of CREB1 and THRB in different celltypes 
matrix_AUC_Scaled1<-matrix_AUC_Scaled[row.names(matrix_AUC_Scaled)%in%c("CREB1","THRB"),]
matrix_AUC_Scaled1<-as.data.frame(t(matrix_AUC_Scaled1))
matrix_AUC_Scaled1$cell_type<-c(rep("CNNM2+",12269),rep("CNNM2-",4197))
matrix_AUC_Scaled1$cell_type<-factor(matrix_AUC_Scaled1$cell_type,levels = c("CNNM2+","CNNM2-"))
groups <- list(c("CNNM2+","CNNM2-"))

pdf("The_violinplot_AUC_scores_of_CREB1.pdf",height = 4.2,width = 5.5)
ggplot(matrix_AUC_Scaled1, aes(x=cell_type,y=CREB1,fill=cell_type)) + 
  geom_violin(trim=F,color="black",width=1, position=position_dodge(0.75)) + 
  geom_boxplot(width=0.1,
               position=position_dodge(0.75),
               color="black",
               outlier.color = "black",
               outlier.fill = "black",
               outlier.shape = 19,
               outlier.size = 0.1,
               outlier.stroke = 0.5,
               outlier.alpha = 45,
               notch = F,
               notchwidth = 0.5)+ 
  scale_fill_manual(values = c("#e97371","#5ac6e9"))+
  theme_bw()+ 
  theme(plot.title = element_text(size=12,hjust=0.5), 
        axis.text.x = element_text(size = 11,color = 'black',hjust = 0.5,angle = 0), 
        axis.text.y = element_text(size = 11,color = 'black'),
        legend.text = element_text(size = 11,color = 'black'), 
        legend.title = element_text(size = 11,color = 'black'), 
        panel.grid.major = element_blank(),   
        panel.grid.minor = element_blank())+  
  xlab("")+
  ylab("The scaled activities of CREB1")+
  stat_compare_means(comparisons=groups,
                     method = "t.test",
                     label = "p.signif",
                     size=4.5)
dev.off()


pdf("The_violinplot_AUC_scores_of_THRB.pdf",height = 4.2,width = 5.5)
ggplot(matrix_AUC_Scaled1, aes(x=cell_type,y=THRB,fill=cell_type)) + 
  geom_violin(trim=F,color="black",width=1, position=position_dodge(0.75)) + 
  geom_boxplot(width=0.1,
               position=position_dodge(0.75),
               color="black",
               outlier.color = "black",
               outlier.fill = "black",
               outlier.shape = 19,
               outlier.size = 0.1,
               outlier.stroke = 0.5,
               outlier.alpha = 45,
               notch = F,
               notchwidth = 0.5)+ 
  scale_fill_manual(values = c("#e97371","#5ac6e9"))+
  theme_bw()+ 
  theme(plot.title = element_text(size=12,hjust=0.5), 
        axis.text.x = element_text(size = 11,color = 'black',hjust = 0.5,angle = 0), 
        axis.text.y = element_text(size = 11,color = 'black'),
        legend.text = element_text(size = 11,color = 'black'), 
        legend.title = element_text(size = 11,color = 'black'), 
        panel.grid.major = element_blank(),   
        panel.grid.minor = element_blank())+  
  xlab("")+
  ylab("The scaled activities of THRB")+
  stat_compare_means(comparisons=groups,
                     method = "t.test",
                     label = "p.signif",
                     size=4.5)
dev.off()








