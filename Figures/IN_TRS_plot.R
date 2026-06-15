############### IN TRS plot ##################
library(Seurat)
library(Nebulosa)
library(BiocFileCache)
library(paletteer)
library(scCustomize)
setwd("D:/Project/SingleCell_MDD/SingleCell_analysis/6_IN_TRS_plot/")
single_data_IN<-readRDS("D:/Project/SingleCell_MDD/SingleCell_analysis/annotation/Inhibitory_neuron_subtypes_reannotation.rds")

############## 读入scPagwas的结果 #############
#load("D:/Project/SingleCell_MDD/GWAS_data/scPagwas/TRS_plot/scPagwas_singlecell_result_block_annotation_hg37_1/ieu_b_102/ieu_b_102_result.Rdata")
############# 提取出抑制性神经元的结果
#scPagwas_result<-a[match(colnames(single_data_IN),row.names(a)),]

############ 将结果加入单细胞数据中
#single_data_IN$scPagwas.gPAS.score<-scPagwas_result$scPagwas.gPAS.score
#single_data_IN$scPagwas.TRS.Score<-scPagwas_result$scPagwas.TRS.Score
#single_data_IN$Random_Correct_BG_adjp<-scPagwas_result$Random_Correct_BG_adjp

########### 提取出对应的信息，绘图 #################
metadata<-single_data_IN@meta.data
metadata<-metadata[,c("Inhibitory_neuron_subtypes","scPagwas.TRS.Score")]
metadata1<-aggregate(metadata$scPagwas.TRS.Score,by=list(metadata$Inhibitory_neuron_subtypes),FUN=median)
metadata1$color<-c("#e0bc58","#64abc0","#fab37f","#e98741","#8fc0dc","#967568")
metadata1<-metadata1[order(metadata1$x,decreasing = T),]
metadata$Inhibitory_neuron_subtypes<-factor(metadata$Inhibitory_neuron_subtypes,levels = metadata1$Group.1)

########## boxplot ##################
library(ggplot2)
library(cowplot)

pdf("IN_TRS_boxplot.pdf",height = 5,width = 6)
ggplot(metadata, aes(x=Inhibitory_neuron_subtypes, y=scPagwas.TRS.Score)) + 
  stat_boxplot(geom = "errorbar",width=0.05, size=0.5,position=position_dodge(0.6),color= "black")+
  theme_bw()+
  geom_boxplot(position = position_dodge(0.6),
               size = 0.5,
               width = 0.8,
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
  theme(axis.text.x = element_text(size = 10,color = 'black',vjust = 1,hjust=1,angle = 35), # 调整x轴坐标文字
        axis.text.y = element_text(size = 10,color = 'black'), # 调整y轴坐标文字
        legend.text = element_text(size = 10,color = 'black'), # 调整legend字体大小
        legend.title = element_text(size = 10,color = 'black'), # 调整legend title大小
        panel.grid.major =  element_blank(),
        panel.grid.minor = element_blank())
dev.off()

############ feature plot ###################
embedding<-single_data_IN@reductions$umap@cell.embeddings
embedding<-cbind(embedding,metadata$scPagwas.TRS.Score)
embedding<-as.data.frame(embedding)
colnames(embedding)[3]<-"scPagwas.TRS.Score"

pdf("IN_TRS_featureplot.pdf",height = 6)
ggplot(embedding,aes(x=UMAP_1,y=UMAP_2,color=scPagwas.TRS.Score))+
  geom_point(size=0.2,alpha=0.2)+
  #theme_bw()+
  scale_colour_gradient2(low="blue",mid="white",high="red",midpoint=0.58)+
  labs(color="scPagwas.TRS.Score")+
  theme_cowplot()
dev.off()


############## density plot
pdf("scPagwas.TRS.Score_density_plot.pdf",height = 4.5,width = 5)
Plot_Density_Custom(seurat_object =single_data_IN, features = "scPagwas.TRS.Score",reduction = "umap",
                    custom_palette = c("white","#B0CFE4","#FACABC","#E77A77","#DC0000FF"))+
  theme(plot.title = element_text(hjust = 0.5))
dev.off()





############# 保存单细胞数据 #############
#saveRDS(single_data_IN,file = "IN.rds")








