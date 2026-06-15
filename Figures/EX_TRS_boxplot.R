############### IN TRS plot ##################
library(Seurat)
library(Nebulosa)
library(BiocFileCache)
library(paletteer)
library(scCustomize)
setwd("D:/Project/SingleCell_MDD/SingleCell_analysis/5_EX_TRS_plot/")
single_data_EX<-readRDS("D:/Project/SingleCell_MDD/SingleCell_analysis/annotation/Excitatory_neuron_subtypes_reannotation.rds")

############## 读入scPagwas的结果 #############
#load("D:/Project/SingleCell_MDD/GWAS_data/scPagwas/TRS_plot/scPagwas_singlecell_result_block_annotation_hg37_1/ieu_b_102/ieu_b_102_result.Rdata")
############# 提取出抑制性神经元的结果
#scPagwas_result<-a[match(colnames(single_data_EX),row.names(a)),]

############ 将结果加入单细胞数据中
#single_data_EX$scPagwas.gPAS.score<-scPagwas_result$scPagwas.gPAS.score
#single_data_EX$scPagwas.TRS.Score<-scPagwas_result$scPagwas.TRS.Score
#single_data_EX$Random_Correct_BG_adjp<-scPagwas_result$Random_Correct_BG_adjp

########### 提取出对应的信息，绘图 #################
metadata<-single_data_EX@meta.data
metadata<-metadata[,c("Excitatory_neuron_subtypes","scPagwas.TRS.Score")]
metadata1<-aggregate(metadata$scPagwas.TRS.Score,by=list(metadata$Excitatory_neuron_subtypes),FUN=median)
metadata1$colors<-c("#f9766e","#fbbab6","#e1c548","#5fa664","#abd0a7","#ca6a6b","#e5b5b5","#bac4d0")
metadata1<-metadata1[order(metadata1$x,decreasing = T),]
metadata$Excitatory_neuron_subtypes<-factor(metadata$Excitatory_neuron_subtypes,levels = metadata1$Group.1)

########## boxplot ##################
library(ggplot2)
library(cowplot)

pdf("EX_TRS_boxplot.pdf",height = 5,width = 6)
ggplot(metadata, aes(x=Excitatory_neuron_subtypes, y=scPagwas.TRS.Score)) + 
  stat_boxplot(geom = "errorbar",width=0.05, size=0.5,position=position_dodge(0.6),color= "black")+
  theme_bw()+
  geom_boxplot(position = position_dodge(0.6),
               size = 0.5,
               width = 0.85,
               fill = metadata1$colors,
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
        panel.grid.minor = element_blank())
dev.off()

############ feature plot ###################
embedding<-single_data_EX@reductions$umap@cell.embeddings
embedding<-cbind(embedding,metadata$scPagwas.TRS.Score)
embedding<-as.data.frame(embedding)
colnames(embedding)[3]<-"scPagwas.TRS.Score"

pdf("EX_TRS_featureplot.pdf",height = 6)
ggplot(embedding,aes(x=UMAP_1,y=UMAP_2,color=scPagwas.TRS.Score))+
  geom_point(size=0.2,alpha=0.2)+
  #theme_classic()+
  scale_color_distiller(palette = "Spectral")+
  labs(color="scPagwas.TRS.Score")+
  theme_cowplot()
dev.off()

############# 保存单细胞数据 #############
#saveRDS(single_data_EX,file = "EX.rds")


pdf("scPagwas.TRS.Score_density_plot.pdf",height = 4.5,width = 5)
Plot_Density_Custom(seurat_object =single_data_EX, features = "scPagwas.TRS.Score",reduction = "umap",
                    custom_palette = c("#B0CFE4","#FACABC","#E77A77","#DC0000FF"))+
  theme(plot.title = element_text(hjust = 0.5))
dev.off()





