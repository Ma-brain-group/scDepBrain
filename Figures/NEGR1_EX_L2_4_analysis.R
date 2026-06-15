############# Excitatory neurons Ex-L2/4 five methods overlap gene ##############
library(Seurat)
library(SeuratObject)
library(ggplot2)

############# set the workspace path
setwd("D:/Project/SingleCell_MDD/SingleCell_analysis/Exitatory_neurons_subtype_analysis/Five_methods_overlap_gene/")

############# method1: the DEG between Ex-L2/4 and other broad celltypes #################
############# import the singlecell data
single_data<-readRDS("/share2/pub/chenchg/chenchg/SingleCell/Brain/SingleCell_data/singlecell_data/MDD_singlecell_data_reannotation.rds")
DEG<-FindMarkers(single_data,assay = "RNA",group.by = "final_anno",ident.1 = "Ex-L2/4",only.pos = T)
DEG<-DEG[DEG$p_val_adj<0.05,]
write.csv(DEG,file="FindMarkers_Ex-L2_4_vs_other_broad_celltypes_gene.csv",quote=F)
DEG<-read.csv("FindMarkers_Ex-L2_4_vs_other_broad_celltypes_gene.csv",row.names = 1)


############# method2: the DEG between Ex-L2/4 and other Excitatory neurons subtypes ####################
############# import the singlecell data
single_data_Excitatory_neurons<-readRDS("/share2/pub/chenchg/chenchg/SingleCell/Brain/SingleCell_data/singlecell_data/Excitatory_neuron_subtypes_reannotation.rds")
DEG1<-FindMarkers(single_data_Excitatory_neurons,assay = "RNA",group.by = "Excitatory_neuron_subtypes",ident.1 = "Ex-L2/4",only.pos = T)
DEG1<-DEG1[DEG1$p_val_adj<0.05,]
write.csv(DEG1,file="FindMarkers_Ex-L2_4_vs_Exs_gene.csv",quote=F)
DEG1<-read.csv("FindMarkers_Ex-L2_4_vs_Exs_gene.csv",row.names = 1)



############# method3: the DEG between TRS top10% and bottom10% ###################
single_data_Excitatory_neurons<-readRDS("D:/Project/SingleCell_MDD/SingleCell_analysis/annotation/Excitatory_neuron_subtypes_reannotation.rds")
DefaultAssay(single_data_Excitatory_neurons)<-"RNA"

############# subset the Ex-L2/4
Idents(single_data_Excitatory_neurons)<-single_data_Excitatory_neurons$Excitatory_neuron_subtypes
single_data_EX12<-subset(single_data_Excitatory_neurons,idents="Ex-L2/4")
single_data_EX12<-NormalizeData(single_data_EX12)

############# get the TRS top10% and bottom10% singlecell data
threshold_up<-quantile(single_data_EX12$scPagwas.TRS.Score,0.9)
threshold_down<-quantile(single_data_EX12$scPagwas.TRS.Score,0.1)
single_data_EX12@meta.data$TRS_type<-NULL
single_data_EX12@meta.data$TRS_type[single_data_EX12@meta.data$scPagwas.TRS.Score>=threshold_up]<-"positive"
single_data_EX12@meta.data$TRS_type[single_data_EX12@meta.data$scPagwas.TRS.Score<threshold_down]<-"negative"
Idents(single_data_EX12)<-single_data_EX12$TRS_type

############# 456 cells
single_data_EX12<-subset(single_data_EX12,idents = c("positive","negative"))
DEG2<-FindMarkers(single_data_EX12,assay = "RNA",group.by = "TRS_type",ident.1 = "positive",ident.2 = "negative",only.pos = T)
DEG2<-DEG2[DEG2$p_val_adj<0.05,]
write.csv(DEG2,file="FindMarkers_TRS_top10_vs_bottom10_gene.csv",quote=F)


############## method4: the result of magma 
magma_risk_gene<-read.csv("D:/Project/SingleCell_MDD/SingleCell_analysis/Manhattan_Plot/magma_fdr_0.05.csv",row.names = 1)


############## method5: the result of fusion 
fusion_risk_gene<-read.csv("D:/Project/SingleCell_MDD/GWAS_analysis/Fusion/fusion_result_Brain_Frontal_Cortex_BA9/Brain_Frontal_Cortex_BA9_result.txt_P_0.05.csv")
#fusion_risk_gene1<-read.csv("D:/Project/SingleCell_MDD/GWAS_analysis/Fusion/fusion_result_Brain_Cortex/Brain_Cortex_result.txt_P_0.05.csv")


############## four methods overlap genes
genes<-Reduce(intersect, list(row.names(DEG),row.names(DEG1),row.names(DEG2),magma_risk_gene$SYMBOL[1:30]))



############## visualize the result ####################
############## Figure1 the top 30 genes in magma
magma_risk_gene1<-magma_risk_gene[1:30,c("SYMBOL","fdr")]
magma_risk_gene1$`-log10(fdr)`<--log10(magma_risk_gene1$fdr)
magma_risk_gene1$Z<-"Z"
################## change the character to factor
magma_risk_gene1<-magma_risk_gene1[order(magma_risk_gene1$`-log10(fdr)`,decreasing = T),]
magma_risk_gene1$SYMBOL<-factor(magma_risk_gene1$SYMBOL,levels = rev(magma_risk_gene1$SYMBOL))

pdf("The_result_of_top30_magma_gene.pdf")
ggplot(magma_risk_gene1, aes(x=Z,y=SYMBOL,fill = -log10(fdr))) + 
  geom_tile(color="white",linetype = 1)+
  coord_equal()+
  theme_minimal()+
  scale_fill_gradient2(low="#F6A45D",midpoint = 4.8,high="#D53D4E") + 
  theme(axis.text.x = element_text(size = 10,color = 'black',angle=90,vjust = 1,hjust = 1), # 调整x轴坐标文字
        axis.text.y = element_text(size = 10,color = 'black'), # 调整y轴坐标文字
        legend.text = element_text(size = 10,color = 'black'),
        legend.title = element_text(size = 10,color = 'black'),
        panel.grid.major =  element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = 'right')+
  xlab("")+
  ylab("")
dev.off()



############## Figure2 the overlap gene between top 30 genes in magma and the DEG between EX12 and other broad celltypes
intersect_gene1<-intersect(magma_risk_gene$SYMBOL[1:30],row.names(DEG))

################## construct a dataframe
dataframe1<-data.frame(gene=magma_risk_gene1$SYMBOL,type=magma_risk_gene1$SYMBOL%in%intersect_gene1)
dataframe1[dataframe1$type=="TRUE",3]<-"1"
dataframe1[dataframe1$type=="FALSE",3]<-"0"
colnames(dataframe1)[3]<-"type1"
dataframe1$Z<-"Z"
dataframe1$gene<-factor(dataframe1$gene,levels = rev(magma_risk_gene1$SYMBOL))


pdf("The_result_of_the_overlap_gene_between_magma_top30_and_DEG_between_EX2_4_and_other_broad_celltypes.pdf")
ggplot(dataframe1, aes(x=Z,y=gene,fill = type1)) + 
  geom_tile(color="white",linetype = 1)+
  coord_equal()+
  theme_minimal()+
  scale_fill_manual(values = c("#F1F2F1","#ACAF37")) + 
  theme(axis.text.x = element_text(size = 10,color = 'black',angle=90,vjust = 1,hjust = 1), # 调整x轴坐标文字
        axis.text.y = element_text(size = 10,color = 'black'), # 调整y轴坐标文字
        legend.text = element_text(size = 10,color = 'black'),
        legend.title = element_text(size = 10,color = 'black'),
        panel.grid.major =  element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = 'right')+
  xlab("")+
  ylab("")
dev.off()



############## Figure3 the overlap gene between top 20 genes in magma and the DEG between EX12 and other Excitatory neurons subtypes
intersect_gene2<-intersect(magma_risk_gene$SYMBOL[1:30],row.names(DEG1))

################## construct a dataframe
dataframe2<-data.frame(gene=magma_risk_gene1$SYMBOL,type=magma_risk_gene1$SYMBOL%in%intersect_gene2)
dataframe2[dataframe2$type=="TRUE",3]<-"1"
dataframe2[dataframe2$type=="FALSE",3]<-"0"
colnames(dataframe2)[3]<-"type1"
dataframe2$Z<-"Z"
dataframe2$gene<-factor(dataframe2$gene,levels = rev(magma_risk_gene1$SYMBOL))


pdf("The_result_of_the_overlap_gene_between_magma_top30_and_DEG_between_EX2_4_and_other_Excitatory_neurons_subtypes.pdf")
ggplot(dataframe2, aes(x=Z,y=gene,fill = type1)) + 
  geom_tile(color="white",linetype = 1)+
  coord_equal()+
  theme_minimal()+
  scale_fill_manual(values = c("#F1F2F1","#45BACD")) + 
  theme(axis.text.x = element_text(size = 10,color = 'black',angle=90,vjust = 1,hjust = 1), # 调整x轴坐标文字
        axis.text.y = element_text(size = 10,color = 'black'), # 调整y轴坐标文字
        legend.text = element_text(size = 10,color = 'black'),
        legend.title = element_text(size = 10,color = 'black'),
        panel.grid.major =  element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = 'right')+
  xlab("")+
  ylab("")
dev.off()



################## Figure4 the overlap gene between top 20 genes in magma and the DEG between TRS top10% and bottom10%
intersect_gene3<-intersect(magma_risk_gene$SYMBOL[1:30],row.names(DEG2))

################## construct a dataframe
dataframe3<-data.frame(gene=magma_risk_gene1$SYMBOL,type=magma_risk_gene1$SYMBOL%in%intersect_gene3)
dataframe3[dataframe3$type=="TRUE",3]<-"1"
dataframe3[dataframe3$type=="FALSE",3]<-"0"
colnames(dataframe3)[3]<-"type1"
dataframe3$Z<-"Z"
dataframe3$gene<-factor(dataframe3$gene,levels = rev(magma_risk_gene1$SYMBOL))


pdf("The_result_of_the_overlap_gene_between_magma_top30_and_DEG_between_TRS_top10_and_bottom10.pdf")
ggplot(dataframe3, aes(x=Z,y=gene,fill = type1)) + 
  geom_tile(color="white",linetype = 1)+
  coord_equal()+
  theme_minimal()+
  scale_fill_manual(values = c("#F1F2F1","#F28E3D")) + 
  theme(axis.text.x = element_text(size = 10,color = 'black',angle=90,vjust = 1,hjust = 1), # 调整x轴坐标文字
        axis.text.y = element_text(size = 10,color = 'black'), # 调整y轴坐标文字
        legend.text = element_text(size = 10,color = 'black'),
        legend.title = element_text(size = 10,color = 'black'),
        panel.grid.major =  element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = 'right')+
  xlab("")+
  ylab("")
dev.off()



################## Figure5 Summarize the results of the figure1-4
dataframe_merge<-cbind(dataframe1[,c("gene","type1")],dataframe2[,c("type1")],dataframe3[,c("type1")])
colnames(dataframe_merge)<-c("gene","sum1","sum2","sum3")
dataframe_merge[,2:4]<-apply(dataframe_merge[,2:4],2,as.numeric)
dataframe_merge$times<-rowSums(dataframe_merge[,2:4])+1
dataframe_merge<-dataframe_merge[order(dataframe_merge$times,decreasing = T),]
dataframe_merge$times<-as.character(dataframe_merge$times)
dataframe_merge$Z<-"z"
dataframe_merge$gene<-factor(dataframe_merge$gene,levels = rev(dataframe_merge$gene))
dataframe_merge$times<-factor(dataframe_merge$times,levels = c("4","3","2","1"))

pdf("Summarize_the_results_of_the_figure1-4.pdf")
ggplot(dataframe_merge, aes(x=Z,y=gene,fill=times)) + 
  geom_tile(color="white",linetype = 1)+
  coord_equal()+
  theme_minimal()+
  scale_fill_manual(values = c("#ED4437","#E68FAC","#EBD57C","#89C8E8")) + 
  theme(axis.text.x = element_text(size = 10,color = 'black',angle=90,vjust = 1,hjust = 1), # 调整x轴坐标文字
        axis.text.y = element_text(size = 10,color = 'black'), # 调整y轴坐标文字
        legend.text = element_text(size = 10,color = 'black'),
        legend.title = element_text(size = 10,color = 'black'),
        panel.grid.major =  element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = 'right')+
  xlab("")+
  ylab("")
dev.off()


















