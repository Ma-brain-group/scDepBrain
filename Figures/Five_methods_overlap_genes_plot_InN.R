############## Five methods overlap genes plot #################
library(Seurat)
library(readr)
library(dplyr)
library(stringr)
library(CMplot)
library(clusterProfiler)
library(org.Hs.eg.db)
library(Seurat)
library(ggplot2)
library(reshape2)

setwd("D:/Project/SingleCell_MDD/SingleCell_analysis/Manhattan_Plot/")

############### Figure1 The gene between fusion and magma in magma result
############### import the magma result
magma<-read.table("D:/Project/SingleCell_MDD/GWAS_combine_singlecell_analysis/1_GWAS_analysis_for_broad_celltype/magma/ieu_b_102.genes.out",
                  header = T)
magma$fdr<-p.adjust(magma$P,method = "fdr")

############### change the gene ID from ENTREZID to gene symbol
eg2 <- bitr(magma$GENE,fromType = 'ENTREZID',
            toType = c('SYMBOL'),
            OrgDb='org.Hs.eg.db')
colnames(magma)[1]<-"ENTREZID"
magma<-merge(eg2,magma,by="ENTREZID")
magma1<-magma[magma$fdr<0.05,]
magma1<-magma1[order(magma1$ZSTAT,decreasing = T),]
write.csv(magma1,file = "magma_fdr_0.05.csv",quote = F)


################## import the result of fusion
fusion<-read.csv("D:/Project/SingleCell_MDD/GWAS_analysis/Fusion/fusion_result_Brain_Frontal_Cortex_BA9/Brain_Frontal_Cortex_BA9_result.txt_P_0.05.csv")
intersect_gene<-intersect(fusion$SYMBOL,magma1$SYMBOL)


################## get the common risk gene from magma and fusion result
magma2<-magma1[match(intersect_gene,magma1$SYMBOL),]
magma2<-magma2[,c("SYMBOL","fdr")]
magma2$`-log10(fdr)`<--log10(magma2$fdr)
magma2<-magma2[order(magma2$`-log10(fdr)`,decreasing = T),]
magma2$Z<-"Z"
################## change the character to factor
magma2$SYMBOL<-factor(magma2$SYMBOL,levels = rev(magma2$SYMBOL))


pdf("The_result_of_magma.pdf")
ggplot(magma2, aes(x=Z,y=SYMBOL,fill = -log10(fdr))) + 
  geom_tile(color="white",linetype = 1)+
  coord_equal()+
  theme_minimal()+
  scale_fill_gradient2(low="#F6A45D",midpoint = 1.3,high="#D53D4E") + 
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



################## Figure1 The gene between fusion and magma in fusion result
fusion1<-fusion[match(intersect_gene,fusion$SYMBOL),]
fusion1<-fusion1[,c("SYMBOL","TWAS.P")]
fusion1$`-log10(P)`<--log10(fusion1$TWAS.P)
fusion1<-fusion1[order(fusion1$`-log10(P)`,decreasing = T),]
fusion1$Z<-"Z"
################## change the character to factor
fusion1$SYMBOL<-factor(fusion1$SYMBOL,levels = rev(fusion1$SYMBOL))


pdf("The_result_of_fusion.pdf")
ggplot(fusion1, aes(x=Z,y=SYMBOL,fill = `-log10(P)`)) + 
  geom_tile(color="white",linetype = 1)+
  coord_equal()+
  theme_minimal()+
  scale_fill_gradient2(low="#F8D3C4",midpoint = 1.2,high="#B21F2A") + 
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



################## Figure3 The result of different expressed gene between TRS top10% vs TRS bottom10%
top10_vs_bottom10<-read.csv("D:/Project/SingleCell_MDD/SingleCell_analysis/8_top10_vs_bottom10_In6_marker/findmarker_top10_vs_bottom10_In6.csv",row.names = 1)
intersect_gene1<-intersect(row.names(top10_vs_bottom10),intersect_gene)

################## construct a dataframe
dataframe<-data.frame(gene=intersect_gene,type=intersect_gene%in%intersect_gene1)
dataframe[dataframe$type=="TRUE",3]<-"1"
dataframe[dataframe$type=="FALSE",3]<-"0"
colnames(dataframe)[3]<-"type1"
dataframe$Z<-"Z"
#dataframe$gene<-factor(dataframe$gene,levels = rev(magma_fusion1$SYMBOL))

pdf("The_result_TRS_top10_vs_TRS_bottom10.pdf")
ggplot(dataframe, aes(x=Z,y=gene,fill = type1)) + 
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



################## Figure4 The result of different expressed gene between In6 and other celltypes
In6_vs_other_celltypes_marker<-read.csv("D:/Project/SingleCell_MDD/SingleCell_analysis/7_In6_vs_other_celltype_marker/FindMarkers_In6_vs_other_celltype_gene.csv",row.names = 1)
intersect_gene2<-intersect(row.names(In6_vs_other_celltypes_marker),intersect_gene)

################## construct a dataframe
dataframe1<-data.frame(gene=intersect_gene,type=intersect_gene%in%intersect_gene2)
dataframe1[dataframe1$type=="TRUE",3]<-"1"
dataframe1[dataframe1$type=="FALSE",3]<-"0"
colnames(dataframe1)[3]<-"type1"
dataframe1$Z<-"Z"
#dataframe1$gene<-factor(dataframe1$gene,levels = rev(magma_fusion1$SYMBOL))

pdf("The_result_In6_vs_other_celltypes.pdf")
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




################## Figure5 The result of different expressed gene between In6 and other Ins
In6_vs_other_Ins_marker<-read.csv("D:/Project/SingleCell_MDD/SingleCell_analysis/7_In6_vs_other_celltype_marker/FindMarkers_In6_vs_Ins_gene.csv",row.names = 1)
intersect_gene3<-intersect(row.names(In6_vs_other_Ins_marker),intersect_gene)

################## construct a dataframe
dataframe2<-data.frame(gene=intersect_gene,type=intersect_gene%in%intersect_gene3)
dataframe2[dataframe2$type=="TRUE",3]<-"1"
dataframe2[dataframe2$type=="FALSE",3]<-"0"
colnames(dataframe2)[3]<-"type1"
dataframe2$Z<-"Z"
#dataframe2$gene<-factor(dataframe2$gene,levels = rev(magma_fusion1$SYMBOL))

pdf("The_result_In6_vs_other_Ins.pdf")
ggplot(dataframe2, aes(x=Z,y=gene,fill = type1)) + 
  geom_tile(color="white",linetype = 1)+
  coord_equal()+
  theme_minimal()+
  scale_fill_manual(values = c("grey","#45BACD")) + 
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
dataframe_merge<-cbind(dataframe[,c("gene","type1")],dataframe1[,c("type1")],dataframe2[,c("type1")])
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





