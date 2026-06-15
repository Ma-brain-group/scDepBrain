############### The correlation between scPagwas and scDRS ###################
library(Seurat)
library(ggplot2)
library(cowplot)
library(ggrepel)
library(ggpubr)
setwd("D:/Project/SingleCell_MDD/GWAS_combine_singlecell_analysis/1_GWAS_analysis_for_broad_celltype/scDRS/")

############### import the singlecell data
single_data<-readRDS("D:/Project/SingleCell_MDD/SingleCell_analysis/annotation/MDD_singlecell_data_reannotation_simple.rds")

############### get the scPagwas score 
TRS<-FetchData(single_data,vars = c("anno","scPagwas.TRS.Score"))
TRS1<-aggregate(TRS$scPagwas.TRS.Score,by=list(TRS$anno),FUN=median)
colnames(TRS1)<-c("cell_type","TRS")

############### get the scDRS score
scRDS_score = read.table(gzfile("MDD.full_score.gz"), header = T)
#write.csv(scRDS_score[,1:6],file = "scDRS_singlecell_result.csv",quote = F)
single_data$scRDS_score<-scRDS_score$norm_score
scDRS<-FetchData(single_data,vars = c("anno","scRDS_score"))
scDRS1<-aggregate(scDRS$scRDS_score,by=list(scDRS$anno),FUN=median)
colnames(scDRS1)<-c("cell_type","scDRS")

###############  merge the data
temp<-cbind(TRS1,scDRS1)
temp[,3]<-NULL
cor.test(temp$TRS,temp$scDRS) 

############### plot the correlation of celltype
pdf("Celltypes_correlation_between_scPagwas_TRS_score_and_scDRS_score.pdf",height = 5,width = 6)
ggplot(temp,aes(TRS,scDRS))+
  geom_point(color="#2BA5B9",size=3.5)+
  geom_smooth(method = 'lm', formula = y ~ x, se = T,color="#D8AD63")+
  stat_cor(data=temp, method = "pearson")+
  theme_bw()+
  xlab("scPagwas TRS score")+
  ylab("scDRS score")+
  geom_text_repel(aes(label = cell_type),size = 3.5)+
  theme(axis.text.x = element_text(size = 10,color = 'black',vjust = 0.5), 
    axis.text.y = element_text(size = 10,color = 'black'), 
    legend.text = element_text(size = 10,color = 'black'), 
    legend.title = element_text(size = 10,color = 'black'), 
    panel.grid = element_blank(),
    plot.title = element_text (hjust = 0.5),
    plot.margin=unit(c(1,1,1,1),'lines'),
    legend.position = 'right')
dev.off()


################ plot the correlation of singlecell
temp1<-FetchData(single_data,vars = c("scPagwas.TRS.Score","scRDS_score"))
cor.test(temp1$scPagwas.TRS.Score,temp1$scRDS_score)

pdf("Singlecell_correlation_between_scPagwas_TRS_score_and_scDRS_score.pdf",height = 5,width = 6)
ggplot(temp1,aes(scPagwas.TRS.Score,scRDS_score))+
  geom_point(color="#2BA5B9",size=1)+
  geom_smooth(method = 'lm', formula = y ~ x, se = T,color="#D8AD63")+
  stat_cor(data=temp1, method = "pearson")+
  theme_bw()+
  xlab("scPagwas TRS score")+
  ylab("scDRS score")+
  theme(axis.text.x = element_text(size = 10,color = 'black',vjust = 0.5), 
        axis.text.y = element_text(size = 10,color = 'black'), 
        legend.text = element_text(size = 10,color = 'black'), 
        legend.title = element_text(size = 10,color = 'black'), 
        panel.grid = element_blank(),
        plot.title = element_text (hjust = 0.5),
        plot.margin=unit(c(1,1,1,1),'lines'),
        legend.position = 'right')
dev.off()


