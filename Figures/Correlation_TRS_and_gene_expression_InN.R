##################### correlation between gene expression and TRS score ########################
rm(list = ls())
library(Seurat)
library(Nebulosa)
library(BiocFileCache)
library(paletteer)
library(scCustomize)
library(ggplot2)
library(dplyr)
library(cowplot) 
library(ggpubr)
setwd("D:/Project/SingleCell_MDD/Figure/Figure3/")

##################### In all cells #####################
single_data<-readRDS("D:/Project/SingleCell_MDD/SingleCell_analysis/annotation/MDD_singlecell_data_reannotation_simple.rds")
## get the gene expression
expression<-single_data@assays$RNA@data
expression<-expression[row.names(expression)%in%c("CNNM2","ESYT2","MAP3K7"),]
expression<-as.data.frame(t(as.matrix(expression)))
expression$scPagwas.TRS.Score<-single_data$scPagwas.TRS.Score
## get the pearson correlation and P value
cor.test(expression$CNNM2,expression$scPagwas.TRS.Score)  ### cor 0.3700493  p-value < 2.2e-16
cor.test(expression$ESYT2,expression$scPagwas.TRS.Score)  ### cor 0.444034 p-value < 2.2e-16
cor.test(expression$MAP3K7,expression$scPagwas.TRS.Score)  ### cor 0.2013972 p-value < 2.2e-16

#################### scatter plot
expression$title <- "Correlation between CNNM2 expression and scPagwas TRS Score in all cells"
pdf("correlation_between_CNNM2_expression_and_TRS_in_all_cells.pdf",height = 6)
ggplot(expression,aes(CNNM2,scPagwas.TRS.Score))+
  geom_point(color="#2BA5B9")+
  geom_smooth(method = 'lm', formula = y ~ x, se = T,color="#D8AD63")+
  stat_cor(data=expression, method = "pearson")+  
  theme_bw()+
  theme(axis.text.x = element_text(size = 10,color = 'black',vjust = 0.5), 
        axis.text.y = element_text(size = 10,color = 'black'), 
        legend.text = element_text(size = 10,color = 'black'), 
        legend.title = element_text(size = 10,color = 'black'), 
        strip.text.x = element_text(size = 10,color= 'black'), 
        panel.grid = element_blank (), 
        plot.margin=unit(c(1,1,1,1),'lines'),
        legend.position = 'right')+
  xlab("CNNM2 expression")+
  ylab("scPagwas TRS Score")+
  facet_grid(. ~ title) 
dev.off()


###################### bins boxplot
###################### get bins
### CNNM2
binning <- expression %>% mutate(rank=ntile(expression$CNNM2,5))
xx <- aggregate(binning$CNNM2,by = list(binning$rank), median)
binning$rank<-factor(binning$rank,levels = c("1","2","3","4","5"))

pdf("binning_boxplot_CNNM2_expression_and_TRS_in_all_cells.pdf",height = 4.5,width = 5.5)
ggplot(binning, aes(x=rank, y=scPagwas.TRS.Score)) + 
  stat_boxplot(geom = "errorbar",width=0.05, size=0.5,position=position_dodge(0.6),color= "black")+
  theme_bw()+
  geom_boxplot(position = position_dodge(0.6),
               size = 0.5,
               width = 0.7,
               fill = c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3'),
               color = "black",
               outlier.color = "black",
               outlier.fill = "black",
               outlier.shape = 19,
               outlier.size = 1,
               outlier.stroke = 0.5,
               outlier.alpha = 45,
               notch = F,
               notchwidth = 0.5)+
  xlab("CNNM2 expression gradient quintile bins")+
  ylab("scPagwas.TRS.Score")+
  theme(axis.text.x = element_text(size = 10,color = 'black',vjust = 1,hjust = 0.5,angle = 0),
        axis.text.y = element_text(size = 10,color = 'black'), 
        legend.text = element_text(size = 10,color = 'black'), 
        legend.title = element_text(size = 10,color = 'black'), 
        panel.grid.major =  element_blank(),
        panel.grid.minor = element_blank())+
  annotate("text", x = Inf, y = Inf, label = "P-value < 2.2e-16", 
           hjust = 3.7, vjust = 2, size = 4, color = "black")
dev.off()

### ESYT2
binning <- expression %>% mutate(rank=ntile(expression$ESYT2,5))
xx <- aggregate(binning$ESYT2,by = list(binning$rank), median)
binning$rank<-factor(binning$rank,levels = c("1","2","3","4","5"))
pdf("binning_boxplot_ESYT2_expression_and_TRS_in_all_cells.pdf",height = 4.5,width = 5.5)
ggplot(binning, aes(x=rank, y=scPagwas.TRS.Score)) + 
  stat_boxplot(geom = "errorbar",width=0.05, size=0.5,position=position_dodge(0.6),color= "black")+
  theme_bw()+
  geom_boxplot(position = position_dodge(0.6),
               size = 0.5,
               width = 0.7,
               fill = c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3'),
               color = "black",
               outlier.color = "black",
               outlier.fill = "black",
               outlier.shape = 19,
               outlier.size = 1,
               outlier.stroke = 0.5,
               outlier.alpha = 45,
               notch = F,
               notchwidth = 0.5)+
  xlab("ESYT2 expression gradient quintile bins")+
  ylab("scPagwas.TRS.Score")+
  theme(axis.text.x = element_text(size = 10,color = 'black',vjust = 1,hjust = 0.5,angle = 0),
        axis.text.y = element_text(size = 10,color = 'black'), 
        legend.text = element_text(size = 10,color = 'black'), 
        legend.title = element_text(size = 10,color = 'black'), 
        panel.grid.major =  element_blank(),
        panel.grid.minor = element_blank())+
  annotate("text", x = Inf, y = Inf, label = "P-value < 2.2e-16", 
           hjust = 3.7, vjust = 2, size = 4, color = "black")
dev.off()


### MAP3K7
binning <- expression %>% mutate(rank=ntile(expression$MAP3K7,5))
xx <- aggregate(binning$MAP3K7,by = list(binning$rank), median)
binning$rank<-factor(binning$rank,levels = c("1","2","3","4","5"))
pdf("binning_boxplot_MAP3K7_expression_and_TRS_in_all_cells.pdf",height = 4.5,width = 5.5)
ggplot(binning, aes(x=rank, y=scPagwas.TRS.Score)) + 
  stat_boxplot(geom = "errorbar",width=0.05, size=0.5,position=position_dodge(0.6),color= "black")+
  theme_bw()+
  geom_boxplot(position = position_dodge(0.6),
               size = 0.5,
               width = 0.7,
               fill = c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3'),
               color = "black",
               outlier.color = "black",
               outlier.fill = "black",
               outlier.shape = 19,
               outlier.size = 1,
               outlier.stroke = 0.5,
               outlier.alpha = 45,
               notch = F,
               notchwidth = 0.5)+
  xlab("MAP3K7 expression gradient quintile bins")+
  ylab("scPagwas.TRS.Score")+
  theme(axis.text.x = element_text(size = 10,color = 'black',vjust = 1,hjust = 0.5,angle = 0),
        axis.text.y = element_text(size = 10,color = 'black'), 
        legend.text = element_text(size = 10,color = 'black'), 
        legend.title = element_text(size = 10,color = 'black'), 
        panel.grid.major =  element_blank(),
        panel.grid.minor = element_blank())+
  annotate("text", x = Inf, y = Inf, label = "P-value < 2.2e-16", 
           hjust = 3.7, vjust = 2, size = 4, color = "black")
dev.off()



################## In inhibitory neurons ######################
single_data<-readRDS("D:/Project/SingleCell_MDD/SingleCell_analysis/annotation/Inhibitory_neuron_subtypes_reannotation.rds")
## get the gene expression
expression<-single_data@assays$RNA@data
expression<-expression[row.names(expression)%in%c("CNNM2","ESYT2","MAP3K7"),]
expression<-as.data.frame(t(as.matrix(expression)))
expression$scPagwas.TRS.Score<-single_data$scPagwas.TRS.Score
## get the pearson correlation and P value
cor.test(expression$CNNM2,expression$scPagwas.TRS.Score)  ### cor 0.4701042  p-value < 2.2e-16
cor.test(expression$ESYT2,expression$scPagwas.TRS.Score)  ### cor 0.4059577 p-value < 2.2e-16
cor.test(expression$MAP3K7,expression$scPagwas.TRS.Score)  ### cor 0.2321215 p-value < 2.2e-16

expression$title <- "Correlation between CNNM2 expression and scPagwas TRS Score in Inhibitory neurons"
pdf("correlation_between_CNNM2_expression_and_TRS_in_INs.pdf",height = 6)
ggplot(expression,aes(CNNM2,scPagwas.TRS.Score))+
  geom_point(color="#2BA5B9")+
  geom_smooth(method = 'lm', formula = y ~ x, se = T,color="#D8AD63")+
  stat_cor(data=expression, method = "pearson")+   
  theme_bw()+
  theme(axis.text.x = element_text(size = 10,color = 'black',vjust = 0.5), 
        axis.text.y = element_text(size = 10,color = 'black'), 
        legend.text = element_text(size = 10,color = 'black'),
        legend.title = element_text(size = 10,color = 'black'),
        strip.text.x = element_text(size = 10,color= 'black'), 
        panel.grid = element_blank (), 
        plot.margin=unit(c(1,1,1,1),'lines'),
        legend.position = 'right')+
  xlab("CNNM2 expression")+
  ylab("scPagwas TRS Score")+
  facet_grid(. ~ title)
dev.off()


###################### bins boxplot
###################### get bins
### CNNM2
binning <- expression %>% mutate(rank=ntile(expression$CNNM2,5))
xx <- aggregate(binning$CNNM2,by = list(binning$rank), median)
binning$rank<-factor(binning$rank,levels = c("1","2","3","4","5"))

pdf("binning_boxplot_CNNM2_expression_and_TRS_in_Inhibitory_neurons.pdf",height = 4.5,width = 5.5)
ggplot(binning, aes(x=rank, y=scPagwas.TRS.Score)) + 
  stat_boxplot(geom = "errorbar",width=0.05, size=0.5,position=position_dodge(0.6),color= "black")+
  theme_bw()+
  geom_boxplot(position = position_dodge(0.6),
               size = 0.5,
               width = 0.7,
               fill = c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3'),
               color = "black",
               outlier.color = "black",
               outlier.fill = "black",
               outlier.shape = 19,
               outlier.size = 1,
               outlier.stroke = 0.5,
               outlier.alpha = 45,
               notch = F,
               notchwidth = 0.5)+
  xlab("CNNM2 expression gradient quintile bins")+
  ylab("scPagwas.TRS.Score")+
  theme(axis.text.x = element_text(size = 10,color = 'black',vjust = 1,hjust = 0.5,angle = 0),
        axis.text.y = element_text(size = 10,color = 'black'), 
        legend.text = element_text(size = 10,color = 'black'), 
        legend.title = element_text(size = 10,color = 'black'), 
        panel.grid.major =  element_blank(),
        panel.grid.minor = element_blank())+
  annotate("text", x = Inf, y = Inf, label = "P-value < 2.2e-16", 
           hjust = 3.7, vjust = 2, size = 4, color = "black")
dev.off()

### ESYT2
binning <- expression %>% mutate(rank=ntile(expression$ESYT2,5))
xx <- aggregate(binning$ESYT2,by = list(binning$rank), median)
binning$rank<-factor(binning$rank,levels = c("1","2","3","4","5"))
pdf("binning_boxplot_ESYT2_expression_and_TRS_in_Inhibitory_neurons.pdf",height = 4.5,width = 5.5)
ggplot(binning, aes(x=rank, y=scPagwas.TRS.Score)) + 
  stat_boxplot(geom = "errorbar",width=0.05, size=0.5,position=position_dodge(0.6),color= "black")+
  theme_bw()+
  geom_boxplot(position = position_dodge(0.6),
               size = 0.5,
               width = 0.7,
               fill = c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3'),
               color = "black",
               outlier.color = "black",
               outlier.fill = "black",
               outlier.shape = 19,
               outlier.size = 1,
               outlier.stroke = 0.5,
               outlier.alpha = 45,
               notch = F,
               notchwidth = 0.5)+
  xlab("ESYT2 expression gradient quintile bins")+
  ylab("scPagwas.TRS.Score")+
  theme(axis.text.x = element_text(size = 10,color = 'black',vjust = 1,hjust = 0.5,angle = 0),
        axis.text.y = element_text(size = 10,color = 'black'), 
        legend.text = element_text(size = 10,color = 'black'), 
        legend.title = element_text(size = 10,color = 'black'), 
        panel.grid.major =  element_blank(),
        panel.grid.minor = element_blank())+
  annotate("text", x = Inf, y = Inf, label = "P-value < 2.2e-16", 
           hjust = 3.7, vjust = 2, size = 4, color = "black")
dev.off()


### MAP3K7
binning <- expression %>% mutate(rank=ntile(expression$MAP3K7,5))
xx <- aggregate(binning$MAP3K7,by = list(binning$rank), median)
binning$rank<-factor(binning$rank,levels = c("1","2","3","4","5"))
pdf("binning_boxplot_MAP3K7_expression_and_TRS_in_Inhibitory_neurons.pdf",height = 4.5,width = 5.5)
ggplot(binning, aes(x=rank, y=scPagwas.TRS.Score)) + 
  stat_boxplot(geom = "errorbar",width=0.05, size=0.5,position=position_dodge(0.6),color= "black")+
  theme_bw()+
  geom_boxplot(position = position_dodge(0.6),
               size = 0.5,
               width = 0.7,
               fill = c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3'),
               color = "black",
               outlier.color = "black",
               outlier.fill = "black",
               outlier.shape = 19,
               outlier.size = 1,
               outlier.stroke = 0.5,
               outlier.alpha = 45,
               notch = F,
               notchwidth = 0.5)+
  xlab("MAP3K7 expression gradient quintile bins")+
  ylab("scPagwas.TRS.Score")+
  theme(axis.text.x = element_text(size = 10,color = 'black',vjust = 1,hjust = 0.5,angle = 0),
        axis.text.y = element_text(size = 10,color = 'black'), 
        legend.text = element_text(size = 10,color = 'black'), 
        legend.title = element_text(size = 10,color = 'black'), 
        panel.grid.major =  element_blank(),
        panel.grid.minor = element_blank())+
  annotate("text", x = Inf, y = Inf, label = "P-value < 2.2e-16", 
           hjust = 3.7, vjust = 2, size = 4, color = "black")
dev.off()


################### In In_PVALB ############################
single_data<-readRDS("D:/Project/SingleCell_MDD/SingleCell_analysis/11_CNNM2_positive_vs_CNNM2_negnative_IN6/CNNM2_IN6.rds")
## get the gene expression
expression<-single_data@assays$RNA@data
expression<-expression[row.names(expression)%in%c("CNNM2","ESYT2","MAP3K7"),]
expression<-as.data.frame(t(as.matrix(expression)))
expression$scPagwas.TRS.Score<-single_data$scPagwas.TRS.Score
## get the pearson correlation and P value
cor.test(expression$CNNM2,expression$scPagwas.TRS.Score)  ### cor 0.4084575   p-value < 2.2e-16
cor.test(expression$ESYT2,expression$scPagwas.TRS.Score)  ### cor 0.3712518  p-value < 2.2e-16
cor.test(expression$MAP3K7,expression$scPagwas.TRS.Score)  ### cor 0.2426953  p-value < 2.2e-16

#### scatter plot
expression$title <- "Correlation between CNNM2 expression and scPagwas TRS Score in In_PVALB"
pdf("correlation_between_CNNM2_expression_and_TRS_in_IN6.pdf",height = 6)
ggplot(expression,aes(CNNM2,scPagwas.TRS.Score))+
  geom_point(color="#2BA5B9")+
  geom_smooth(method = 'lm', formula = y ~ x, se = T,color="#D8AD63")+
  stat_cor(data=expression, method = "pearson")+   
  theme_bw()+
  theme(axis.text.x = element_text(size = 10,color = 'black',vjust = 0.5), 
        axis.text.y = element_text(size = 10,color = 'black'), 
        legend.text = element_text(size = 10,color = 'black'), 
        legend.title = element_text(size = 10,color = 'black'), 
        strip.text.x = element_text(size = 10,color= 'black'), 
        panel.grid = element_blank (), 
        plot.margin=unit(c(1,1,1,1),'lines'),
        legend.position = 'right')+
  xlab("CNNM2 expression")+
  ylab("scPagwas TRS Score")+
  facet_grid(. ~ title)
dev.off()


###################### bins boxplot
###################### get bins
### CNNM2
binning <- expression %>% mutate(rank=ntile(expression$CNNM2,5))
xx <- aggregate(binning$CNNM2,by = list(binning$rank), median)
binning$rank<-factor(binning$rank,levels = c("1","2","3","4","5"))

pdf("binning_boxplot_CNNM2_expression_and_TRS_in_In_PVALB.pdf",height = 4.5,width = 5.5)
ggplot(binning, aes(x=rank, y=scPagwas.TRS.Score)) + 
  stat_boxplot(geom = "errorbar",width=0.05, size=0.5,position=position_dodge(0.6),color= "black")+
  theme_bw()+
  geom_boxplot(position = position_dodge(0.6),
               size = 0.5,
               width = 0.7,
               fill = c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3'),
               color = "black",
               outlier.color = "black",
               outlier.fill = "black",
               outlier.shape = 19,
               outlier.size = 1,
               outlier.stroke = 0.5,
               outlier.alpha = 45,
               notch = F,
               notchwidth = 0.5)+
  xlab("CNNM2 expression gradient quintile bins")+
  ylab("scPagwas.TRS.Score")+
  theme(axis.text.x = element_text(size = 10,color = 'black',vjust = 1,hjust = 0.5,angle = 0),
        axis.text.y = element_text(size = 10,color = 'black'), 
        legend.text = element_text(size = 10,color = 'black'), 
        legend.title = element_text(size = 10,color = 'black'), 
        panel.grid.major =  element_blank(),
        panel.grid.minor = element_blank())+
  annotate("text", x = Inf, y = Inf, label = "P-value < 2.2e-16", 
           hjust = 3.7, vjust = 2, size = 4, color = "black")
dev.off()

### ESYT2
binning <- expression %>% mutate(rank=ntile(expression$ESYT2,5))
xx <- aggregate(binning$ESYT2,by = list(binning$rank), median)
binning$rank<-factor(binning$rank,levels = c("1","2","3","4","5"))
pdf("binning_boxplot_ESYT2_expression_and_TRS_in_In_PVALB.pdf",height = 4.5,width = 5.5)
ggplot(binning, aes(x=rank, y=scPagwas.TRS.Score)) + 
  stat_boxplot(geom = "errorbar",width=0.05, size=0.5,position=position_dodge(0.6),color= "black")+
  theme_bw()+
  geom_boxplot(position = position_dodge(0.6),
               size = 0.5,
               width = 0.7,
               fill = c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3'),
               color = "black",
               outlier.color = "black",
               outlier.fill = "black",
               outlier.shape = 19,
               outlier.size = 1,
               outlier.stroke = 0.5,
               outlier.alpha = 45,
               notch = F,
               notchwidth = 0.5)+
  xlab("ESYT2 expression gradient quintile bins")+
  ylab("scPagwas.TRS.Score")+
  theme(axis.text.x = element_text(size = 10,color = 'black',vjust = 1,hjust = 0.5,angle = 0),
        axis.text.y = element_text(size = 10,color = 'black'), 
        legend.text = element_text(size = 10,color = 'black'), 
        legend.title = element_text(size = 10,color = 'black'), 
        panel.grid.major =  element_blank(),
        panel.grid.minor = element_blank())+
  annotate("text", x = Inf, y = Inf, label = "P-value < 2.2e-16", 
           hjust = 3.7, vjust = 2, size = 4, color = "black")
dev.off()


### MAP3K7
binning <- expression %>% mutate(rank=ntile(expression$MAP3K7,5))
xx <- aggregate(binning$MAP3K7,by = list(binning$rank), median)
binning$rank<-factor(binning$rank,levels = c("1","2","3","4","5"))
pdf("binning_boxplot_MAP3K7_expression_and_TRS_in_In_PVALB.pdf",height = 4.5,width = 5.5)
ggplot(binning, aes(x=rank, y=scPagwas.TRS.Score)) + 
  stat_boxplot(geom = "errorbar",width=0.05, size=0.5,position=position_dodge(0.6),color= "black")+
  theme_bw()+
  geom_boxplot(position = position_dodge(0.6),
               size = 0.5,
               width = 0.7,
               fill = c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3'),
               color = "black",
               outlier.color = "black",
               outlier.fill = "black",
               outlier.shape = 19,
               outlier.size = 1,
               outlier.stroke = 0.5,
               outlier.alpha = 45,
               notch = F,
               notchwidth = 0.5)+
  xlab("MAP3K7 expression gradient quintile bins")+
  ylab("scPagwas.TRS.Score")+
  theme(axis.text.x = element_text(size = 10,color = 'black',vjust = 1,hjust = 0.5,angle = 0),
        axis.text.y = element_text(size = 10,color = 'black'), 
        legend.text = element_text(size = 10,color = 'black'), 
        legend.title = element_text(size = 10,color = 'black'), 
        panel.grid.major =  element_blank(),
        panel.grid.minor = element_blank())+
  annotate("text", x = Inf, y = Inf, label = "P-value < 2.2e-16", 
           hjust = 3.7, vjust = 2, size = 4, color = "black")
dev.off()




