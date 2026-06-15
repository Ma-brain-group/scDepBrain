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
setwd("D:/Project/SingleCell_MDD/Figure/Figure4/")

##################### In all cells #####################
single_data<-readRDS("D:/Project/SingleCell_MDD/SingleCell_analysis/annotation/MDD_singlecell_data_reannotation_simple.rds")
## get the gene expression
expression<-single_data@assays$RNA@data
expression<-expression[row.names(expression)%in%c("GRM5","NEGR1","PCLO"),]
expression<-as.data.frame(t(as.matrix(expression)))
expression$scPagwas.TRS.Score<-single_data$scPagwas.TRS.Score
## get the pearson correlation and P value
cor.test(expression$NEGR1,expression$scPagwas.TRS.Score)  ### cor 0.7759995  p-value < 2.2e-16
cor.test(expression$PCLO,expression$scPagwas.TRS.Score)  ### cor 0.6902418 p-value < 2.2e-16
cor.test(expression$GRM5,expression$scPagwas.TRS.Score)  ### cor 0.802419 p-value < 2.2e-16

#################### scatter plot
expression$title <- "Correlation between NEGR1 expression and scPagwas TRS Score in all cells"
pdf("correlation_between_NEGR1_expression_and_TRS_in_all_cells.pdf",height = 6)
ggplot(expression,aes(NEGR1,scPagwas.TRS.Score))+
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
  xlab("NEGR1 expression")+
  ylab("scPagwas TRS Score")+
  facet_grid(. ~ title) 
dev.off()


###################### bins boxplot
###################### get bins
### NEGR1
binning <- expression %>% mutate(rank=ntile(expression$NEGR1,5))
xx <- aggregate(binning$NEGR1,by = list(binning$rank), median)
binning$rank<-factor(binning$rank,levels = c("1","2","3","4","5"))

pdf("binning_boxplot_NEGR1_expression_and_TRS_in_all_cells.pdf",height = 4.5,width = 5.5)
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
  xlab("NEGR1 expression gradient quintile bins")+
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

### GRM5
binning <- expression %>% mutate(rank=ntile(expression$GRM5,5))
xx <- aggregate(binning$GRM5,by = list(binning$rank), median)
binning$rank<-factor(binning$rank,levels = c("1","2","3","4","5"))
pdf("binning_boxplot_GRM5_expression_and_TRS_in_all_cells.pdf",height = 4.5,width = 5.5)
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
  xlab("GRM5 expression gradient quintile bins")+
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


### PCLO
binning <- expression %>% mutate(rank=ntile(expression$PCLO,5))
xx <- aggregate(binning$PCLO,by = list(binning$rank), median)
binning$rank<-factor(binning$rank,levels = c("1","2","3","4","5"))
pdf("binning_boxplot_PCLO_expression_and_TRS_in_all_cells.pdf",height = 4.5,width = 5.5)
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
  xlab("PCLO expression gradient quintile bins")+
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



################## In Excitatory neuron ######################
single_data<-readRDS("D:/Project/SingleCell_MDD/SingleCell_analysis/annotation/Excitatory_neuron_subtypes_reannotation.rds")
## get the gene expression
expression<-single_data@assays$RNA@data
expression<-expression[row.names(expression)%in%c("NEGR1","GRM5","PCLO"),]
expression<-as.data.frame(t(as.matrix(expression)))
expression$scPagwas.TRS.Score<-single_data$scPagwas.TRS.Score
## get the pearson correlation and P value
cor.test(expression$NEGR1,expression$scPagwas.TRS.Score)  ### cor 0.8151592  p-value < 2.2e-16
cor.test(expression$GRM5,expression$scPagwas.TRS.Score)  ### cor 0.8052811 p-value < 2.2e-16
cor.test(expression$PCLO,expression$scPagwas.TRS.Score)  ### cor 0.6018046 p-value < 2.2e-16

expression$title <- "Correlation between NEGR1 expression and scPagwas TRS Score in Excitatory neurons"
pdf("correlation_between_NEGR1_expression_and_TRS_in_Exs.pdf",height = 6)
ggplot(expression,aes(NEGR1,scPagwas.TRS.Score))+
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
  xlab("NEGR1 expression")+
  ylab("scPagwas TRS Score")+
  facet_grid(. ~ title)
dev.off()


###################### bins boxplot
###################### get bins
### NEGR1
binning <- expression %>% mutate(rank=ntile(expression$NEGR1,5))
xx <- aggregate(binning$NEGR1,by = list(binning$rank), median)
binning$rank<-factor(binning$rank,levels = c("1","2","3","4","5"))

pdf("binning_boxplot_NEGR1_expression_and_TRS_in_Excitatory_neurons.pdf",height = 4.5,width = 5.5)
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
  xlab("NEGR1 expression gradient quintile bins")+
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

### GRM5
binning <- expression %>% mutate(rank=ntile(expression$GRM5,5))
xx <- aggregate(binning$GRM5,by = list(binning$rank), median)
binning$rank<-factor(binning$rank,levels = c("1","2","3","4","5"))
pdf("binning_boxplot_GRM5_expression_and_TRS_in_Excitatory_neurons.pdf",height = 4.5,width = 5.5)
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
  xlab("GRM5 expression gradient quintile bins")+
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


### PCLO
binning <- expression %>% mutate(rank=ntile(expression$PCLO,5))
xx <- aggregate(binning$PCLO,by = list(binning$rank), median)
binning$rank<-factor(binning$rank,levels = c("1","2","3","4","5"))
pdf("binning_boxplot_PCLO_expression_and_TRS_in_Excitatory_neurons.pdf",height = 4.5,width = 5.5)
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
  xlab("PCLO expression gradient quintile bins")+
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


################### In Ex-L2/4 ############################
single_data<-readRDS("D:/Project/SingleCell_MDD/SingleCell_analysis/Exitatory_neurons_subtype_analysis/Ex_L2_4_NEGR_analysis/NEGR1_EX_L2_4.rds")
## get the gene expression
expression<-single_data@assays$RNA@data
expression<-expression[row.names(expression)%in%c("NEGR1","GRM5","PCLO"),]
expression<-as.data.frame(t(as.matrix(expression)))
expression$scPagwas.TRS.Score<-single_data$scPagwas.TRS.Score
## get the pearson correlation and P value
cor.test(expression$NEGR1,expression$scPagwas.TRS.Score)  ### cor 0.8561865   p-value < 2.2e-16
cor.test(expression$GRM5,expression$scPagwas.TRS.Score)  ### cor 0.8114676  p-value < 2.2e-16
cor.test(expression$PCLO,expression$scPagwas.TRS.Score)  ### cor 0.6151362  p-value < 2.2e-16

#### scatter plot
expression$title <- "Correlation between NEGR1 expression and scPagwas TRS Score in Ex-L2/4"
pdf("correlation_between_NEGR1_expression_and_TRS_in_Ex-L2-4.pdf",height = 6)
ggplot(expression,aes(NEGR1,scPagwas.TRS.Score))+
  geom_point(color="#2BA5B9")+
  geom_smooth(method = 'lm', formula = y ~ x, se = T,color="#D8AD63")+
  stat_cor(data=expression, method = "pearson")+   #相关性检验的R包
  theme_bw()+
  theme(axis.text.x = element_text(size = 10,color = 'black',vjust = 0.5), # 调整x轴坐标文字
        axis.text.y = element_text(size = 10,color = 'black'), # 调整y轴坐标文字
        legend.text = element_text(size = 10,color = 'black'), # 调整legend字体大小
        legend.title = element_text(size = 10,color = 'black'), # 调整legend title大小
        strip.text.x = element_text(size = 10,color= 'black'), ## 调整分面中的字体大小
        panel.grid = element_blank (), # 去除网格线
        plot.margin=unit(c(1,1,1,1),'lines'),#设置边距
        legend.position = 'right')+
  xlab("NEGR1 expression")+
  ylab("scPagwas TRS Score")+
  facet_grid(. ~ title)
dev.off()


###################### bins boxplot
###################### get bins
### NEGR1
binning <- expression %>% mutate(rank=ntile(expression$NEGR1,5))
xx <- aggregate(binning$NEGR1,by = list(binning$rank), median)
binning$rank<-factor(binning$rank,levels = c("1","2","3","4","5"))

pdf("binning_boxplot_NEGR1_expression_and_TRS_in_Ex-L2-4.pdf",height = 4.5,width = 5.5)
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
  xlab("NEGR1 expression gradient quintile bins")+
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

### GRM5
binning <- expression %>% mutate(rank=ntile(expression$GRM5,5))
xx <- aggregate(binning$GRM5,by = list(binning$rank), median)
binning$rank<-factor(binning$rank,levels = c("1","2","3","4","5"))
pdf("binning_boxplot_GRM5_expression_and_TRS_in_Ex-L2-4.pdf",height = 4.5,width = 5.5)
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
  xlab("GRM5 expression gradient quintile bins")+
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


### PCLO
binning <- expression %>% mutate(rank=ntile(expression$PCLO,5))
xx <- aggregate(binning$PCLO,by = list(binning$rank), median)
binning$rank<-factor(binning$rank,levels = c("1","2","3","4","5"))
pdf("binning_boxplot_PCLO_expression_and_TRS_in_Ex-L2-4.pdf",height = 4.5,width = 5.5)
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
  xlab("PCLO expression gradient quintile bins")+
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
