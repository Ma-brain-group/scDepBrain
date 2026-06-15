##################### calculate the number of Ex-L2/4 in different brain region #################
library(Seurat)
library(dplyr)
library(stringr)
library(ggplot2)
library(scRNAtoolVis)
setwd("D:/Project/SingleCell_MDD/SingleCell_analysis/Exitatory_neurons_subtype_analysis/Ex_L2_4_brain_region/")

##################### the number of brain region cells in Ex-L2/4
single_data_Ex_L2_4<-readRDS("D:/Project/SingleCell_MDD/SingleCell_analysis/Exitatory_neurons_subtype_analysis/Ex_L2_4_NEGR_analysis/NEGR1_EX_L2_4.rds")
table(single_data_Ex_L2_4$tissue.region)
Idents(single_data_Ex_L2_4)<-single_data_Ex_L2_4$tissue.region
cell_number<-as.data.frame(table(single_data_Ex_L2_4$tissue.region))
cell_number$brain_region_colors<-c("#ED4437","#E1884A","#8ACC72","#1F78B4","#89C8E8","#B3446C","#EBD57C","#E68FAC","#CAA2F4",
                                   "#96873B","#B49D99","#B37557","#FC9A9A","#6A3D9A")
cell_number<-cell_number[order(cell_number$Freq,decreasing = T),]
cell_number$Var1<-factor(cell_number$Var1,levels = cell_number$Var1)

pdf("The_number_of_Ex_L2_4_in_different_brain_region.pdf",height = 4.5,width = 5)
ggplot(cell_number,aes(Var1,Freq))+
  geom_col(aes(fill=Var1),width = 0.6)+
  scale_fill_manual(values = cell_number$brain_region_colors)+  
  geom_text(aes(label=Freq),size = 3,vjust=-0.5) + 
  theme_bw()+
  xlab(" ")+
  ylab("Cell counts")+
  theme(axis.text.x = element_text(size = 10,color = 'black',vjust = 1,hjust = 1,angle = 35), 
        axis.text.y = element_text(size = 10,color = 'black'), 
        legend.text = element_text(size = 10,color = 'black'), 
        legend.title = element_text(size = 10,color = 'black'), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
dev.off()


##################### the number of brain region cells in 30w cells
single_data_all<-readRDS("D:/Project/SingleCell_MDD/SingleCell_analysis/annotation/MDD_singlecell_data_reannotation_simple.rds")
table(single_data_all$tissue.region)
Idents(single_data_all)<-single_data_all$tissue.region
cell_number<-as.data.frame(table(single_data_all$tissue.region))
cell_number$brain_region_colors<-c("#ED4437","#E1884A","#8ACC72","#1F78B4","#89C8E8","#B3446C","#EBD57C","#E68FAC","#CAA2F4",
                                   "#96873B","#B49D99","#B37557","#FC9A9A","#6A3D9A")
cell_number<-cell_number[order(cell_number$Freq,decreasing = T),]
cell_number$Var1<-factor(cell_number$Var1,levels = cell_number$Var1)

pdf("The_number_of_all_cells_in_different_brain_region.pdf",height = 4.5,width = 5)
ggplot(cell_number,aes(Var1,Freq))+
  geom_col(aes(fill=Var1),width = 0.6)+
  scale_fill_manual(values = cell_number$brain_region_colors)+  
  geom_text(aes(label=Freq),size = 3,vjust=-0.5) + 
  theme_bw()+
  xlab(" ")+
  ylab("Cell counts")+
  theme(axis.text.x = element_text(size = 10,color = 'black',vjust = 1,hjust = 1,angle = 35), 
        axis.text.y = element_text(size = 10,color = 'black'), 
        legend.text = element_text(size = 10,color = 'black'), 
        legend.title = element_text(size = 10,color = 'black'), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
dev.off()


################ calculate the percent of In_PVALB in different brain region
cell_percent<-data.frame(brain_region=names(table(single_data_Ex_L2_4$tissue.region)/table(single_data_all$tissue.region)),
                         percent=as.numeric(table(single_data_Ex_L2_4$tissue.region)/table(single_data_all$tissue.region))*100)
cell_percent$brain_region_colors<-c("#ED4437","#E1884A","#8ACC72","#1F78B4","#89C8E8","#B3446C","#EBD57C","#E68FAC","#CAA2F4",
                                    "#96873B","#B49D99","#B37557","#FC9A9A","#6A3D9A")
cell_percent<-cell_percent[order(cell_percent$percent,decreasing = T),]
cell_percent$brain_region<-factor(cell_percent$brain_region,levels = cell_percent$brain_region)

pdf("The_percent_of_Ex_L2_4_in_different_brain_region.pdf",height = 4.5,width = 5)
ggplot(cell_percent,aes(brain_region,percent))+
  geom_col(aes(fill=brain_region),width = 0.6)+
  scale_fill_manual(values = cell_percent$brain_region_colors)+  
  geom_text(aes(label=round(percent, 2)),size = 3,vjust=-0.5) + 
  theme_bw()+
  xlab(" ")+
  ylab("Cell percent (%)")+
  theme(axis.text.x = element_text(size = 10,color = 'black',vjust = 1,hjust = 1,angle = 35), 
        axis.text.y = element_text(size = 10,color = 'black'), 
        legend.text = element_text(size = 10,color = 'black'), 
        legend.title = element_text(size = 10,color = 'black'), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
dev.off()
