############### CNNM2 is the DEG between M1C and other brain region #################
library(Seurat)
library(dplyr)
library(stringr)
library(ggplot2)
library(scRNAtoolVis)
setwd("D:/Project/SingleCell_MDD/SingleCell_analysis/CNNM2_Brain_region/")

############### In all cells 
single_data_all<-readRDS("D:/Project/SingleCell_MDD/SingleCell_analysis/annotation/MDD_singlecell_data_reannotation_simple.rds")
DimPlot(single_data_all,group.by = "final_anno",reduction = "umap",raster = T,label = T,repel = T)
Idents(single_data_all)<-single_data_all$tissue.region
marker_in_all_cells<-FindMarkers(single_data_all,assay = "RNA",slot = "data",group.by = "tissue.region",ident.1 = "M1C",only.pos = F)
marker_in_all_cells<-marker_in_all_cells[marker_in_all_cells$p_val_adj<0.05,]
write.csv(marker_in_all_cells,file = "brain_region_marker_in_all_cells.csv",quote = F)
marker_in_all_cells<-read.csv("brain_region_marker_in_all_cells.csv",row.names = 1)

############### In all Inhibitory neurons
single_data_Inhibitory_neurons<-readRDS("D:/Project/SingleCell_MDD/SingleCell_analysis/annotation/Inhibitory_neuron_subtypes_reannotation.rds")
DimPlot(single_data_Inhibitory_neurons,group.by = "Inhibitory_neuron_subtypes",reduction = "umap",raster = T,label = T,repel = T)
Idents(single_data_Inhibitory_neurons)<-single_data_Inhibitory_neurons$tissue.region
marker_in_all_Inhibitory_neurons<-FindMarkers(single_data_Inhibitory_neurons,assay = "RNA",slot = "data",group.by = "tissue.region",ident.1 = "M1C",only.pos = F)
marker_in_all_Inhibitory_neurons<-marker_in_all_Inhibitory_neurons[marker_in_all_Inhibitory_neurons$p_val_adj<0.05,]
write.csv(marker_in_all_Inhibitory_neurons,file = "brain_region_marker_in_all_Inhibitory_neurons.csv",quote = F)
marker_in_all_Inhibitory_neurons<-read.csv("brain_region_marker_in_all_Inhibitory_neurons.csv",row.names = 1)

############### In In_PVALB
single_data_In_PVALB<-readRDS("D:/Project/SingleCell_MDD/SingleCell_analysis/11_CNNM2_positive_vs_CNNM2_negnative_IN6/CNNM2_IN6.rds")
#DimPlot(single_data_In_PVALB,group.by = "final_anno",reduction = "umap",raster = T,label = T,repel = T)
Idents(single_data_In_PVALB)<-single_data_In_PVALB$tissue.region
marker_in_In_PVALB<-FindMarkers(single_data_In_PVALB,assay = "RNA",slot = "data",group.by = "tissue.region",ident.1 = "M1C",only.pos = F)
marker_in_In_PVALB<-marker_in_In_PVALB[marker_in_In_PVALB$p_val_adj<0.05,]
write.csv(marker_in_In_PVALB,file = "brain_region_marker_in_In_PVALB.csv",quote = F)
marker_in_In_PVALB<-read.csv("brain_region_marker_in_In_PVALB.csv",row.names = 1)

############### add two another columns, one is gene name, the other is celltypes
marker_in_all_cells$cluster<-"In all cells"
marker_in_all_cells$gene<-row.names(marker_in_all_cells)

marker_in_all_Inhibitory_neurons$cluster<-"In inhibitory neurons"
marker_in_all_Inhibitory_neurons$gene<-row.names(marker_in_all_Inhibitory_neurons)

marker_in_In_PVALB$cluster<-"In In_PVALB"
marker_in_In_PVALB$gene<-row.names(marker_in_In_PVALB)

############### merge all result
marker_merge<-rbind(marker_in_all_cells,marker_in_all_Inhibitory_neurons,marker_in_In_PVALB)
marker_merge$cluster<-factor(marker_merge$cluster,levels = c("In all cells","In inhibitory neurons","In In_PVALB"))
mygene <- c("CNNM2")

pdf("jjVolcano_CNNM2_in_M1C_and_other_brain_region.pdf",height = 5,width = 6)
jjVolcano(diffData = marker_merge,
          tile.col = corrplot::COL2('RdBu', 15)[c(4,8,12)],
          myMarkers = mygene,
          #topGeneN = 5,
          size  = 3,
          fontface = 'italic')
dev.off()

pdf("CNNM2_in_M1C_and_other_brain_region.pdf",height = 5,width = 6)
markerVocalno(markers = marker_merge,
              ownGene = mygene,
              labelCol = ggsci::pal_npg()(9))+
  theme(axis.text.x = element_text (size = 10,color = 'black',angle = 0,hjust = 0.5))
dev.off()


##################### the number of brain region cells in In_PVALB
single_data_In_PVALB<-readRDS("D:/Project/SingleCell_MDD/SingleCell_analysis/11_CNNM2_positive_vs_CNNM2_negnative_IN6/CNNM2_IN6.rds")
table(single_data_In_PVALB$tissue.region)
Idents(single_data_In_PVALB)<-single_data_In_PVALB$tissue.region
cell_number<-as.data.frame(table(single_data_In_PVALB$tissue.region))
cell_number$brain_region_colors<-c("#ED4437","#E1884A","#8ACC72","#1F78B4","#89C8E8","#B3446C","#EBD57C","#E68FAC","#CAA2F4",
                                   "#96873B","#B49D99","#B37557","#FC9A9A","#6A3D9A")
cell_number<-cell_number[order(cell_number$Freq,decreasing = T),]
cell_number$Var1<-factor(cell_number$Var1,levels = cell_number$Var1)

pdf("The_number_of_In_PVALB_in_different_brain_region.pdf",height = 4.5,width = 5)
ggplot(cell_number,aes(Var1,Freq))+
  geom_col(aes(fill=Var1),width = 0.6)+
  scale_fill_manual(values = cell_number$brain_region_colors)+  ##指定颜色
  geom_text(aes(label=Freq),size = 3,vjust=-0.5) + 
  theme_bw()+
  xlab(" ")+
  ylab("Cell counts")+
  theme(axis.text.x = element_text(size = 10,color = 'black',vjust = 1,hjust = 1,angle = 35), # 调整x轴坐标文字
        axis.text.y = element_text(size = 10,color = 'black'), # 调整y轴坐标文字
        legend.text = element_text(size = 10,color = 'black'), # 调整legend字体大小
        legend.title = element_text(size = 10,color = 'black'), # 调整legend title大小
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
  scale_fill_manual(values = cell_number$brain_region_colors)+  ##指定颜色
  geom_text(aes(label=Freq),size = 3,vjust=-0.5) + 
  theme_bw()+
  xlab(" ")+
  ylab("Cell counts")+
  theme(axis.text.x = element_text(size = 10,color = 'black',vjust = 1,hjust = 1,angle = 35), # 调整x轴坐标文字
        axis.text.y = element_text(size = 10,color = 'black'), # 调整y轴坐标文字
        legend.text = element_text(size = 10,color = 'black'), # 调整legend字体大小
        legend.title = element_text(size = 10,color = 'black'), # 调整legend title大小
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
dev.off()


################ calculate the percent of In_PVALB in different brain region
cell_percent<-data.frame(brain_region=names(table(single_data_In_PVALB$tissue.region)/table(single_data_all$tissue.region)),
                         percent=as.numeric(table(single_data_In_PVALB$tissue.region)/table(single_data_all$tissue.region))*100)
cell_percent$brain_region_colors<-c("#ED4437","#E1884A","#8ACC72","#1F78B4","#89C8E8","#B3446C","#EBD57C","#E68FAC","#CAA2F4",
                                   "#96873B","#B49D99","#B37557","#FC9A9A","#6A3D9A")
cell_percent<-cell_percent[order(cell_percent$percent,decreasing = T),]
cell_percent$brain_region<-factor(cell_percent$brain_region,levels = cell_percent$brain_region)

pdf("The_percent_of_In_PVALB_in_different_brain_region.pdf",height = 4.5,width = 5)
ggplot(cell_percent,aes(brain_region,percent))+
  geom_col(aes(fill=brain_region),width = 0.6)+
  scale_fill_manual(values = cell_percent$brain_region_colors)+  ##指定颜色
  geom_text(aes(label=round(percent, 2)),size = 3,vjust=-0.5) + 
  theme_bw()+
  xlab(" ")+
  ylab("Cell percent (%)")+
  theme(axis.text.x = element_text(size = 10,color = 'black',vjust = 1,hjust = 1,angle = 35), # 调整x轴坐标文字
        axis.text.y = element_text(size = 10,color = 'black'), # 调整y轴坐标文字
        legend.text = element_text(size = 10,color = 'black'), # 调整legend字体大小
        legend.title = element_text(size = 10,color = 'black'), # 调整legend title大小
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
dev.off()














