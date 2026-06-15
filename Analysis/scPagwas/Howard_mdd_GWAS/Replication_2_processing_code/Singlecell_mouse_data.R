##################### mouse singlecell plot ####################
library(Seurat)
library(ggplot2)
setwd("D:/Project/SingleCell_MDD/Figure/Singlecell_Mouse/annotation/")

##################### import singlecell data
single_data<-readRDS("D:/Project/SingleCell_MDD/Mouse_data/merge_result/mouse_merge1.rds")
cluster_colors<-c("#4ea64a","#8e4c99","#e47faf","#b698c5","#d5231d","#e88f18","#3777ac")

##################### DimPlot for cell_type
pdf("DimPlot_for_cell_type.pdf",height = 4.5,width = 5.5)
DimPlot(single_data,group.by ="cell_type",reduction = "umap",label = T,repel=T,raster=FALSE,cols = cluster_colors)
dev.off()

##################### DimPlot for datasets batch
Idents(single_data)<-single_data$datasets
single_data$datasets[WhichCells(single_data,idents = "10X")]<-"BICCN2020"
pdf("DimPlot_for_datasets_batch.pdf",height = 4.5,width = 5.5)
DimPlot(single_data,group.by ="datasets",reduction = "umap",label = F,raster=FALSE,cols = cluster_colors)
dev.off()

##################### DotPlot for cell_type
pdf("DotPlot_for_cell_type.pdf",width = 5.5,height = 5)
DotPlot(single_data,features=unique(rev(c("SLC17A7","SATB2",##Excitatory neurons
                                          "GAD1","GAD2", ##Inhibitory neurons
                                          #"RELN","GLCE", ##Purkinje neurons
                                          "MOBP","MBP","MOG","PLP1", ##Oligodendrocytes
                                          "PDGFRA","MEGF11", ##OPCs
                                          "GJA1","SLC4A4",##Astrocytes
                                          "APBB1IP","P2RY12","FYB", ##Microglia
                                          "PECAM1","FLT1","CLDN5" ##Endothelial cells 
                                                   ))),group.by="cell_type",assay="RNA")+
  scale_colour_gradient2(low="#3A71AA",mid="white",high="#B22028",midpoint=0)+
  theme_bw()+
  coord_flip()+
  theme(axis.text.x = element_text(size = 10,color = 'black',angle = 35,hjust = 1,vjust = 1), 
        axis.text.y = element_text(size = 10,color = 'black'), 
        legend.text = element_text(size = 10,color = 'black'), 
        legend.title = element_text(size = 10,color = 'black'), 
        panel.grid.major =  element_blank(),
        legend.position = 'right')
dev.off()


##################### calculate the cell numbers for each celltypes
cell_number<-as.data.frame(table(single_data$cell_type))
cell_number<-cell_number[order(cell_number$Freq,decreasing = T),]
cell_number$Var1<-factor(cell_number$Var1,levels = cell_number$Var1)

pdf("The_cell_numbers_in_Singlecell_mouse_data.pdf",height = 4.5,width = 5)
ggplot(cell_number,aes(Var1,Freq))+
  geom_col(aes(fill=Var1),width = 0.6)+
  scale_fill_manual(values = c("#4ea64a","#8e4c99","#e47faf","#d5231d","#b698c5","#3777ac","#e88f18"))+  
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


##################### calculate the cell numbers for each brain region
cell_number_brain_region<-as.data.frame(table(single_data$brain_region))
cell_number_brain_region<-cell_number_brain_region[order(cell_number_brain_region$Freq,decreasing = T),]
cell_number_brain_region$Var1<-factor(cell_number_brain_region$Var1,levels = cell_number_brain_region$Var1)

pdf("The_cell_numbers_in_brain_region_Singlecell_mouse_data.pdf",height = 4.5,width = 5)
ggplot(cell_number_brain_region,aes(Var1,Freq))+
  geom_col(aes(fill=Var1),width = 0.6)+
  scale_fill_manual(values = c("#CAA2F4","#1F78B4","#ED4437","#B49D99","#E1884A"))+  
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






##################### scPagwas for cell_type and brain region
setwd("D:/Project/SingleCell_MDD/Figure/Singlecell_Mouse/scPagwas/")
scPagwas_result<-read.csv("singlecell_scPagwas_score_pvalue_Result1.csv",row.names = 1)
single_data$scPagwas.TRS.Score<-scPagwas_result$scPagwas.TRS.Score
#### get the plot data
scPagwas_result1<-FetchData(single_data,vars = c("scPagwas.TRS.Score","cell_type","brain_region"))

##################### boxplot for celltype
scPagwas_result1_celltype<-aggregate(scPagwas_result1$scPagwas.TRS.Score,by=list(scPagwas_result1$cell_type),median)
scPagwas_result1_celltype<-scPagwas_result1_celltype[order(scPagwas_result1_celltype$x,decreasing = T),]
scPagwas_result1$cell_type<-factor(scPagwas_result1$cell_type,levels = scPagwas_result1_celltype$Group.1)

pdf("Singlecell_mouse_scPagwas_ecore_for_celltype.pdf",height = 4.5,width = 5.5)
ggplot(scPagwas_result1, aes(x=cell_type, y=scPagwas.TRS.Score)) + 
  stat_boxplot(geom = "errorbar",width=0.05, size=0.5,position=position_dodge(0.6),color= "black")+
  theme_bw()+
  geom_boxplot(position = position_dodge(0.6),
               size = 0.5,
               width = 0.7,
               fill = c("#8e4c99","#4ea64a","#b698c5","#d5231d","#e47faf","#e88f18","#3777ac"),
               color = "black",
               outlier.color = "black",
               outlier.fill = "black",
               outlier.shape = 19,
               outlier.size = 1,
               outlier.stroke = 0.5,
               outlier.alpha = 45,
               notch = F,
               notchwidth = 0.5)+
  xlab("Brain region")+
  ylab("scPagwas.TRS.Score")+
  theme(axis.text.x = element_text(size = 10,color = 'black',vjust = 1,hjust = 1,angle = 35),
        axis.text.y = element_text(size = 10,color = 'black'), 
        legend.text = element_text(size = 10,color = 'black'), 
        legend.title = element_text(size = 10,color = 'black'), 
        panel.grid.major =  element_blank(),
        panel.grid.minor = element_blank())
dev.off()



##################### boxplot for brain region
scPagwas_result1_brain_region<-aggregate(scPagwas_result1$scPagwas.TRS.Score,by=list(scPagwas_result1$brain_region),median)
scPagwas_result1_brain_region<-scPagwas_result1_brain_region[order(scPagwas_result1_brain_region$x,decreasing = T),]
scPagwas_result1$brain_region<-factor(scPagwas_result1$brain_region,levels = scPagwas_result1_brain_region$Group.1)

pdf("Singlecell_mouse_scPagwas_ecore_for_brain_region.pdf",height = 4.5,width = 5.5)
ggplot(scPagwas_result1, aes(x=brain_region, y=scPagwas.TRS.Score)) + 
  stat_boxplot(geom = "errorbar",width=0.05, size=0.5,position=position_dodge(0.6),color= "black")+
  theme_bw()+
  geom_boxplot(position = position_dodge(0.6),
               size = 0.5,
               width = 0.65,
               fill = c("#CAA2F4","#1F78B4","#ED4437","#E1884A","#B49D99"),
               color = "black",
               outlier.color = "black",
               outlier.fill = "black",
               outlier.shape = 19,
               outlier.size = 1,
               outlier.stroke = 0.5,
               outlier.alpha = 45,
               notch = F,
               notchwidth = 0.5)+
  xlab("Brain region")+
  ylab("scPagwas.TRS.Score")+
  theme(axis.text.x = element_text(size = 10,color = 'black',vjust = 1,hjust = 1,angle = 35),
        axis.text.y = element_text(size = 10,color = 'black'), 
        legend.text = element_text(size = 10,color = 'black'), 
        legend.title = element_text(size = 10,color = 'black'), 
        panel.grid.major =  element_blank(),
        panel.grid.minor = element_blank())
dev.off()

