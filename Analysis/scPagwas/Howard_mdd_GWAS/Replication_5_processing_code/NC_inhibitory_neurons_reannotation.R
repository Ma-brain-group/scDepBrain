################# nature communication singlecell data reannotation ####################
library(Seurat)
library(ggplot2)
library(harmony)
library(dplyr)
setwd("D:/Project/SingleCell_MDD/SingleCell_analysis/16_validation_NC/2_INs_subset/")
single_data<-readRDS("D:/Project/SingleCell_MDD/SingleCell_analysis/16_validation_NC/1_NC_data_process/GSE213982/2023_NC.rds")

################ 160222 cells
single_data
DimPlot(single_data,group.by = "sample",reduction = "umap",label = F,raster = F)
DimPlot(single_data,group.by = "batch",reduction = "umap",label = F,raster = F)
DimPlot(single_data,group.by = "broad_cell_type",reduction = "umap",label = T,raster = F)
DimPlot(single_data,group.by = "subtype",reduction = "umap",label = T,raster = F)


################ change the raw broad_cell_type to fit for our discovery singlecell data
table(single_data$broad_cell_type)
Idents(single_data)<-single_data$broad_cell_type
single_data$broad_cell_type[WhichCells(single_data,idents = "ExN")]<-"Excitatory neurons"
single_data$broad_cell_type[WhichCells(single_data,idents = "InN")]<-"Inhibitory neurons"
single_data$broad_cell_type[WhichCells(single_data,idents = "End")]<-"Endothelial cells"
single_data$broad_cell_type[WhichCells(single_data,idents = "Ast")]<-"Astrocytes"
single_data$broad_cell_type[WhichCells(single_data,idents = "Oli")]<-"Oligodendrocytes"
single_data$broad_cell_type[WhichCells(single_data,idents = "OPC")]<-"OPCs"
single_data$broad_cell_type[WhichCells(single_data,idents = "Mic")]<-"Microglia"
DimPlot(single_data,group.by = "broad_cell_type",reduction = "umap",label = T,repel = T,raster = F)

################ remove the Mix celltypes
single_data<-subset(single_data,idents = "Mix",invert=T)
table(single_data$broad_cell_type)


############### the expression of CNNM2 in broad celltypes
pdf("the_expression_of_CNNM2_in_broad_celltypes.pdf",height = 5)
DotPlot(single_data,features=c("CNNM2"),group.by="broad_cell_type",assay="RNA")+
  scale_colour_gradient2(low="#3A71AA",mid="white",high="#B22028",midpoint=0)+
  theme_bw()+
  coord_flip()+
  theme(axis.text.x = element_text (size = 10,color = 'black',angle = 35,vjust = 1, hjust = 1),
        axis.text.y = element_text(size = 10,color = 'black'),
        legend.text = element_text(size = 10,color = 'black'),
        legend.title = element_text(size = 10,color = 'black'))
dev.off()


############### validation using statistical methods
############### method1 FindMarkers
Idents(single_data)<-single_data$broad_cell_type
DEG_findmarker<-FindMarkers(single_data,assay = "RNA",slot = "data",group.by = "broad_cell_type",ident.1 = "Inhibitory neurons",
                            only.pos = T,logfc.threshold = 0.2)
DEG_findmarker<-DEG_findmarker[DEG_findmarker$p_val_adj<0.05,]
DEG_findmarker[row.names(DEG_findmarker)=="CNNM2",]
write.csv(DEG_findmarker,file = "findmarker_INs_vs_broad_cell_type_DEG.csv",quote = F)

############### method2 COSG 
library(COSG)
library(dplyr)
DEG_COSG <- cosg(single_data,
                 groups='all',
                 assay='RNA',
                 slot='data',
                 mu=1,
                 n_genes_user=500)
DEG_COSG<-as.data.frame(DEG_COSG)
#write.csv(DEG_COSG,file = "COSG_DEG.csv",quote = F)


################ subset the Inhibitory neurons and reannotation
single_data_INs <- subset(single_data,idents = "Inhibitory neurons")
single_data_INs <- NormalizeData(single_data_INs, normalization.method = "LogNormalize",scale.factor = 10000)
single_data_INs <- FindVariableFeatures(single_data_INs, selection.method = "vst", nfeatures = 2000)
single_data_INs <- ScaleData(single_data_INs, features = VariableFeatures(single_data_INs))
single_data_INs <- RunPCA(single_data_INs, features = VariableFeatures(single_data_INs))
single_data_INs <- RunHarmony(single_data_INs, group.by.vars="sample", plot_convergence = TRUE)
ElbowPlot(single_data_INs)
single_data_INs<-FindNeighbors(single_data_INs,reduction = "harmony", dims = 1:15)
single_data_INs<-FindClusters(single_data_INs,resolution = 0.1)
single_data_INs<-RunUMAP(single_data_INs,reduction = "harmony", dims = 1:15)
DimPlot(single_data_INs,group.by = "RNA_snn_res.0.1",reduction = "umap",label = T,repel = T)
saveRDS(single_data_INs,file = "NC_inhibitory_neurons_reannotation.rds")


################ import the singlecell data of Inhibitory neurons
single_data_INs<-readRDS("NC_inhibitory_neurons_reannotation.rds")
DimPlot(single_data_INs,group.by = "RNA_snn_res.0.1",reduction = "umap",label = T,repel = T)



################ annotate the Inhibitory neurons subtypes to fit for discovery singlecell data
DotPlot(single_data_INs,features=rev(c("VIP", ## In_VIP
                                       "SST", ## In_SST
                                       "PVALB", ## In_PVALB
                                       "LAMP5", ## In_LAMP5
                                       "RELN" ## In_Reelin
                                     )),group.by="RNA_snn_res.0.1",assay="RNA")+
  scale_colour_gradient2(low="#3A71AA",mid="white",high="#B22028",midpoint=0)+
  theme_bw()+
  coord_flip()+
  theme (axis.text.x = element_text (angle = 35, hjust = 1))

################ 
marker5<-FindMarkers(single_data_INs,group.by = "RNA_snn_res.0.1",ident.1 = "5",only.pos = T)
marker5<-marker5[marker5$p_val_adj<0.05,]

Idents(single_data_INs)<-single_data_INs$seurat_clusters
single_data_INs$RNA_snn_res.0.1_cell_type<-single_data_INs$seurat_clusters
single_data_INs$RNA_snn_res.0.1_cell_type<-recode(single_data_INs$RNA_snn_res.0.1_cell_type,
                                                 "0"="In_PVALB",
                                                 "1"="In_VIP",
                                                 "2"="In_SST",
                                                 "3"="In_LAMP5",
                                                 "4"="In_Reelin",
                                                 "5"="In_mix",
                                                 "6"="In_LAMP5",
                                                 "7"="In_PVALB",
                                                 "8"="In_Reelin",
                                                 "9"="In_SST"
                                                 )

################ change the character to factor 
single_data_INs$RNA_snn_res.0.1_cell_type<-factor(single_data_INs$RNA_snn_res.0.1_cell_type,levels = c("In_VIP","In_SST","In_PVALB","In_LAMP5",
                                                                                                       "In_Reelin","In_mix"))

################ plot the result
color<-c("#396046","#A04F1D","#5D5495","#A3236D","#658E2D","#9FB8CE","#3C5A93",
         "#CC8482","#A11F1F","#525252","#9EBE70","#4E8730")

pdf("NC_inhibitory_neurons_reannotation_dimplot.pdf",height =6,width = 6.5)
DimPlot(single_data_INs,group.by = "RNA_snn_res.0.1_cell_type",reduction = "umap",cols = color)                                                                                                       
dev.off()                                                                                                       
                                                                                                       
pdf("NC_inhibitory_neurons_reannotation_marker_gene_dotplot.pdf",height = 6,width = 6)
DotPlot(single_data_INs,features=rev(c("VIP", ## In_VIP
                                       "SST", ## In_SST
                                       "PVALB", ## In_PVALB
                                       "LAMP5", ## In_LAMP5
                                       "RELN" ## In_Reelin
                                       )),group.by="RNA_snn_res.0.1_cell_type",assay="RNA")+
  scale_colour_gradient2(low="#3A71AA",mid="white",high="#B22028",midpoint=0)+
  theme_bw()+
  coord_flip()+
  theme(axis.text.x = element_text (size = 10,color = 'black',angle = 35,vjust = 1, hjust = 1),
        axis.text.y = element_text(size = 10,color = 'black'),
        legend.text = element_text(size = 10,color = 'black'),
        legend.title = element_text(size = 10,color = 'black'))
dev.off()

pdf("NC_CNNM2_in_inhibitory_neurons_dotplot.pdf",height = 4)
DotPlot(single_data_INs,features="CNNM2",group.by="RNA_snn_res.0.1_cell_type",assay="RNA")+
  scale_colour_gradient2(low="#3A71AA",mid="white",high="#B22028",midpoint=-0.5)+
  theme_bw()+
  coord_flip()+
  theme(axis.text.x = element_text (size = 10,color = 'black',angle = 35,vjust = 1, hjust = 1),
        axis.text.y = element_text(size = 10,color = 'black'),
        legend.text = element_text(size = 10,color = 'black'),
        legend.title = element_text(size = 10,color = 'black'))
dev.off()


############## 
#Idents(single_data_INs)<-single_data_INs$RNA_snn_res.0.1_cell_type
#single_data_INs<-subset(single_data_INs,idents = "In_mix",invert=T)
#marker<-FindMarkers(single_data_INs,assay = "RNA",slot = "data",group.by = "RNA_snn_res.0.1_cell_type",ident.1 = "In_PVALB",logfc.threshold = 0,min.pct = 0)
#marker<-marker[marker$p_val_adj<0.05,]
#marker[row.names(marker)=="CNNM2",]

############## save the result
saveRDS(single_data,file = "D:/Project/SingleCell_MDD/SingleCell_analysis/16_validation_NC/1_NC_data_process/GSE213982/2023_NC.rds")
saveRDS(single_data_INs,file = "NC_inhibitory_neurons_reannotation.rds")




