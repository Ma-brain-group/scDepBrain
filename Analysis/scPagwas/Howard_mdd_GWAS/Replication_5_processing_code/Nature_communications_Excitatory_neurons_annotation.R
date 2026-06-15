#################### Nature communications Excitatory neurons annotation ###################
library(Seurat)
library(harmony)
library(ggplot2)
library(dplyr)
setwd("D:/Project/SingleCell_MDD/SingleCell_analysis/Exitatory_neurons_subtype_analysis/NC_Excitatory_neurons_annotation/")

#################### import the whole singlecell data
single_data<-readRDS("D:/Project/SingleCell_MDD/SingleCell_analysis/16_validation_NC/1_NC_data_process/GSE213982/2023_NC.rds")
### subset the excitatory neurons
Idents(single_data)<-single_data$broad_cell_type
single_data_excitatory_neurons<-subset(single_data,idents = "Excitatory neurons")
DimPlot(single_data_excitatory_neurons,reduction = "umap",group.by = "subtype")

### recluster the Excitatory neurons
single_data_excitatory_neurons <- NormalizeData(single_data_excitatory_neurons, normalization.method = "LogNormalize", scale.factor = 10000)
single_data_excitatory_neurons <- FindVariableFeatures(single_data_excitatory_neurons, selection.method = "vst", nfeatures = 3000)
single_data_excitatory_neurons <- ScaleData(single_data_excitatory_neurons, verbose = FALSE)
single_data_excitatory_neurons <- RunPCA(single_data_excitatory_neurons, npcs = 20, verbose = FALSE)
single_data_excitatory_neurons <- RunHarmony(single_data_excitatory_neurons, group.by.vars="sample", plot_convergence = TRUE) 
single_data_excitatory_neurons <- RunUMAP(single_data_excitatory_neurons, reduction = "harmony", dims = 1:20)
single_data_excitatory_neurons <- RunTSNE(single_data_excitatory_neurons, reduction = "harmony", dims = 1:20)
single_data_excitatory_neurons <- FindNeighbors(single_data_excitatory_neurons, reduction = "harmony", dims = 1:20)
single_data_excitatory_neurons <- FindClusters(single_data_excitatory_neurons, resolution = 0.2)

pdf("DimPlot_RNA_snn_res.0.2.pdf")
DimPlot(single_data_excitatory_neurons,reduction = "umap",group.by = "seurat_clusters",label = T)
dev.off()

### celltypes annotation
genes <- c("NRGN", ## Ex_NRGN
           "CUX2", ## Ex-L2/3
           "RASGRF2", ## Ex-L2/3
           "RORB", ## Ex-L4
           "HTR2C",
           "FEZF2", ## Ex-5
           "ETV1", ## Ex-5/6
           "RXFP1", ## Ex-5/6
           "FOXP2", ## Ex-6
           "NR4A2", ## Ex-L6/6b
           #"SYNPR", ## Ex-L6
           "TLE4", ## Ex-L5/6
           "NTNG2" ## Ex-L6
           #"ADRA2A" ## Ex-L6b
)

pdf("DotPlot_RNA_snn_res.0.2.pdf")
DotPlot(single_data_excitatory_neurons,features=unique(rev(genes)),group.by="seurat_clusters",assay="RNA")+
  scale_colour_gradient2(low="#3A71AA",mid="white",high="#B22028",midpoint=0)+
  theme_bw()+
  coord_flip()+
  theme(axis.text.x = element_text (size = 10,color = 'black',angle = 35,vjust = 1, hjust = 1),
        axis.text.y = element_text(size = 10,color = 'black'),
        legend.text = element_text(size = 10,color = 'black'),
        legend.title = element_text(size = 10,color = 'black')
)
dev.off()

single_data_excitatory_neurons$Excitatory_neuron_subtypes<-single_data_excitatory_neurons$seurat_clusters
single_data_excitatory_neurons$Excitatory_neuron_subtypes<-recode(single_data_excitatory_neurons$Excitatory_neuron_subtypes,
                                                                  "0"="Ex-L2/4",
                                                                  "1"="Ex-NRGN",
                                                                  "2"="Ex-L2/3",
                                                                  "3"="Ex-L4/6",
                                                                  "4"="Ex-L5/6",
                                                                  "5"="Ex-L4/6",
                                                                  "6"="Ex-L2/3",
                                                                  "7"="Ex-L5/6",
                                                                  "8"="Ex-L5/6",
                                                                  "9"="Ex-L2/3",
                                                                  "10"="Ex-L2/3",
                                                                  "11"="Ex-L4/6",
                                                                  "12"="Ex-L5/6",
                                                                  "13"="Ex-L4/6",
                                                                  "14"="Ex-L6",
                                                                  "15"="Ex-L5/6",
                                                                  "16"="Ex-NRGN")

single_data_excitatory_neurons$Excitatory_neuron_subtypes<-factor(single_data_excitatory_neurons$Excitatory_neuron_subtypes,levels = c("Ex-NRGN",
                                                                  "Ex-L2/3","Ex-L2/4","Ex-L4/6","Ex-L5/6","Ex-L6"))
DimPlot(single_data_excitatory_neurons,reduction = "umap",group.by = "Excitatory_neuron_subtypes",label = T,repel = T)

pdf("DotPlot_Excitatory_neuron_subtypes.pdf",height = 4.5,width = 5)
DotPlot(single_data_excitatory_neurons,features=unique(rev(genes)),group.by="Excitatory_neuron_subtypes",assay="RNA")+
  scale_colour_gradient2(low="#3A71AA",mid="white",high="#B22028",midpoint=0)+
  theme_bw()+
  coord_flip()+
  theme(axis.text.x = element_text (size = 10,color = 'black',angle = 35,vjust = 1, hjust = 1),
        axis.text.y = element_text(size = 10,color = 'black'),
        legend.text = element_text(size = 10,color = 'black'),
        legend.title = element_text(size = 10,color = 'black')
  )
dev.off()


colorlist<-c("#f9766e","#fbbab6","#e1c548","#5fa664","#ca6a6b","#e5b5b5")
pdf("DimPlot_Excitatory_neuron_subtypes.pdf",height = 4.5,width = 5)
DimPlot(single_data_excitatory_neurons,reduction = "umap",group.by = "Excitatory_neuron_subtypes",cols = colorlist,label = T,repel = T)
dev.off()

########### save the result
saveRDS(single_data_excitatory_neurons,file = "NC_Excitatory_neuron_subtypes.rds")


########### subset the Ex-L2/4 and save the rds
Idents(single_data_excitatory_neurons)<-single_data_excitatory_neurons$Excitatory_neuron_subtypes
single_data_Ex_L2_4<-subset(single_data_excitatory_neurons,idents="Ex-L2/4")

########### define NEGR1 type
expression<-single_data_Ex_L2_4@assays$RNA@data
expression<-expression[row.names(expression)%in%c("NEGR1"),]
expression<-as.data.frame(as.matrix(expression))

colnames(expression)<-"NEGR1"
quantile(expression$NEGR1,0.3)
quantile(expression$NEGR1,0.7)
expression[expression$NEGR1<2.721594,2]<-"NEGR1-"
expression[expression$NEGR1>3.282143,2]<-"NEGR1+"

single_data_Ex_L2_4$NEGR1_type<-expression[,2]
table(single_data_Ex_L2_4$NEGR1_type)

########### save the result
saveRDS(single_data_Ex_L2_4,file = "NC_NEGR1_EX_L2_4.rds")











