

#----Step 1 human brain cerebral organoid merge .sh
```sh
#PBS -N Brain_organoid
#PBS -q fat
#PBS -l nodes=fat03
#PBS -l mem=450gb
#PBS -l ncpus=5
cd /share2/pub/chenchg/chenchg/SingleCell/Brain/organoid/Sample_merge/
source activate R4.2.2
Rscript Brain_organoid.R
```

#----Step 2 human brain cerebral organoid merge .r
```R
############### Brain organoid singlecell data merge ################
rm(list = ls())
library(Seurat)
library(harmony)
library(scater)
library(stringr)
library(ggplot2)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(openxlsx)
set.seed(123)
setwd("/share2/pub/chenchg/chenchg/SingleCell/Brain/organoid/Sample_merge/")

############### import the singlecell data of different samples
scrna.list<-list()
for (i in str_split_fixed(list.files("/share2/pub/chenchg/chenchg/SingleCell/Brain/organoid/Sample_rds/"),pattern = "\\.",n=2)[,1]) {
  scrna.list[[i]]<-readRDS(paste0("/share2/pub/chenchg/chenchg/SingleCell/Brain/organoid/Sample_rds/",i,".rds"))
  scrna.list[[i]]$sample<-i
}
############### merge all samples
pancreas_merged <- merge(scrna.list[[1]], y = scrna.list[2:length(scrna.list)], project = "merged", merge.data = TRUE)

############### change the gene ID from ENSG to gene symbol
count<-pancreas_merged@assays$RNA@counts
count<-as.matrix(count)
symbol<-scrna.list$`HMO-1`@assays[["RNA"]]@meta.features[["gene_symbols"]]
row.names(count)<-symbol
metadata<-pancreas_merged@meta.data

############### Create SeuratObject and run Seurat pipeline
single_data <- CreateSeuratObject(counts=count,min.cells = 0,min.features = 0,meta.data = metadata,project = "merge")
single_data <- NormalizeData(single_data, normalization.method = "LogNormalize", scale.factor = 10000)
single_data <- CellCycleScoring(single_data, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes, set.ident = TRUE)
single_data <- FindVariableFeatures(single_data, assay = single_data@active.assay, selection.method = "vst", nfeatures = 3000)
single_data <- ScaleData(single_data)
single_data <- RunPCA(object = single_data, assay = single_data@active.assay,  npcs = 50)
single_data <- RunHarmony(object = single_data,
                          assay.use = single_data@active.assay,
                          reduction.use="pca",
                          dims.use = 1:50,
                          group.by.vars = "sample",
                          plot_convergence = TRUE)
single_data <- RunUMAP(object = single_data, assay = single_data@active.assay, reduction = "harmony", dims = 1:50)
single_data <- RunTSNE(object = single_data, assay = single_data@active.assay, reduction = "harmony", dims = 1:50)

############## FindClusters
single_data <- FindNeighbors(object = single_data, assay = single_data@active.assay, reduction = "harmony", dims = 1:50)
single_data <- FindClusters(object = single_data, resolution = 0.1)
single_data <- FindClusters(object = single_data, resolution = 0.2)
single_data <- FindClusters(object = single_data, resolution = 0.3)
single_data <- FindClusters(object = single_data, resolution = 0.4)
single_data <- FindClusters(object = single_data, resolution = 0.5)

############## save the merged rds
saveRDS(single_data, file = "Brain_organoid.rds")


############## celltype annotation
setwd("D:/Project/SingleCell_MDD/organoid/Sample_merge/")
single_data<-readRDS("Brain_organoid.rds")
table(single_data$sample)
DimPlot(single_data,group.by="sample",reduction = "umap",raster=FALSE)

############## import the marker gene
marker_gene<-read.xlsx("D:/Project/SingleCell_MDD/大脑单细胞数据标准marker基因.xlsx")

##### DimPlot for origional celltypes
table(single_data$annotation2)
pdf("DimPlot_for_annotation2.pdf",height = 6,width = 6.5)
DimPlot(single_data,group.by="annotation2",reduction = "umap",raster=FALSE,label = T,repel = T)
dev.off()

##### DimPlot for seurat clusters
table(single_data$RNA_snn_res.0.3)
single_data$RNA_snn_res.0.3<-factor(single_data$RNA_snn_res.0.3,levels = rep(0:20))
pdf("DimPlot_for_seurat_clusters.pdf",height = 6,width = 6.5)
DimPlot(single_data,group.by="RNA_snn_res.0.3",reduction = "umap",raster=FALSE,label = T,repel = T)
dev.off()

##### DotPlot for seurat clusters
pdf("DotPlot_for_seurat_clusters.pdf",height = 8,width = 6.5)
DotPlot(single_data,features=unique(rev(c("SLC17A7", "SATB2","CBLN2", "LDB2", "NEUROD2", "NEUROD6","FOXP2", "SLC17A6", "GBX2", "PBX3", "RBFOX1",
                                          "GAD1", "GAD2","RELN", "GLCE", "DLX2", "PDE4DIP", "SLC32A1", "CALB2",
                                          "MOBP", "MBP", "MOG", "PLP1","SOX10", ## Oligodendrocytes
                                          "PDGFRA", "MEGF11", "CSPG4","COL20A1", "PMP2", ## OPCs
                                          "GJA1", "SLC4A4", "NTSR2", "AQP4", "GFAP", "ALDH1L1", "GJB6", "SLCO1C1", ## Astrocytes
                                          "APBB1IP", "P2RY12", "CX3CR1", "C1QA", "TREM2", "TLR4", "AIF1", "PTPRC", ## Microglia
                                          "PECAM1", "FLT1", "RAMP2","ITM2A", ## Endothelial cells
                                          "EPCAM","KRT8","KRT18", ## Epithelial cells
                                          "COL1A1", "COL1A2", "ACTA2", "DCN", ## Fibroblasts
                                          "OTX2", "TRPM3", "RSPO2", "TTR", ## IPCs
                                          "TOP2A", "MKI67", "CENPF", "UBE2C", ## NPCs
                                          "GLI3", "PAX6", "SOX9", "HOPX", "FAM107A", "NES" ## RG
                                          ))),group.by="RNA_snn_res.0.3",assay="RNA")+
  scale_colour_gradient2(low="#3A71AA",mid="white",high="#B22028",midpoint=0)+
  theme_bw()+
  coord_flip()+
  theme(axis.text.x = element_text (size = 10,color = 'black',angle = 35,vjust = 1, hjust = 1),
        axis.text.y = element_text(size = 10,color = 'black'), 
        legend.text = element_text(size = 10,color = 'black'), 
        legend.title = element_text(size = 10,color = 'black') 
  )
dev.off()



single_data$cell_type<-single_data$RNA_snn_res.0.3
single_data$cell_type<-recode(single_data$cell_type,
                              "0"="Inhibitory neurons",
                              "1"="Excitatory neurons",
                              "2"="Inhibitory neurons",
                              "3"="Astrocytes",
                              "4"="RGs",
                              "5"="NPCs",
                              "6"="Endothelial cells",
                              "7"="IPCs",
                              "8"="Excitatory neurons",
                              "9"="RGs",
                              "10"="Inhibitory neurons",
                              "11"="Excitatory neurons",
                              "12"="Endothelial cells",
                              "13"="Astrocytes",
                              "14"="Epithelial cells",
                              "15"="RGs",
                              "16"="Endothelial cells",
                              "17"="Epithelial cells",
                              "18"="Astrocytes",
                              "19"="Astrocytes",
                              "20"="Microglia")


single_data$cell_type<-factor(single_data$cell_type,levels = c("Excitatory neurons","Inhibitory neurons","Astrocytes","Microglia","Endothelial cells",
                                                               "Epithelial cells","IPCs","NPCs","RGs"))

cell_type_color<-c("#4ea64a","#8e4c99","#d5231d","#e88f18","#3777ac","#e47faf","#b698c5","#a05528","#B3446C")
pdf("DimPlot_for_cell_type.pdf",height = 5.5,width = 6.3)
DimPlot(single_data,group.by="cell_type",reduction = "umap",cols=cell_type_color,raster=FALSE,label = T,repel = T)
dev.off()

pdf("DimPlot_for_sample.pdf",height = 5.5,width = 6.7)
DimPlot(single_data,group.by="sample",reduction = "umap",raster=FALSE,label = F)
dev.off()


pdf("DotPlot_for_cell_type.pdf",height = 5.5,width = 6)
DotPlot(single_data,features=unique(rev(c("SLC17A7", "SATB2", "LDB2", "NEUROD2", "NEUROD6",
                                          "GAD1", "GAD2","RELN",
                                          "GJA1", "SLC4A4", "NTSR2", "AQP4", ## Astrocytes
                                          "C1QA", "TREM2", "TLR4", "AIF1", "PTPRC", ## Microglia
                                          "PECAM1", "FLT1", "RAMP2","ITM2A", ## Endothelial cells
                                          "EPCAM","KRT8","KRT18", ## Epithelial cells
                                          "OTX2", "TRPM3", "RSPO2", "TTR", ## IPCs
                                          "TOP2A", "MKI67", "CENPF", "UBE2C", ## NPCs
                                          "GLI3", "PAX6", "SOX9", "FAM107A", "NES" ## RG
))),group.by="cell_type",assay="RNA")+
  scale_colour_gradient2(low="#3A71AA",mid="white",high="#B22028",midpoint=0)+
  theme_bw()+
  coord_flip()+
  theme(axis.text.x = element_text (size = 10,color = 'black',angle = 35,vjust = 1, hjust = 1),
        axis.text.y = element_text(size = 10,color = 'black'), 
        legend.text = element_text(size = 10,color = 'black'), 
        legend.title = element_text(size = 10,color = 'black') 
  )
dev.off()

########### save the result
saveRDS(single_data,file = "Brain_organoid.rds")


```

#----Step 3 human brain cerebral organoid Seurat format plot .r
```R
################# Brain organoid scPagwas plot ####################
library(ggplot2)
library(Seurat)
library(ggplot2)
library(ggvenn) 
library(cowplot)
library(ggpubr)
library(Nebulosa)
library(BiocFileCache)
library(paletteer)
library(scCustomize)
setwd("D:/Project/SingleCell_MDD/organoid/scPagwas/")

################# import the singlecell data
single_data<-readRDS("D:/Project/SingleCell_MDD/organoid/Sample_merge/Brain_organoid.rds")
single_data
################# import the scPagwas result
scPagwas<-read.csv("singlecell_scPagwas_score_pvalue_Result.csv",row.names = 1)

################# import the metadata
identical(row.names(scPagwas),colnames(single_data))
single_data$scPagwas.TRS.Score<-scPagwas$scPagwas.TRS.Score
metadata<-FetchData(single_data,vars = c("cell_type","scPagwas.TRS.Score"))

################# boxplot for celltypes
metadata1<-aggregate(metadata$scPagwas.TRS.Score,by=list(metadata$cell_type),median)
metadata1$color<-c("#4ea64a","#8e4c99","#d5231d","#e88f18","#3777ac","#e47faf","#b698c5","#a05528","#B3446C")
metadata1<-metadata1[order(metadata1$x,decreasing = T),]
metadata$cell_type<-factor(metadata$cell_type,levels = metadata1$Group.1)

pdf("boxplot_scPagwas_score_for_cell_type.pdf",height = 5,width = 6)
ggplot(metadata, aes(x=cell_type, y=scPagwas.TRS.Score)) + 
  stat_boxplot(geom = "errorbar",width=0.05, size=0.5,position=position_dodge(0.6),color= "black")+
  theme_bw()+
  geom_boxplot(position = position_dodge(0.6),
               size = 0.5,
               width = 0.85,
               fill = metadata1$color,
               color = "black",
               outlier.color = "black",
               outlier.fill = "black",
               outlier.shape = 19,
               outlier.size = 1,
               outlier.stroke = 0.5,
               outlier.alpha = 45,
               notch = F,
               notchwidth = 0.5)+
  xlab("Cell type")+
  ylab("scPagwas TRS Score")+
  theme(axis.text.x = element_text(size = 10,color = 'black',vjust = 1,hjust = 1,angle = 35), 
        axis.text.y = element_text(size = 10,color = 'black'), 
        legend.text = element_text(size = 10,color = 'black'), 
        legend.title = element_text(size = 10,color = 'black'),
        panel.grid.major =  element_blank(),
        panel.grid.minor = element_blank())
dev.off()


############# density plot for scPagwas.TRS.Score
pdf("Density_plot_for_scPagwas_TRS_Score.pdf",height = 5.5,width = 5.5)
Plot_Density_Custom(seurat_object =single_data, features = "scPagwas.TRS.Score",reduction = "umap",
                    custom_palette = c("#B0CFE4","#FACABC","#E77A77","#DC0000FF"))+
  theme(plot.title = element_text(hjust = 0.5))
dev.off()


############# scPagwas celltypes result
scPagwas_celltype<-read.csv("celltype_result.csv",row.names = 1)
scPagwas_celltype<-scPagwas_celltype[order(scPagwas_celltype$pvalue),]
colnames(metadata1)[1]<-"celltype"
scPagwas_celltype<-merge(scPagwas_celltype,metadata1,by="celltype")
scPagwas_celltype$trans_P<--log10(scPagwas_celltype$pvalue)
scPagwas_celltype<-scPagwas_celltype[order(scPagwas_celltype$trans_P,decreasing = T),]
scPagwas_celltype$celltype<-factor(scPagwas_celltype$celltype,levels = rev(scPagwas_celltype$celltype))

pdf("ggdotchart_for_celltype_result.pdf",height = 4,width = 5)
ggdotchart(scPagwas_celltype, x = "celltype", y = "trans_P",
           color = "celltype",                              
           palette = rev(scPagwas_celltype$color),
           sorting = "descending",                       # Sort value in descending order
           add = "segments",                             # Add segments from y = 0 to dots
           add.params = list(color = "lightgray", size = 2), # Change segment color and size
           group = "celltype",                                # Order by groups
           dot.size = 7,                                 # Large dot size
           label = round(scPagwas_celltype$trans_P,3), # Add mpg values as dot labels
           font.label = list(color = "black", size = 8, vjust = 0.5),               # Adjust label parameters
           ggtheme = theme_pubr())+
  coord_flip()+
  geom_hline(yintercept = -log10(0.05), linetype = 2, color = "lightgray")+
  theme(axis.text.x = element_text (size = 10,color = 'black',angle = 0,hjust = 0.5),
        axis.text.y = element_text(size = 10,color = 'black'),
        legend.text = element_text(size = 10,color = 'black'),
        legend.title = element_text(size = 10,color = 'black'),
        legend.position = 'none'
  )+
  ylab("-log10(P)")+
  scale_y_continuous(limits = c (0, 10), breaks = seq (0, 10, 2))
dev.off()


```



