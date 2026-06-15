library(Seurat)
library(dplyr)
library(plyr)
library(sctransform)
library(harmony)
setwd("/share/pub/qiuf/brain/01-data/healthy/adult/QC")
options(future.globals.maxSize= 8912896000)
#file<-list.files()
#sclist<-list()
#for(i in 1:length(file)){
  
#  sclist[[i]]<-readRDS(file[i])
  
#}
#merge<-merge(sclist[[1]],y=sclist[2:length(sclist)],project="merge",merge.data=T)
#pancreas.features <- SelectIntegrationFeatures(object.list = sclist, nfeatures = 3000)
#VariableFeatures(merge) <- pancreas.features
#merge <- RunPCA(object = merge, assay = "SCT", features = pancreas.features, npcs = 50)
#merge@meta.data$tissue.region<-merge@meta.data$`tissue region`
#saveRDS(merge,"pca.rds")
merge<-readRDS("pca.rds")
#### step3: Harmony
merge <- RunHarmony(object = merge,
                    assay.use = "SCT",
                    reduction = "pca",
                    dims.use = 1:50,
                    group.by.vars = "batch",
                    plot_convergence = TRUE)

#### step4: cell cluster
merge <- FindNeighbors(object = merge, assay = "SCT", reduction = "harmony", dims = 1:50)
merge <- FindClusters(merge,resolution =0.2, verbose = FALSE)
merge <- RunUMAP(object = merge, assay = "SCT", reduction = "harmony", dims = 1:50)
saveRDS(merge,"resolution0.2.rds")
###改变注释
library(Seurat)
library(dplyr)
library(plyr)
setwd("/share/pub/qiuf/brain/01-data/healthy/adult")
data<-readRDS("/share/pub/qiuf/brain/01-data/healthy/adult/resolution0.2_final.rds")
data$tissue.region[data$tissue.region=="CgG"]<-"ACC"
data$tissue.region[data$tissue.region=="occipital Lobe, Visual Cortex"]<-"V1C"
data$tissue.region[data$tissue.region%in%c("M1","M1lm","M1ul")]<-"M1C"
data$tissue.region[data$tissue.region%in%c("S1lm","S1ul")]<-"S1C"
#data$tissue.region[data$orig.ident=="FrontalCortex"]<-"Frontal Cortex"
#data$tissue.region[data$orig.ident=="VisualCortex"]<-"occipital Lobe, Visual Cortex"
saveRDS(data,"/share/pub/qiuf/brain/01-data/healthy/adult/resolution0.2_final1.rds")
