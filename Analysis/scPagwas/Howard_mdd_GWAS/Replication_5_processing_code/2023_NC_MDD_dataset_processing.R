##2023NC 包含男性和女性的单细胞数据 GSE213982
library(Seurat)
library(dplyr)
library(patchwork)
library(harmony)
library(stringr)

##读入数据创建seurat对象
##一般10 X的单细胞数据都用Read10X读入
##必须要包含barcodes.tsv.gz，features.tsv.gz，matrix.mtx.gz这三个文件，不能改名字
setwd("/share2/pub/chenchg/chenchg/SingleCell/Brain/SingleCell_data/2023-NC/GSE213982/")
list.files()
pbmc.data <- Read10X(data.dir = getwd())
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "2023NC", min.cells = 3, min.features = 200)

##进行基本的质控
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,raster=FALSE)

##依据文章中的标准进行质控
pbmc <- subset(pbmc, subset = nFeature_RNA > 350  & nCount_RNA<120000 & percent.mt < 10)
##160222个细胞被保留
pbmc

##这里读入下载的metadata信息
library(openxlsx)
metadata<-read.xlsx("41467_2023_38530_MOESM3_ESM .xlsx",sheet = 1)
metadata<-na.omit(metadata)
colnames(metadata)<-metadata[1,]
metadata<-metadata[-1,]
metadata<-metadata[,c(1:4,6)]

##通过拆分细胞名，合并样本信息
library(stringr)
temp<-str_split_fixed(colnames(pbmc),pattern = "\\.",n=4)
colnames(temp)<-c("Sample","barcodes","broad_cell_type","subtype")
temp<-merge(temp,metadata,by="Sample")

##将metadata加入seurat对象中
##批次信息
pbmc$batch<-temp$Batch
##样本信息
pbmc$sample<-temp$Sample
##性别信息
pbmc$sex<-temp$Sex
##MDD or control
pbmc$phenotype<-temp$Condition
##Chemistry
pbmc$Chemistry<-temp$Chemistry
##broad cell type
pbmc$broad_cell_type<-temp$broad_cell_type
##subtype
pbmc$subtype<-temp$subtype


##标准化
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

##高变基因的选择,即2000 variable features
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

##scaledata
pbmc <- ScaleData(pbmc, features = VariableFeatures(object = pbmc))
gc()

##线性降维，依据之前找出的2000个高变基因进行降维
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc),npcs=100)

##RunHarmony，去批次效应
pbmc <- RunHarmony(object=pbmc,group.by.vars="sample", plot_convergence = TRUE)

##细胞聚类,文中推荐用70个PC，resolution = 0.7，k.param为30
pbmc<-FindNeighbors(pbmc,dims=1:70,k.param=30,reduction = "harmony")
pbmc <- FindClusters(pbmc, resolution = 0.7)

##非线性降维（umap,tsne）
pbmc <- RunUMAP(pbmc,reduction = "harmony",dims = 1:70)

saveRDS(pbmc,file="2023_NC.rds")



























