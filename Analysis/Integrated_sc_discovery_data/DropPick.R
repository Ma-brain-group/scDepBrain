rm(list=ls())
gc()
setwd("E:/brain/01-data/healthy/data")
setwd("E:/brain/01-data/healthy")

#ϸ�����ڻ���
s.genes<-cc.genes.updated.2019$s.genes
g2m.genes<-cc.genes.updated.2019$g2m.genes
##
file<-list.files("GSE140231")
GSE140231<-list()
for(i in 1:length(file)){
  GSE140231[[i]]<-Read10X(paste0("GSE140231","/",file[i]))
  GSE140231[[i]] = CreateSeuratObject(counts = GSE140231[[i]], min.cells=3, min.features=200, project=substr(file[i],12,30));
  GSE140231[[i]][["percent.mt"]] <- PercentageFeatureSet(GSE140231[[i]], pattern = "^MT-")
  GSE140231[[i]][["percent.rp"]]  = PercentageFeatureSet(GSE140231[[i]], pattern = "^RP[SL][[:digit:]]")
}
###��Դ����δ�ἰԤ����
GSE140231<-merge(GSE140231[[1]],GSE140231[c(2:length(GSE140231))],project="GSE140231", merge.data = TRUE)#17974cells
#��meta.data
meta.data<-read.delim("E:/brain/01-data/healthy/series_matrix/GSE140231.txt", header=FALSE, row.names=1)
GSE140231@meta.data<-cbind(GSE140231@meta.data,t(meta.data[2:nrow(meta.data),match(GSE140231@meta.data$orig.ident,meta.data[1,])]))
GSE140231@meta.data$sample<-GSE140231@meta.data$orig.ident
GSE140231@meta.data$orig.ident<-"GSE140231"
GSE140231@meta.data$batch<-"GSE140231"
saveRDS(GSE140231,"GSE140231.rds")

#GSE126836
rm(list=ls())
gc()
s.genes<-cc.genes.updated.2019$s.genes
g2m.genes<-cc.genes.updated.2019$g2m.genes
file<-list.files("GSE126836")
GSE126836<-list()
for(i in 1:length(file)){
  list<-list.files(paste0("GSE126836","/",file[i]))
  GSE126836[[i]]<-readMM(paste0("GSE126836","/",file[i],"/",list[3]))
  colnames(GSE126836[[i]])<-t(read.csv(paste0("GSE126836","/",file[i],"/",list[1]),header = T))
  row.names(GSE126836[[i]])<-t(read.csv(paste0("GSE126836","/",file[i],"/",list[2]),header = T))
  GSE126836[[i]] = CreateSeuratObject(counts = GSE126836[[i]], min.cells=3, min.features=200, project=file[i]);
  GSE126836[[i]][["percent.mt"]] <- PercentageFeatureSet(GSE126836[[i]], pattern = "^MT-")
  GSE126836[[i]][["percent.rp"]]  = PercentageFeatureSet(GSE126836[[i]], pattern = "^RP[SL][[:digit:]]")
}
GSE126836<-merge(GSE126836[[1]],GSE126836[c(2:length(GSE126836))],project="GSE126836", merge.data = TRUE)#40453cells
meta.data<-read.delim("E:/brain/01-data/healthy/series_matrix/GSE126836.txt", header=FALSE, row.names=1)
GSE126836@meta.data<-cbind(GSE126836@meta.data,t(meta.data[c(1,3),match(substr(GSE126836@meta.data$orig.ident,1,6),meta.data[2,])]),batch="GSE126836")
GSE126836@meta.data$sample<-GSE126836@meta.data$orig.ident
GSE126836@meta.data$orig.ident<-"GSE126836"
GSE126836@meta.data$age<-"adult"
saveRDS(GSE126836,"GSE126836.rds")
#GSE127774
rm(list=ls())
gc()
s.genes<-cc.genes.updated.2019$s.genes
g2m.genes<-cc.genes.updated.2019$g2m.genes
ACC<-Read10X("GSE127774/ACC")
ACC = CreateSeuratObject(counts = ACC, min.cells=3, min.features=200, project="ACC");
ACC[["percent.mt"]] <- PercentageFeatureSet(ACC, pattern = "^MT-")
ACC[["percent.rp"]]  = PercentageFeatureSet(ACC, pattern = "^RP[SL][[:digit:]]")
#ACC@meta.data$`tissue region`<-"Cingulate Anterior"
ACC@meta.data$`tissue region`<-"ACC"
ACC@meta.data$batch<-"ACC"
ACC@meta.data$orig.ident<-"GSE127774"
ACC@meta.data$sample<-"GSE127774"
ACC@meta.data$age<-"adult"
saveRDS(ACC,"GSE127774.ACC.rds")

rm(list=ls())
gc()
CER<-Read10X("GSE127774/CER")
CER = CreateSeuratObject(counts = CER, min.cells=3, min.features=200, project="CER");
CER[["percent.mt"]] <- PercentageFeatureSet(CER, pattern = "^MT-")
CER[["percent.rp"]]  = PercentageFeatureSet(CER, pattern = "^RP[SL][[:digit:]]")
#CER@meta.data$`tissue region`<-"Cerebellar Grey Matter"
CER@meta.data$`tissue region`<-"CER"
CER@meta.data$batch<-"CER"
CER@meta.data$orig.ident<-"GSE127774"
CER@meta.data$sample<-"GSE127774"
CER@meta.data$age<-"adult"
saveRDS(CER,"GSE127774.CER.rds")

rm(list=ls())
gc()
CN<-Read10X("GSE127774/CN")
CN = CreateSeuratObject(counts = CN, min.cells=3, min.features=200, project="CN")
CN[["percent.mt"]] <- PercentageFeatureSet(CN, pattern = "^MT-")
CN[["percent.rp"]] <- PercentageFeatureSet(CN, pattern = "^RP[SL][[:digit:]]")
#CN@meta.data$`tissue region`<-"Caudate"
CN@meta.data$`tissue region`<-"CN"
CN@meta.data$batch<-"CN"
CN@meta.data$orig.ident<-"GSE127774"
CN@meta.data$sample<-"GSE127774"
CN@meta.data$age<-"adult"
saveRDS(CN,"GSE127774.CN.rds")


#Allen brain
#R.R��Ԥ�����ļ�
#��ȡԤ�����ļ����.rds�ļ�
rm(list=ls())
gc()
library("Seurat")
library("Matrix")
s.genes<-cc.genes.updated.2019$s.genes
g2m.genes<-cc.genes.updated.2019$g2m.genes
M1<-readRDS("ALLEN/M1/M1.rds")
meta.data<-read.csv("ALLEN/M1/metadata.csv")
M1 = CreateSeuratObject(counts = M1, min.cells=3, min.features=200, project="M1",names.delim = "/")
M1@meta.data$`tissue region`<-meta.data[match(colnames(M1),meta.data$sample_name),]$region_label
M1@meta.data<-cbind(M1@meta.data,sex=meta.data[match(colnames(M1),meta.data$sample_name),12])
M1@meta.data$batch<-"Allen-M1"
M1@meta.data$sample<-meta.data[match(colnames(M1),meta.data$sample_name),]$external_donor_name_label
M1@meta.data$sex[M1@meta.data$sex=="M"]<-"Male"
M1@meta.data$sex[M1@meta.data$sex=="F"]<-"Female"
M1@meta.data$age<-60
M1@meta.data$age[M1@meta.data$sex=="Male"]<-50
M1[["percent.mt"]] <- PercentageFeatureSet(M1, pattern = "^MT-")
M1[["percent.rp"]] <- PercentageFeatureSet(M1, pattern = "^RP[SL][[:digit:]]")
saveRDS(M1,"M1.rds")

rm(list=ls())
gc()
Multi<-readRDS("ALLEN/Multiple/Multi.rds")
meta.data<-read.csv("ALLEN/Multiple/metadata.csv")
Multi = CreateSeuratObject(counts = Multi, min.cells=3, min.features=200, project="Multi",names.delim = "/")
Multi@meta.data$`tissue region`<-meta.data[match(colnames(Multi),meta.data$sample_name),]$region_label
Multi@meta.data<-cbind(Multi@meta.data,sex=meta.data[match(colnames(Multi),meta.data$sample_name),18])
Multi@meta.data$batch<-"Allen-Multi"
Multi@meta.data$sample<-meta.data[match(colnames(Multi),meta.data$sample_name),]$external_donor_name_label
Multi@meta.data$sex[Multi@meta.data$sex=="M"]<-"Male"
Multi@meta.data$sex[Multi@meta.data$sex=="F"]<-"Female"
Multi@meta.data$age<-50
Multi@meta.data$age[Multi@meta.data$sample=="H200.1023"]<-43
Multi@meta.data$age[Multi@meta.data$sample=="H200.1030"]<-54

Multi[["percent.mt"]] <- PercentageFeatureSet(Multi, pattern = "^MT")
Multi[["percent.rp"]] <- PercentageFeatureSet(Multi, pattern = "^RP[SL][[:digit:]]")
saveRDS(Multi,"Multi.rds")

#PRJNA434002
rm(list=ls())
gc()
library("Seurat")
library("Matrix")
s.genes<-cc.genes.updated.2019$s.genes
g2m.genes<-cc.genes.updated.2019$g2m.genes
PRJNA434002<-Read10X("PRJNA434002")
PRJNA434002 = CreateSeuratObject(counts = PRJNA434002, min.cells=3, min.features=200, project="PRJNA434002");
PRJNA434002[["percent.mt"]] <- PercentageFeatureSet(PRJNA434002, pattern = "^MT-")
PRJNA434002[["percent.rp"]]  = PercentageFeatureSet(PRJNA434002, pattern = "^RP[SL][[:digit:]]")
meta.data<-read.delim("E:/brain/01-data/healthy/data/series_matrix/PRJNA434002.txt", header=T,check.names = F)
meta.data<-meta.data[,c(1,2,3,5,6,7,8,11)]
a<-meta.data[match(colnames(PRJNA434002),meta.data[,1]),2:ncol(meta.data)]
colnames(a)<-c("anno","sample","tissue region","age","sex","disease","post-mortem interval (hours)")
PRJNA434002@meta.data<-cbind(PRJNA434002@meta.data,a)
PRJNA434002<-PRJNA434002[,PRJNA434002@meta.data$disease=="Control"]#52556cells
PRJNA434002@meta.data$batch<-"PRJNA434002"
saveRDS(PRJNA434002,"PRJNA434002.rds")


#GSE97942
rm(list=ls())
gc()
library("Seurat")
library("Matrix")
library(stringr)
s.genes<-cc.genes.updated.2019$s.genes
g2m.genes<-cc.genes.updated.2019$g2m.genes
file<-list.files("GSE97942")
GSE97942<-list()
name<-str_split_fixed(file,pattern="_",n=3)[,2]
for(i in 1:length(file)){
  GSE97942[[i]]<-read.table(paste0("GSE97942/",file[i]),header = T,row.names = 1)
  GSE97942[[i]]<-as.matrix(GSE97942[[i]])
  GSE97942[[i]]<-Matrix(GSE97942[[i]],sparse = T)
  GSE97942[[i]] = CreateSeuratObject(counts = GSE97942[[i]], min.cells=3, min.features=200, project=name[i],names.delim = "/");
}
GSE97942<-merge(GSE97942[[1]],y=GSE97942[2:length(GSE97942)],project="merge",merge.data = T)
GSE97942[["percent.mt"]] <- PercentageFeatureSet(GSE97942, pattern = "^MT")
GSE97942[["percent.rp"]]  = PercentageFeatureSet(GSE97942, pattern = "^RP[SL][[:digit:]]")
GSE97942@meta.data$batch<-"GSE97942"
meta.data<-read.delim("E:/brain/01-data/healthy/data/series_matrix/GSE97942.txt", check.names = F,header=F,row.names = 1)
a<-str_split_fixed(colnames(GSE97942),pattern="_",n=4)
meta.data[3,]<-as.numeric(meta.data[3,])
a[a[,3]!="",2]<-a[a[,3]!="",3]
GSE97942@meta.data$anno<-a[,1]
GSE97942@meta.data$sample<-a[,2]
GSE97942@meta.data<-cbind(GSE97942@meta.data,t(meta.data[3:4,match(GSE97942@meta.data$sample,meta.data[1,])]))
GSE97942@meta.data$age<-t(meta.data[3,match(GSE97942@meta.data$sample,meta.data[1,])])
saveRDS(GSE97942,"GSE97942.rds")


#PsychENCODE adult
library(Matrix)
library(Seurat)
s.genes<-cc.genes.updated.2019$s.genes
g2m.genes<-cc.genes.updated.2019$g2m.genes
load("E:/brain/01-data/healthy/data/PsychENCODE/adult/Sestan.adultHumanNuclei.Psychencode.Rdata")
row.names(umi.raw)<-unlist(lapply(strsplit(row.names(umi.raw),"\\|"),function(x){unlist(x)[[2]]}))
umi.raw<-Matrix(umi.raw,sparse = T)
umi.raw = CreateSeuratObject(counts = umi.raw, min.cells=3, min.features=200, project="PsychENCODE",names.delim = "/")
umi.raw@meta.data<-cbind(umi.raw@meta.data,anno=meta2[match(colnames(umi.raw),row.names(meta2)),5])
umi.raw@meta.data$orig.ident<-meta2[match(colnames(umi.raw),row.names(meta2)),3]
umi.raw@meta.data$`tissue region`<-"DFC"
umi.raw@meta.data$sample<-str_split_fixed(umi.raw@meta.data$orig.ident,pattern = "DFC",n=2)[,1]
umi.raw@meta.data$batch<-"Psychencode"
umi.raw@meta.data$age[umi.raw@meta.data$sample=="HSB340"]<-"19"
umi.raw@meta.data$age[umi.raw@meta.data$sample=="HSB189"]<-"36"
umi.raw@meta.data$age[umi.raw@meta.data$sample=="HSB106"]<-"64"
umi.raw@meta.data$sex<-"Male"
umi.raw<-umi.raw[,!is.na(umi.raw@meta.data$orig.ident)]
umi.raw[["percent.mt"]] <- PercentageFeatureSet(umi.raw, pattern = "^MT")
umi.raw[["percent.rp"]]  = PercentageFeatureSet(umi.raw, pattern = "^RP[SL][[:digit:]]")
saveRDS(umi.raw,"PsychENCODE.adult.rds")


#phs001836
rm(list=ls())
gc()
library(Matrix)
library(Seurat)
s.genes<-cc.genes.updated.2019$s.genes
g2m.genes<-cc.genes.updated.2019$g2m.genes
load("phs001836/raw_counts_mat.rdata")
meta.data<-read.csv("phs001836/cell_metadata.csv")
raw_counts_mat = CreateSeuratObject(counts = raw_counts_mat, min.cells=3, min.features=200, project="phs001836",names.delim = "/")
raw_counts_mat@meta.data<-cbind(raw_counts_mat@meta.data,donor=meta.data[match(colnames(raw_counts_mat),meta.data$Cell),4])
raw_counts_mat@meta.data$anno<-meta.data[match(colnames(raw_counts_mat),meta.data$Cell),2]
raw_counts_mat@meta.data$`tissue region`<-"neocortex"
raw_counts_mat@meta.data$sample<-raw_counts_mat@meta.data$donor
raw_counts_mat@meta.data$batch<-"phs001836"
raw_counts_mat@meta.data$type<-"fetal"
raw_counts_mat[["percent.mt"]] <- PercentageFeatureSet(raw_counts_mat, pattern = "^MT")
raw_counts_mat[["percent.rp"]]  = PercentageFeatureSet(raw_counts_mat, pattern = "^RP[SL][[:digit:]]")#33986
raw_counts_mat<-raw_counts_mat[,!is.na(raw_counts_mat@meta.data$anno)]
VlnPlot(raw_counts_mat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rp"), ncol = 4)#
raw_counts_mat<-subset(raw_counts_mat,subset=nFeature_RNA>200&nCount_RNA<7500&percent.mt<5&percent.rp<20)#33342
raw_counts_mat <- NormalizeData(raw_counts_mat, normalization.method = "LogNormalize", scale.factor = 10000)
raw_counts_mat <- CellCycleScoring(raw_counts_mat,  s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
raw_counts_mat<-SCTransform(raw_counts_mat, vars.to.regress = c('nFeature_RNA', 'nCount_RNA',"percent.mt","percent.rp","S.Score", "G2M.Score"), verbose = FALSE)
saveRDS(raw_counts_mat,"phs001836.fetal.rds")


#phs000424.v8.p1
##hip pfc
rm(list=ls())
gc()
library(Matrix)
library(Seurat)
s.genes<-cc.genes.updated.2019$s.genes
g2m.genes<-cc.genes.updated.2019$g2m.genes
phs000424.v8.p1<-read.delim("phs000424.v8.p1/GTEx_droncseq_hip_pcf.umi_counts.txt",row.names = 1)
phs000424.v8.p1<-as.matrix(phs000424.v8.p1)
phs000424.v8.p1<-Matrix(phs000424.v8.p1,sparse = T)
phs000424.v8.p1<-CreateSeuratObject(counts = phs000424.v8.p1, min.cells=3, min.features=200, project="phs000424.v8.p1",names.delim = "/")
#a<-unlist(lapply(strsplit(colnames(phs000424.v8.p1),"_"),function(x){unlist(x)[[1]]}))
meta.data<-read.delim("phs000424.v8.p1/meta.data.txt")
phs000424.v8.p1@meta.data$anno<-meta.data[match(colnames(phs000424.v8.p1),meta.data$Cell.ID),3]
phs000424.v8.p1@meta.data$cluster<-meta.data[match(colnames(phs000424.v8.p1),meta.data$Cell.ID),2]
phs000424.v8.p1@meta.data$batch<-"phs000424.v8.p1"
phs000424.v8.p1@meta.data$type<-"adult"
phs000424.v8.p1[["percent.mt"]] <- PercentageFeatureSet(phs000424.v8.p1, pattern = "^MT-")
phs000424.v8.p1[["percent.rp"]]  = PercentageFeatureSet(phs000424.v8.p1, pattern = "^RP[SL][[:digit:]]")#
VlnPlot(phs000424.v8.p1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rp"), ncol = 4)#
phs000424.v8.p1<-subset(phs000424.v8.p1,subset=nFeature_RNA<3000&nCount_RNA<5000&percent.mt<5&percent.rp<10)#14652
phs000424.v8.p1 <- NormalizeData(phs000424.v8.p1, normalization.method = "LogNormalize", scale.factor = 10000)
phs000424.v8.p1 <- CellCycleScoring(phs000424.v8.p1,  s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
phs000424.v8.p1<-SCTransform(phs000424.v8.p1, vars.to.regress = c('nFeature_RNA', 'nCount_RNA',"percent.mt","percent.rp","S.Score", "G2M.Score"), verbose = FALSE)
saveRDS(phs000424.v8.p1,"phs000424.v8.p1.rds")


#darmanis
rm(list=ls())
gc()
library(Matrix)
library(Seurat)
s.genes<-cc.genes.updated.2019$s.genes
g2m.genes<-cc.genes.updated.2019$g2m.genes
darmanis<-readRDS("brainscRNAseqdatabase/darmanis.rds")
darmanis <- as.Seurat(darmanis, counts = "counts", data = "logcounts")
darmanis@meta.data$`tissue region`<-darmanis@meta.data$Tissue
darmanis@meta.data$batch<-darmanis@meta.data$plate
darmanis@meta.data$orig.ident<-"darmanis"
darmanis@meta.data$sample<-darmanis@meta.data$individual
darmanis[["percent.mt"]] <- PercentageFeatureSet(darmanis, pattern = "^MT")
darmanis[["percent.rp"]]  = PercentageFeatureSet(darmanis, pattern = "^RP[SL][[:digit:]]")#
darmanis@assays$RNA<-darmanis@assays$originalexp
DefaultAssay(darmanis)<-"RNA"
darmanis@meta.data$nCount_RNA<-darmanis@meta.data$nCount_originalexp
darmanis@meta.data$nFeature_RNA<-darmanis@meta.data$nFeature_originalexp
VlnPlot(darmanis, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rp"), ncol = 4)#
darmanis <- NormalizeData(darmanis, normalization.method = "LogNormalize", scale.factor = 10000)
darmanis <- CellCycleScoring(darmanis,  s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
darmanis<-SCTransform(darmanis, vars.to.regress = c('nFeature_RNA', 'nCount_RNA',"percent.mt","percent.rp","S.Score", "G2M.Score"), verbose = FALSE)
saveRDS(darmanis,"darmanis.rds")#466

#lake
rm(list=ls())
gc()
library(Matrix)
library(Seurat)
s.genes<-cc.genes.updated.2019$s.genes
g2m.genes<-cc.genes.updated.2019$g2m.genes
lake<-readRDS("brainscRNAseqdatabase/lake.rds")
lake <- as.Seurat(lake, counts = "normcounts", data = "logcounts")
lake@meta.data$`tissue region`<-lake@meta.data$Source
lake@meta.data$orig.ident<-"lake"
lake@meta.data$sample<-"unknow"
lake[["percent.mt"]] <- PercentageFeatureSet(lake, pattern = "^MT")
lake[["percent.rp"]]  = PercentageFeatureSet(lake, pattern = "^RP[SL][[:digit:]]")#
lake@assays$RNA<-lake@assays$originalexp
DefaultAssay(lake)<-"RNA"
lake@meta.data$nCount_RNA<-lake@meta.data$nCount_originalexp
lake@meta.data$nFeature_RNA<-lake@meta.data$nFeature_originalexp
VlnPlot(lake, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rp"), ncol = 4)#
lake<-subset(lake,subset=percent.mt<20)#2991
lake <- NormalizeData(lake, normalization.method = "LogNormalize", scale.factor = 10000)
lake <- CellCycleScoring(lake,  s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
lake<-SCTransform(lake, vars.to.regress = c('nFeature_RNA', 'nCount_RNA',"percent.mt","percent.rp","S.Score", "G2M.Score"), verbose = FALSE)
saveRDS(lake,"lake.rds")

#GSE119212

GSE119212<-Read10X("GSE119212")
GSE119212<-CreateSeuratObject(counts = GSE119212, min.cells=3, min.features=200, project="GSE119212",names.delim = "/")






#33w
setwd("E:/brain/01-data/healthy/SCTdata")
file<-list.files()
data<-readRDS(file[4])
data@meta.data$batch<-"PsychENCODE"
saveRDS(data,"PsychENCODE.adult.rds")
data<-FindVariableFeatures(data,selection.method = "vst", nfeatures = 2000)
data <- RunPCA(object = data, assay = "SCT", features = VariableFeatures(data), npcs = 15)
data <- FindNeighbors(object = data, assay = "SCT", reduction = "pca", dims = 1:15)
data <- FindClusters(data,resolution =0.5, verbose = FALSE)
data <- RunUMAP(object = data, assay = "SCT", reduction = "pca", dims = 1:15)
data@meta.data$tissue.region<-data@meta.data$`tissue region`
DimPlot(data, reduction = "umap",group.by="tissue.region", pt.size = .1,label = T,raster=FALSE)
CER<-readRDS("GSE127774.CER.rds")
CER<-SCTransform(CER, vars.to.regress = c('nFeature_RNA', 'nCount_RNA',"percent.mt","percent.rp","S.Score", "G2M.Score"), verbose = FALSE)
saveRDS(CER,"GSE127774.CER.rds")
file<-list.files("E:/brain/01-data/healthy/SCTdata")

