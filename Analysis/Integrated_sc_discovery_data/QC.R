library("Seurat")
library(scater)
s.genes<-cc.genes.updated.2019$s.genes
g2m.genes<-cc.genes.updated.2019$g2m.genes
library(sctransform)
library(stringr)
setwd("/share/pub/qiuf/brain/01-data/healthy/adult/QC")
file<-list.files()
name<-str_split_fixed(file,pattern=".rds",n=2)[,1]
for(i in 1:length(file)){
data<-readRDS(file[i])
a<-isOutlier(data@meta.data$nFeature_RNA, log=TRUE)
b<-isOutlier(data@meta.data$nCount_RNA, log=TRUE)
c<-isOutlier(data@meta.data$percent.mt,type="higher")
d<-isOutlier(data@meta.data$percent.rp,type="higher")
data<-data[,!(a|b|c|d)]
data<- NormalizeData(data, normalization.method = "LogNormalize", scale.factor = 10000)
data<- CellCycleScoring(data,  s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
data<-SCTransform(data, vars.to.regress = c('nFeature_RNA', 'nCount_RNA',"percent.mt","percent.rp","S.Score", "G2M.Score"), verbose = FALSE)
print(name[i])
print(dim(data))
saveRDS(data,paste0(name[i],".rds"))
}

