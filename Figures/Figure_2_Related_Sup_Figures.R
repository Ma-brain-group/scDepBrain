################# FigureS1 plot ###################
library(Seurat)
library(ggplot2)
library(corrplot)
library(cowplot)
library(dplyr)
library(plyr)
library(Nebulosa)
library(BiocFileCache)
library(paletteer)
library(scCustomize)
library(pheatmap)
library(reshape2)
setwd("D:/Project/SingleCell_MDD/Figure/FigureS1/")

################# import the singlecell data
single_data<-readRDS("D:/Project/SingleCell_MDD/SingleCell_analysis/annotation/MDD_singlecell_data_reannotation_simple.rds")

################# dimplot for brain region
table(single_data$tissue.region)
brain_region_colors<-c("#ED4437","#E1884A","#8ACC72","#1F78B4","#89C8E8","#B3446C","#EBD57C","#E68FAC","#CAA2F4",
                       "#96873B","#B49D99","#B37557","#FC9A9A","#6A3D9A")
pdf("dimplot_for_brain_region.pdf",height = 5,width = 5.3)
DimPlot(single_data,group.by = "tissue.region",reduction = "umap",raster = T,cols = brain_region_colors)+
  ggtitle("Brain region")
dev.off() 

################# dimplot for age
table(single_data$age)
my36colors <-c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
               '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
               '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
               '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
               '#712820', '#DCC1DD', '#CCE0F5',  '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
               '#968175')
single_data$age<-factor(single_data$age,levels = c("4", "6", "8", "12", "13", "14", "15", "19", "20", "21", "22", "23", "35", 
                                                   "36", "43", "48", "49", "50", "54", "55", "56", "59", "60", "64", "70", "adult"))
pdf("dimplot_for_age.pdf",height = 5,width = 6)
DimPlot(single_data,group.by = "age",reduction = "umap",raster = T,cols = my36colors)+
  ggtitle("Age")
dev.off() 


################ dimplot for cluster
pdf("dimplot_for_Cell_clusters.pdf",height = 5,width = 5.7)
DimPlot(single_data,group.by = "SCT_snn_res.0.2",reduction = "umap",raster = T,cols = my36colors)+
  ggtitle("Cell clusters")
dev.off()


################ dimplot for sample
table(single_data$sample)
pdf("dimplot_for_sample.pdf",height = 5,width = 12)
DimPlot(single_data,group.by = "sample",reduction = "umap",raster = T)+
  ggtitle("Sample")
dev.off() 


############### dimplot for celltypes in different datasets
table(single_data$orig.ident)
Idents(single_data)<-single_data$orig.ident
single_data$orig.ident[WhichCells(single_data,idents = c("HSB106DFC","HSB189DFC","HSB340DFC"))]<-"PsychENCODE"
single_data$orig.ident[WhichCells(single_data,idents = c("CerebellarHem","FrontalCortex","VisualCortex"))]<-"GSE97942"

table(single_data$final_anno)
cell_type_color<-c("#d5231d","#3777ac","#fbbab6","#e1c548","#5fa664","#abd0a7","#ca6a6b","#e5b5b5","#f9766e","#bac4d0","#8fc0dc","#e98741","#fab37f","#967568","#64abc0","#e0bc58","#e88f18","#e47faf","#b698c5","#a05528")
pdf("dimplot_for_celltypes_split_by_datasets.pdf",height = 7.5,width = 14)
DimPlot(single_data,group.by = "final_anno",reduction = "umap",split.by="orig.ident",pt.size = 2,cols = cell_type_color,raster = T,ncol = 4)+
  ggtitle("Datasets")
dev.off()


################ Cells per dataset and brain region
cell_count<-as.data.frame(table(single_data$tissue.region,single_data$orig.ident))
cell_count$Freq1<-cell_count$Freq/1000
cell_count$Freq1<-round(cell_count$Freq1,1)
cell_count$Freq1<-factor(cell_count$Freq1)
cell_count<-cell_count[cell_count$Freq1!=0,]

pdf("Cells_per_dataset_and_brain_region.pdf",height = 4,width = 5)
ggplot(data = cell_count, aes(x = Var1, y = Var2, size = Freq1, label = Freq1)) +
  geom_point(color = "#00bfc4") +  
  theme_bw() +
  geom_text(vjust = 0.4, color = "black", size = 3.5) + 
  scale_size_continuous(range = c(4, 10)) +  
  theme(axis.text.x = element_text(size = 10, color = 'black', hjust = 1, vjust = 1, angle = 35), 
        axis.text.y = element_text(size = 10, color = 'black'),
        legend.text = element_text(size = 10, color = 'black'),
        legend.title = element_text(size = 10, color = 'black'), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none")+
  xlab("Brain region")+
  ylab("Datasets")
dev.off()



################# Brain region correlation 
Idents(single_data)<-single_data$tissue.region
Cellratio <- prop.table(table(Idents(single_data), single_data$anno), 2)
Cellratio <- as.data.frame(Cellratio)
Cellratio$Var1<-factor(Cellratio$Var1,levels = c("A1C","ACC","CER","CN","CTX","DFC","FC","LA","M1C","MTG","PFC","S1C","SN","V1C"))

ggplot(Cellratio) + 
  geom_bar(aes(x =Var2, y= Freq, fill = Var1),stat = "identity",width = 0.7,size = 0.5,colour = NA)+ 
  scale_fill_manual(values = brain_region_colors)+
  theme_bw()+
  theme(axis.text.x = element_text(size = 10, color = 'black', hjust = 1, vjust = 1, angle = 35), 
        axis.text.y = element_text(size = 10, color = 'black'),
        legend.text = element_text(size = 10, color = 'black'),
        legend.title = element_text(size = 10, color = 'black'), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  ylab("Percentage")+
  xlab("")

#### calculate the correlation matrix
percent_matrix<-matrix(Cellratio$Freq,nrow = table(Cellratio$Var2),ncol = table(Cellratio$Var1))
colnames(percent_matrix)<-names(table(Cellratio$Var2))
rownames(percent_matrix)<-Cellratio$Var1[1:14]
percent_matrix<-t(percent_matrix)
cor_matrix<-cor(percent_matrix,method = "spearman")

pdf("Correlation_for_different_brain_region.pdf",height = 5,width = 5)
corrplot(corr=cor_matrix,method="color",order="hclust",hclust.method="ward.D2",sig.level=0.05,col=COL2('RdBu', 200)[200:1],addrect = 2)
dev.off()


################ celltypes percent in different age period
##区间左闭右开
##P1:4-8  child
##P2:12-20 adolescent
##P3:20-40 young people
##P4:40-60 Middle-aged
##P5:>=60 senium
Idents(single_data)<-single_data$age
single_data<-subset(single_data,idents = "adult",invert=T)
table(single_data$age)

## add another column for age period
single_data$age1<-single_data$age
single_data$age1[WhichCells(single_data,idents = c("4","6","8"))]<-"P1"
single_data$age1[WhichCells(single_data,idents = c("12","13","14","15","19"))]<-"P2"
single_data$age1[WhichCells(single_data,idents = c("20","21","22","23","35","36"))]<-"P3"
single_data$age1[WhichCells(single_data,idents = c("43","48","49","50","54","55","56","59"))]<-"P4"
single_data$age1[WhichCells(single_data,idents = c("60","64","70"))]<-"P5"
temp<-FetchData(single_data,vars=c("anno","sample","age1"))

##比例=mean(Pn阶段某个样本某一细胞类型数目/Pn阶段该样本所有细胞类型数)
## 计算Pn阶段某样本某个细胞类型的数目
dot<-ddply(temp,.(anno,sample),count)

## 计算Pn阶段某个样本提供的全部细胞数目
total_cells<-aggregate(freq ~ sample + age1,data = dot,sum)

## 合并数据框
merged_df <- merge(dot, total_cells, by = c("sample", "age1"))

## 用Pn阶段某样本某个细胞类型的数目/Pn阶段某个样本提供的全部细胞数目
merged_df$proportion <- merged_df$freq.x / merged_df$freq.y

## 取均值操作
mean_proportion <- aggregate(proportion ~ anno + age1, data = merged_df, FUN = mean)

pdf("Celltypes_percent_with_age_period.pdf",height = 4,width=6.5)
ggplot(data=mean_proportion,aes(x=age1,y=proportion,color=anno,group=anno))+
  theme_bw()+
  geom_line(size = 1)+
  scale_color_manual(values = c("#d5231d","#3777ac","#4ea64a","#8e4c99","#e88f18","#e47faf","#b698c5","#a05528"))+
  theme(axis.text.x = element_text(size = 10, color = 'black', hjust = 0.5, vjust = 1, angle = 0), 
        axis.text.y = element_text(size = 10, color = 'black'),
        legend.text = element_text(size = 10, color = 'black'),
        legend.title = element_text(size = 10, color = 'black'))+
  xlab("Period")+
  ylab("Proportion")
dev.off()


#aaa<-dot[dot$anno=="Inhibitory neurons"&dot$age1=="P5",]
#bbb<-dot[dot$age1=="P5"&dot$sample%in%aaa$sample,]
#ccc<-aggregate(bbb$freq,by=list(bbb$sample),sum)
#mean(aaa$freq/ccc$x)

############## For cell subtypes
temp<-FetchData(single_data,vars=c("final_anno","sample","age1"))
dot<-ddply(temp,.(final_anno,sample),count)
total_cells<-aggregate(freq ~ sample + age1,data = dot,sum)
merged_df <- merge(dot, total_cells, by = c("sample", "age1"))
merged_df$proportion <- merged_df$freq.x / merged_df$freq.y
mean_proportion <- aggregate(proportion ~ final_anno + age1, data = merged_df, FUN = mean)
###only plot cell subtypes 
mean_proportion<-mean_proportion[!mean_proportion$final_anno%in%c("Astrocytes","Endothelial cells","Microglia","Oligodendrocytes","OPCs","Purkinje neurons"),]


mean_proportion1<-mean_proportion[mean_proportion$final_anno%in%c("Ex-L2/3","Ex-L2/4","Ex-L4/6","Ex-L5","Ex-L5/6","Ex-L6","Ex-NRGN","Ex_mix"),]


pdf("Excitatory_neurons_subtypes_percent_with_age_period.pdf",height = 4,width=6.5)
ggplot(data=mean_proportion1,aes(x=age1,y=proportion,color=final_anno,group=final_anno))+
  theme_bw()+
  geom_line(size = 1)+
  scale_color_manual(values = c("#fbbab6","#e1c548","#5fa664","#abd0a7","#ca6a6b","#e5b5b5","#f9766e","#bac4d0"))+
  theme(axis.text.x = element_text(size = 10, color = 'black', hjust = 0.5, vjust = 1, angle = 0), 
        axis.text.y = element_text(size = 10, color = 'black'),
        legend.text = element_text(size = 10, color = 'black'),
        legend.title = element_text(size = 10, color = 'black'))+
  xlab("Period")+
  ylab("Proportion")
dev.off()


mean_proportion2<-mean_proportion[mean_proportion$final_anno%in%c("In_CALM1","In_LAMP5","In_PVALB","In_SHANK2","In_SST","In_VIP"),]

pdf("Inhibitory_neurons_subtypes_percent_with_age_period.pdf",height = 4,width=6.5)
ggplot(data=mean_proportion2,aes(x=age1,y=proportion,color=final_anno,group=final_anno))+
  theme_bw()+
  geom_line(size = 1)+
  scale_color_manual(values = c("#8fc0dc","#e98741","#fab37f","#967568","#64abc0","#e0bc58"))+
  theme(axis.text.x = element_text(size = 10, color = 'black', hjust = 0.5, vjust = 1, angle = 0), 
        axis.text.y = element_text(size = 10, color = 'black'),
        legend.text = element_text(size = 10, color = 'black'),
        legend.title = element_text(size = 10, color = 'black'))+
  xlab("Period")+
  ylab("Proportion")
dev.off()

############# density plot for marker genes
###Excitatory neurons
pdf("Density_Plot_SLC17A7_for_Excitatory_neurons.pdf",height = 5,width = 5.2)
Plot_Density_Custom(seurat_object =single_data, features = "SLC17A7",reduction = "umap",
                    custom_palette = c("#B0CFE4","#FACABC","#E77A77","#DC0000FF"))+
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

pdf("Density_Plot_SATB2_for_Excitatory_neurons.pdf",height = 5,width = 5.2)
Plot_Density_Custom(seurat_object =single_data, features = "SATB2",reduction = "umap",
                    custom_palette = c("#B0CFE4","#FACABC","#E77A77","#DC0000FF"))+
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

###Inhibitory neurons
pdf("Density_Plot_GAD1_for_Inhibitory_neurons.pdf",height = 5,width = 5.2)
Plot_Density_Custom(seurat_object =single_data, features = "GAD1",reduction = "umap",
                    custom_palette = c("#B0CFE4","#FACABC","#E77A77","#DC0000FF"))+
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

pdf("Density_Plot_GAD2_for_Inhibitory_neurons.pdf",height = 5,width = 5.2)
Plot_Density_Custom(seurat_object =single_data, features = "GAD2",reduction = "umap",
                    custom_palette = c("#B0CFE4","#FACABC","#E77A77","#DC0000FF"))+
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

###Purkinje neurons
pdf("Density_Plot_RELN_for_Purkinje_neurons.pdf",height = 5,width = 5.2)
Plot_Density_Custom(seurat_object =single_data, features = "RELN",reduction = "umap",
                    custom_palette = c("#B0CFE4","#FACABC","#E77A77","#DC0000FF"))+
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

###Endothelial cells
pdf("Density_Plot_FLT1_for_Endothelial_cells.pdf",height = 5,width = 5.2)
Plot_Density_Custom(seurat_object =single_data, features = "FLT1",reduction = "umap",
                    custom_palette = c("#B0CFE4","#FACABC","#E77A77","#DC0000FF"))+
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

###Astrocytes
pdf("Density_Plot_SLC1A2_for_Astrocytes.pdf",height = 5,width = 5.2)
Plot_Density_Custom(seurat_object =single_data, features = "SLC1A2",reduction = "umap",
                    custom_palette = c("#B0CFE4","#FACABC","#E77A77","#DC0000FF"))+
  theme(plot.title = element_text(hjust = 0.5))
dev.off()


###Oligodendrocytes
pdf("Density_Plot_MOBP_for_Oligodendrocytes.pdf",height = 5,width = 5.2)
Plot_Density_Custom(seurat_object =single_data, features = "MOBP",reduction = "umap",
                    custom_palette = c("#B0CFE4","#FACABC","#E77A77","#DC0000FF"))+
  theme(plot.title = element_text(hjust = 0.5))
dev.off()


###OPCs
pdf("Density_Plot_PCDH15_for_OPCs.pdf",height = 5,width = 5.2)
Plot_Density_Custom(seurat_object =single_data, features = "PCDH15",reduction = "umap",
                    custom_palette = c("#B0CFE4","#FACABC","#E77A77","#DC0000FF"))+
  theme(plot.title = element_text(hjust = 0.5))
dev.off()


###Microglia
pdf("Density_Plot_FYB_for_Microglia.pdf",height = 5,width = 5.2)
Plot_Density_Custom(seurat_object =single_data, features = "FYB",reduction = "umap",
                    custom_palette = c("#B0CFE4","#FACABC","#E77A77","#DC0000FF"))+
  theme(plot.title = element_text(hjust = 0.5))
dev.off()


############# dotplot for broad celltypes
genes <- c("SLC17A7", "SATB2", ## Excitatory neurons
           "GAD1", "GAD2", ## Inhibitory neurons
           "RELN", ## Purkinje neurons
           "FLT1", "DUSP1", "COBLL1", ## Endothelial cells
           "SLC1A2", "SLC1A3", "SLC4A4", ## Astrocytes
           "MOBP", "MBP", "MOG", ## Oligodendrocytes
           "PCDH15", "PDGFRA", ## OPCs
           "APBB1IP", "P2RY12", "FYB" ## Microglia
)

single_data$anno<-factor(single_data$anno,levels = c("Excitatory neurons","Inhibitory neurons","Purkinje neurons","Endothelial cells",
                                                     "Astrocytes","Oligodendrocytes","OPCs","Microglia"))

pdf("Dotplot_for_broad_celltypes.pdf",height = 5,width = 5)
DotPlot(single_data,features=unique(rev(genes)),group.by="anno")+
  scale_colour_gradient2(low="#3A71AA",mid="white",high="#B22028",midpoint=0)+
  theme_bw()+
  coord_flip()+
  theme(axis.text.x = element_text (size = 10,color = 'black',angle = 35,vjust = 1, hjust = 1),
        axis.text.y = element_text(size = 10,color = 'black'), 
        legend.text = element_text(size = 10,color = 'black'), 
        legend.title = element_text(size = 10,color = 'black') 
  )
dev.off()



############# layer specific gene mapping
###Inhibitory neurons
In<-readRDS("D:/Project/SingleCell_MDD/SingleCell_analysis/annotation/Inhibitory_neuron_subtypes_reannotation.rds")
Idents(In)<-In$Inhibitory_neuron_subtypes
data_avr <- AverageExpression(In,slot="data")
data_avr<-as.data.frame(data_avr$SCT)

ln_layer_gene<-data.frame(layer=c("L1/2","L1/2","L1/2/6","L2/3","L3","L4/5","L5/6","L5/6","L6","L6b"),
                          gene=c("CXCL14","CHRNA7","CNR1","LAMP5","SV2C","SULF1","PDE1A","TOX","SYNPR","ADRA2A"))
data_avr1<-data_avr[row.names(data_avr)%in%ln_layer_gene$gene,]
data_avr1<-data_avr1[c("CXCL14","CHRNA7","CNR1","LAMP5","SV2C","SULF1","PDE1A","TOX","SYNPR","ADRA2A"),]
## scale
data_avr1<-t(scale(t(data_avr1)))

## melt the data
data<-melt(data_avr1)
colnames(data)<-c('gene','sample','value')
data$gene<-factor(data$gene,levels = rev(c("CXCL14","CHRNA7","CNR1","LAMP5","SV2C","SULF1","PDE1A","TOX","SYNPR","ADRA2A")))
data_group<-data.frame(genegroup=c("L1/2","L1/2","L1/2/6","L2/3","L3","L4/5","L5/6","L5/6","L6","L6b"),
                       gene=c("CXCL14","CHRNA7","CNR1","LAMP5","SV2C","SULF1","PDE1A","TOX","SYNPR","ADRA2A"))

pdf("Inhibitory_neurons_layer_specific_gene.pdf",width = 6,height = 5)
ggplot(data, aes(sample, gene)) + 
  geom_tile(aes(fill = value))+
  scale_fill_gradient2(low="#67ADB7",high="#af2157") + 
  coord_equal()+
  scale_y_discrete(position = 'right') +
  scale_x_discrete(position = 'top') +
  geom_text(data = data_group,aes(x=0,y=gene,label=genegroup),size=3.9)+
  coord_cartesian(clip = 'off') +
  theme_minimal() +
  theme(axis.title.x=element_blank(), 
        axis.title.y=element_blank(), 
        axis.ticks.x=element_blank(), 
        axis.ticks.y=element_blank(), 
        axis.text.x = element_text(size = 11,color = 'black',angle = 90,hjust = 0),
        axis.text.y = element_text(size = 10,color = 'black',hjust = 1), 
        legend.text = element_text(size = 10,color = 'black'),
        legend.title = element_text(size = 10,color = 'black'),
        panel.grid.major =  element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin=unit(c(2,2,2,8),'lines'),
        legend.position = 'right')
dev.off()


#### Excitatory neurons
EX<-readRDS("D:/Project/SingleCell_MDD/SingleCell_analysis/annotation/Excitatory_neuron_subtypes_reannotation.rds")
Idents(EX)<-EX$Excitatory_neuron_subtypes
data_avr <- AverageExpression(EX,slot="data")
data_avr<-as.data.frame(data_avr$SCT)

ex_layer_gene<-data.frame(layer=c("L2","L2/3/4","L2","L2/3","L3/4","L4","L3/5/6","L5/6","L6","L5/6","L5/6","L5/6","L6","L4c/6","L6","L6","L6/6b","L6b"),
                          gene=c("LAMP5","CUX2","GLRA3","CARTPT","PRSS12","RORB","TPBG","TOX","FOXP2","ETV1","RPRM","RXFP1","TLE4","GRIK4","NTNG2","OPRK1",
                                 "NR4A2","ADRA2A"))
data_avr1<-data_avr[row.names(data_avr)%in%ex_layer_gene$gene,]
data_avr1<-data_avr1[c("LAMP5","CUX2","GLRA3","CARTPT","PRSS12","RORB","TPBG","TOX","FOXP2","ETV1","RPRM","RXFP1","TLE4","GRIK4","NTNG2","OPRK1",
                       "NR4A2","ADRA2A"),]

## scale
data_avr1<-t(scale(t(data_avr1)))
## melt
data<-melt(data_avr1)
colnames(data)<-c('gene','sample','value')
data$gene<-factor(data$gene,levels = rev(c("LAMP5","CUX2","GLRA3","CARTPT","PRSS12","RORB","TPBG","TOX","FOXP2","ETV1",
                                           "RPRM","RXFP1","TLE4","GRIK4","NTNG2","OPRK1","NR4A2","ADRA2A")))

data_group<-data.frame(genegroup=c("L2","L2/3/4","L2","L2/3","L3/4","L4","L3/5/6","L5/6","L6","L5/6","L5/6","L5/6","L6","L4c/6","L6","L6","L6/6b","L6b"),
                       gene=c("LAMP5","CUX2","GLRA3","CARTPT","PRSS12","RORB","TPBG","TOX","FOXP2","ETV1","RPRM","RXFP1",
                              "TLE4","GRIK4","NTNG2","OPRK1","NR4A2","ADRA2A"))


pdf("Excitatory_neurons_layer_specific_gene.pdf",width = 6,height = 5)
ggplot(data, aes(sample, gene)) + 
  geom_tile(aes(fill = value))+
  scale_fill_gradient2(low="#67ADB7",high="#af2157") + 
  coord_equal()+
  scale_y_discrete(position = 'right') +
  scale_x_discrete(position = 'top') +
  geom_text(data = data_group,aes(x=0,y=gene,label=genegroup),size=3.9)+
  coord_cartesian(clip = 'off') +
  theme_minimal() + 
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.ticks.x=element_blank(), 
        axis.ticks.y=element_blank(),
        axis.text.x = element_text(size = 11,color = 'black',angle = 90,hjust = 0),
        axis.text.y = element_text(size = 10,color = 'black'), 
        legend.text = element_text(size = 10,color = 'black'),
        legend.title = element_text(size = 10,color = 'black'),
        panel.grid.major =  element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin=unit(c(2,2,2,10),'lines'),
        legend.position = 'right')
dev.off()










