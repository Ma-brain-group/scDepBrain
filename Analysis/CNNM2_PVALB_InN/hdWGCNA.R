################## scWGCNA ###################
library(Seurat)
library(hdWGCNA)
library(WGCNA)
library(tidyverse)
library(cowplot)
library(patchwork)
library(dplyr)
library(UCell)
library(sctransform)
library(igraph)
theme_set(theme_cowplot())

################## set random seed for reproducibility
set.seed(123)
################## optionally enable multithreading
enableWGCNAThreads(nThreads = 8)

################## set the workspace path
setwd("D:/Project/SingleCell_MDD/SingleCell_analysis/scWGCNA/")

################## import the In6 singlecell data with annotation
single_data<-readRDS("D:/Project/SingleCell_MDD/SingleCell_analysis/11_CNNM2_positive_vs_CNNM2_negnative_IN6/CNNM2_IN6.rds")
DefaultAssay(single_data)<-"RNA"
Idents(single_data)<-single_data$CNNM2_type

################## calculate the DEG between CNNM2+ and CNNM2-
DEG<-FindMarkers(single_data,assay = "RNA",group.by = "CNNM2_type",ident.1 = "CNNM2+",ident.2 = "CNNM2-",only.pos = F)
DEG<-DEG[DEG$p_val_adj<0.05,]

################# SetupForWGCNA
single_data <- SetupForWGCNA(
  single_data,
  gene_select = "custom", # the gene selection approach
  features = row.names(DEG),
  wgcna_name = "tutorial" # the name of the hdWGCNA experiment
)

################ Construct metacells
################ metacells是来自同一来源生物样品的一小群相似细胞的聚集体
single_data <- MetacellsByGroups(
  seurat_obj = single_data,
  group.by = c("CNNM2_type","sample"), # specify the columns in seurat_obj@meta.data to group by
  reduction = 'harmony', # select the dimensionality reduction to perform KNN on
  k = 25, # nearest-neighbors parameter
  max_shared = 10, # maximum number of shared cells between two metacells
  ident.group = "CNNM2_type" # set the Idents of the metacell seurat object
)
################ normalize metacell expression matrix
single_data <- NormalizeMetacells(single_data)


################ Set up the expression matrix
single_data <- SetDatExpr(
  single_data,
  group_name = "CNNM2+", # the name of the group of interest in the group.by column
  group.by='CNNM2_type', # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
  assay = 'RNA', # using RNA assay
  slot = 'data' # using normalized data
)

############# Select soft-power threshold
single_data <- TestSoftPowers(
  single_data,
  networkType = 'signed' # you can also use "unsigned" or "signed hybrid"
)
plot_list <- PlotSoftPowers(single_data)

pdf("soft-power.pdf")
wrap_plots(plot_list, ncol=2)
dev.off()

power_table <- GetPowerTable(single_data)
head(power_table)

############# Construct co-expression network
############# set the softpower
single_data <- ConstructNetwork(
  single_data,
  soft_power = 5,
  tom_name = 'CNNM2+' # name of the topoligical overlap matrix written to disk
)
############# 下游分析中应该忽略灰色module中的基因
PlotDendrogram(single_data, main='CNNM2+ hdWGCNA Dendrogram')

############ inspect the topoligcal overlap matrix (TOM)
TOM <- GetTOM(single_data)


############# Compute harmonized module eigengenes
single_data <- ModuleEigengenes(
  single_data,
  group.by.vars="sample"
)
hMEs <- GetMEs(single_data)


############# Compute module connectivity
single_data <- ModuleConnectivity(
  single_data,
  group.by = 'CNNM2_type', group_name = 'CNNM2+'
)

single_data <- ResetModuleNames(
  single_data,
  new_name = "Depression-M"
)

############# get module gene without grey module
modules <- GetModules(single_data) %>% subset(module != 'grey')
modules[modules$gene_name=="CNNM2",] ##CNNM2 is in module4

############# get the hub gene from each module
hub_df <- GetHubGenes(single_data, n_hubs = 50)

############# save the result
saveRDS(single_data, file='hdWGCNA.rds')



############# 可视化结果 ################
hdWGCNA<-readRDS("hdWGCNA.rds")
MEs <- GetMEs(hdWGCNA, harmonized=TRUE)
modules <- GetModules(hdWGCNA)
mods <- levels(modules$module)
mods <- mods[mods != 'grey']
######### 将ME的信息加入单细胞数据中
hdWGCNA@meta.data <- cbind(hdWGCNA@meta.data, MEs)
data<-FetchData(hdWGCNA,vars = c("CNNM2_type","Depression-M1","Depression-M2","Depression-M3","Depression-M4","Depression-M5",
                                 "Depression-M6","Depression-M7"))
data$CNNM2_type<-factor(data$CNNM2_type,levels = c("CNNM2+","CNNM2-"))
data1<-aggregate(data[,2:8],by=list(data$CNNM2_type),FUN=median)

######### violin plot
pplist = list()
for (i in colnames(data)[2:8]){
data2<-data[,c("CNNM2_type",i)]
colnames(data2)<-c("CNNM2_type","score")
p<-ggplot(data2, aes(x=CNNM2_type, y=score,fill=CNNM2_type)) + 
  geom_violin(trim=FALSE,color="white") + #绘制小提琴图, “color=”设置小提琴图的轮廓线的颜色(以下设为背景为白色，其实表示不要轮廓线)
  #"trim"如果为TRUE(默认值),则将小提琴的尾部修剪到数据范围。如果为FALSE,不修剪尾部。
  geom_boxplot(width=0.1,
               position=position_dodge(0.9),
               color="black",
               outlier.color = "black",
               outlier.fill = "black",
               outlier.shape = 19,
               outlier.size = 0.3,
               outlier.stroke = 0.5,
               outlier.alpha = 45,
               notch = F,
               notchwidth = 0.5)+ #绘制箱线图
  scale_fill_manual(values = c("#E77A77","#89C8E8"))+ #设置填充的颜色
  theme_bw()+ #背景变为白色
  theme(plot.title = element_text(size=12,hjust=0.5), # 标题居中
        axis.text.x = element_text(size = 10,color = 'black'), # 调整x轴坐标文字
        axis.text.y = element_text(size = 10,color = 'black'),
        legend.text = element_text(size = 10,color = 'black'), 
        legend.title = element_text(size = 10,color = 'black'), 
        panel.grid.major = element_blank(),   #不显示网格线
        panel.grid.minor = element_blank(),
        legend.position="none")+  #不显示网格线
  xlab("")+
  ylab(paste0(i,"scores"))
ggsave(filename = paste0(i,"_scores.pdf"),p,height = 4,width = 4)
pplist[[i]] = p
}

pdf("merge1.pdf",width = 6.5,height = 6)
plot_grid(pplist[['Depression-M1']],
          pplist[['Depression-M2']],
          pplist[['Depression-M3']],
          pplist[['Depression-M4']]
)
dev.off()

pdf("merge2.pdf",width = 6.5,height = 6)
plot_grid(pplist[['Depression-M5']],
          pplist[['Depression-M6']],
          pplist[['Depression-M7']]
)
dev.off()


############ dotplot
hdWGCNA@meta.data$CNNM2_type<-factor(hdWGCNA@meta.data$CNNM2_type,levels=c("CNNM2+","CNNM2-"))
pdf("MEs_dotplot.pdf")
DotPlot(hdWGCNA, features=mods, group.by = 'CNNM2_type')+
  scale_colour_gradient2(low="#3A71AA",mid="white",high="#B22028",midpoint=0)+
  theme_bw()+
  coord_flip()+
  theme(axis.text.x = element_text(size = 10,color = 'black',vjust = 1,hjust = 1,angle = 35), 
        axis.text.y = element_text(size = 10,color = 'black'), 
        legend.text = element_text(size = 10,color = 'black'), 
        legend.title = element_text(size = 10,color = 'black'), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin=unit(c(1,1,1,1),'lines'),
        legend.position = 'right')
dev.off()

############ hub gene network
ModuleNetworkPlot(
  hdWGCNA,
  outdir = 'ModuleNetworks',
  vertex.label.cex=0.8
)

ModuleNetworkPlot(
  hdWGCNA, 
  mods = "Depression-M4",
  outdir='ModuleNetworks2', # new folder name
  n_inner = 10, # number of genes in inner ring
  n_outer = 30, # number of genes in outer ring
  n_conns = Inf, # show all of the connections
  plot_size=c(10,10), # larger plotting area
  vertex.label.cex=1 # font size
)

########### 提取出Depression-M4中top100的基因
hub_df <- GetHubGenes(hdWGCNA, n_hubs = 300)
write.csv(hub_df,file = "hub_df_top300.csv",quote = F)
########### Depression-M4中只有132个基因
hub_df_module4 <- hub_df[hub_df$module=="Depression-M4",]
write.csv(hub_df_module4,file = "module4_hub_gene_top132.csv",quote = F)

hub_df_module4 <- hub_df_module4[order(hub_df_module4$kME,decreasing = T),]
hub_df_module4 <- hub_df_module4[1:100,]
write.csv(hub_df_module4,file = "module4_hub_gene_top100.csv",quote = F)

########### 富集分析
library(ggplot2)
library(dplyr)
library(DOSE)
library(org.Hs.eg.db)
library(topGO)
library(clusterProfiler)
library(pathview)
library(enrichplot)

######### KEGG
marker <- bitr(hub_df_module4$gene_name, fromType="SYMBOL",toType="ENTREZID", OrgDb="org.Hs.eg.db")
kk_up <- enrichKEGG(marker$ENTREZID,
                    keyType = "kegg",
                    organism="hsa",
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.05)
kk_up<-setReadable(kk_up, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
temp<-kk_up@result
temp<-temp[temp$p.adjust<0.05,]
write.csv(temp,file = "module4_top100_hub_KEGG.csv",quote = F)

temp$`-log10(fdr)`<- -log10(temp$p.adjust)
pdf("module4_top100_hub_KEGG.pdf",height = 4)
ggplot(data = temp, 
       aes(`-log10(fdr)`, reorder(Description,`-log10(fdr)`))) +
  geom_bar(stat="identity",
           alpha=0.5,
           fill="#FE8D3C",
           color="black",
           width = 0.8) + 
  theme_classic()+
  xlab("-log10(fdr)")+
  ylab("")+
  scale_x_continuous(expand = c(0,0))+
  theme(axis.text.y = element_text(colour = 'black', size = 10),
        axis.line.y = element_line(colour = 'black', linewidth = 0.7),
        axis.title.y = element_text(colour = 'black', size = 12),
        axis.ticks.y = element_line(colour = 'black', linewidth = 0.7),
        axis.line.x = element_line(colour = 'black', linewidth = 0.7),
        axis.text.x = element_text(colour = 'black', size = 10),
        axis.ticks.x = element_line(colour = 'black', linewidth = 0.7),
        axis.title.x = element_text(colour = 'black', size = 12),
        panel.grid = element_blank()
  )
dev.off()



########## GO
ego_up <- enrichGO(marker$ENTREZID,
                   OrgDb = "org.Hs.eg.db",
                   keyType="ENTREZID",
                   pAdjustMethod="BH",
                   ont="BP",
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05,
                   readable = T)

## 对GO terms去除冗余
ego_up <- clusterProfiler::simplify(ego_up, cutoff=0.7, by="p.adjust", select_fun=min)
temp_GO<-ego_up@result
temp_GO<-temp_GO[temp_GO$p.adjust<0.05,]
write.csv(temp_GO,file = "module4_top100_hub_GO.csv",quote = F)


temp_GO1<-temp_GO[c(1:10),]
temp_GO1$`-log10(fdr)`<- -log10(temp_GO1$p.adjust)
pdf("module4_top100_hub_GO.pdf",height = 4)
ggplot(data = temp_GO1, 
       aes(`-log10(fdr)`, reorder(Description,`-log10(fdr)`))) +
  geom_bar(stat="identity",
           alpha=0.5,
           fill="#FE8D3C",
           color="black",
           width = 0.8) + 
  theme_classic()+
  xlab("-log10(fdr)")+
  ylab("")+
  scale_x_continuous(expand = c(0,0))+
  theme(axis.text.y = element_text(colour = 'black', size = 10),
        axis.line.y = element_line(colour = 'black', linewidth = 0.7),
        axis.title.y = element_text(colour = 'black', size = 12),
        axis.ticks.y = element_line(colour = 'black', linewidth = 0.7),
        axis.line.x = element_line(colour = 'black', linewidth = 0.7),
        axis.text.x = element_text(colour = 'black', size = 10),
        axis.ticks.x = element_line(colour = 'black', linewidth = 0.7),
        axis.title.x = element_text(colour = 'black', size = 12),
        panel.grid = element_blank()
)
dev.off()


############## 将module4中的基因映射到PPI网络中 #################
string<-read.table("D:/Project/SingleCell_MDD/SingleCell_analysis/Manhattan_Plot/9606.protein.physical.links.v12.0.txt",header = T)
############## 转换网络中的蛋白名
string_annotation<-read.table("D:/Project/SingleCell_MDD/SingleCell_analysis/Manhattan_Plot/9606.protein.info.v12.0.txt",header = T,sep = "\t")
string_annotation<-string_annotation[,1:2]


colnames(string_annotation)[1]<-"protein1"
string<-merge(string_annotation,string,by="protein1")
string<-string[,2:4]
colnames(string)[1]<-"source"

colnames(string_annotation)[1]<-"protein2"
string<-merge(string_annotation,string,by="protein2")
string<-string[,c(3,2,4)]
colnames(string)[2]<-"target"



############## 映射基因
module4_hub_gene<-read.csv("module4_hub_gene_top100.csv",row.names = 1)
colnames(module4_hub_gene)[1]<-"source"
############## 只需要匹配一列即可
string1<-merge(string,module4_hub_gene,by="source")
### 2672行
string1<-string1[,1:3]

string2 <- string1 %>%
  group_by(source) %>%
  arrange(desc(combined_score)) %>%
  slice_head(n = 10) %>%
  ungroup()

### 若source和target相同，则去掉其中一行
string2 <- string2 %>%
  mutate(identifier = ifelse(source < target, paste(source, target, sep = ","), paste(target, source, sep = ",")))

# 去掉重复的互换行
string3 <- string2 %>%
  distinct(identifier, .keep_all = TRUE) %>%
  select(-identifier)  # 删除标识符列
write.csv(string3,file = "network2.csv",quote = F)




############# 后续不需要，后续为hdWGCNA的绘制网络图的代码
############# hdWGCNA hub gene plot 
############# get the hub gene neywork 
wgcna_name <- hdWGCNA@misc$active_wgcna
MEs <- GetMEs(hdWGCNA, wgcna_name)
modules <- GetModules(hdWGCNA, wgcna_name)
TOM <- GetTOM(hdWGCNA, wgcna_name)

############# set the module and hub gene number
mods <- "Depression-M4"
n_hubs <- 100
hub_list <- lapply(mods, function(cur_mod) {
  cur <- subset(modules, module == cur_mod)
  cur <- cur[, c("gene_name", paste0("kME_", cur_mod))] %>% 
    top_n(n_hubs)
  colnames(cur)[2] <- "var"
  cur %>% arrange(desc(var)) %>% .$gene_name
})

names(hub_list) <- mods

cur_mod <- mods
print(cur_mod)
cur_color <- modules %>% subset(module == cur_mod) %>% 
  .$color %>% unique
n_genes = 100

########### set the edges number
n_conns = 500
cur_kME <- paste0("kME_", cur_mod)
cur_genes <- hub_list[[cur_mod]]
matchind <- match(cur_genes, colnames(TOM))
reducedTOM = TOM[matchind, matchind]
########### choose top 1500 gene for edges
orderind <- order(reducedTOM, decreasing = TRUE)
connections2keep <- orderind[1:n_conns]
reducedTOM <- matrix(0, nrow(reducedTOM), ncol(reducedTOM))

########### get the adjacency matrix
reducedTOM[connections2keep] <- 1
colnames(reducedTOM)<-hub_list$`Depression-M4`
rownames(reducedTOM)<-hub_list$`Depression-M4`


########### change the format of adjacency matrix
result <- data.frame(Row = character(), Column = character(), Value = numeric(), stringsAsFactors = FALSE)

for (i in 1:nrow(reducedTOM)) {
  for (j in 1:ncol(reducedTOM)) {
    result <- rbind(result, data.frame(Row = rownames(reducedTOM)[i], 
                                       Column = colnames(reducedTOM)[j], 
                                       Value = reducedTOM[i, j]))
  }
}

########### remove the dege is 0
result<-result[result$Value!=0,]
table(result$Row)

########### 创建一个辅助列，用于标记对称行
result$Pair <- ifelse(result$Row < result$Column, 
                      paste(result$Row, result$Column, sep = "-"), 
                      paste(result$Column, result$Row, sep = "-"))
########### 去除对称重复行
result1 <- result[!duplicated(result$Pair), ]

########### 去掉辅助列
result1$Pair <- NULL

########### save the network
colnames(result1)<-c("source","target","edge")
write.csv(result1,file = "network.csv",quote = F)




label_center = FALSE
if (label_center) {
  cur_genes[11:25] <- ""
}
gA <- graph.adjacency(as.matrix(reducedTOM[1:10, 1:10]), 
                      mode = "undirected", weighted = TRUE, diag = FALSE)
gB <- graph.adjacency(as.matrix(reducedTOM[11:n_genes, 
                                           11:n_genes]), mode = "undirected", weighted = TRUE, 
                      diag = FALSE)
layoutCircle <- rbind(layout.circle(gA)/2, layout.circle(gB))
g1 <- graph.adjacency(as.matrix(reducedTOM), mode = "undirected", 
                      weighted = TRUE, diag = FALSE)



edge.alpha = 0.25
vertex.label.cex = 1
vertex.size = 6
plot(g1, edge.color = adjustcolor(cur_color, alpha.f = 0.25), 
     edge.alpha = edge.alpha, vertex.color = cur_color, 
     vertex.label = as.character(cur_genes), vertex.label.dist = 1.1, 
     vertex.label.degree = -pi/4, vertex.label.color = "black", 
     vertex.label.family = "Helvetica", vertex.label.font = 3, 
     vertex.label.cex = vertex.label.cex, vertex.frame.color = "black", 
     layout = jitter(layoutCircle), vertex.size = vertex.size, 
     main = paste(cur_mod))









