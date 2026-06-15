################ In_PVALB cellchat ####################
##devtools::install_github("sqjin/CellChat")
library(Seurat)
library(dplyr)
library(SeuratData)
library(patchwork)
library(ggplot2)
library(CellChat)
library(ggalluvial)
library(svglite)
library(ggpubr)
options(stringsAsFactors = FALSE)
setwd("/share2/pub/chenchg/chenchg/SingleCell/Brain/SingleCell_data/In6_cellchat_CNNM2_negnative/")

############ import the singlecell data
single_data_CNNM2_IN6<-readRDS("/share2/pub/chenchg/chenchg/SingleCell/Brain/SingleCell_data/singlecell_data/CNNM2_IN6.rds")
single_data<-readRDS("/share2/pub/chenchg/chenchg/SingleCell/Brain/SingleCell_data/singlecell_data/resolution0.2_final4.rds")

########### add another colunm for CNNM2_type
single_data$CNNM2_type<-single_data$anno
single_data$CNNM2_type[WhichCells(single_data_CNNM2_IN6,idents = c("CNNM2+"))] <- "CNNM2+_In6"
single_data$CNNM2_type[WhichCells(single_data_CNNM2_IN6,idents = c("CNNM2-"))] <- "CNNM2-_In6"
########### define the Inhibitory neurons to other Inhibitory neurons
Idents(single_data)<-single_data$CNNM2_type
single_data$CNNM2_type[WhichCells(single_data,idents = c("Inhibitory neurons"))] <- "Other Inhibitory neurons"

############ create cellchat object
cellchat <- createCellChat(object = single_data, group.by = "CNNM2_type", assay = "RNA")
############# add CellChat ligand-receptor database
CellChatDB <- CellChatDB.human 
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") ##这里也可以不选，直接用默认的也行
cellchat@DB <- CellChatDB.use

############ data preprocess
cellchat <- subsetData(cellchat) # subset the expression data of signaling genes for saving computation cost
#future::plan("multiprocess", workers = 4) # do parallel  
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)  

############# infer the cell-cell interaction
cellchat <- computeCommunProb(cellchat) 
cellchat <- filterCommunication(cellchat, min.cells = 10)

############# infer the cell-cell interaction in signal pathway
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

############# save the cellchat result 
saveRDS(cellchat,file = "cellchat_In6_CNNM2.rds")


############# visualize the result ###################
############# import the cellchat result
setwd("D:/Project/SingleCell_MDD/SingleCell_analysis/In6_CNNM2_cellchat/")
cellchat<-readRDS("cellchat_In6_CNNM2.rds")

############# visualize the cell-cell interaction count and weight
groupSize <- as.numeric(table(cellchat@idents))
mat <- cellchat@net$count
pdf("cellchat_count.pdf",height = 6,width = 6)
#par(mfrow =c(5,2),xpd=T)
for (i in 1:nrow(mat)){
  mat2 <- matrix(0,nrow = nrow(mat),ncol = ncol(mat),dimnames = dimnames(mat))
  mat2[i,] <- mat[i,]
  netVisual_circle(mat2,vertex.weight = groupSize, weight.scale = T, arrow.width = 0.2,
                   arrow.size = 0.1, edge.weight.max = max(mat),title.name = rownames(mat)[i])
}
dev.off()

mat <- cellchat@net$weight
pdf("cellchat_weight.pdf",height = 6,width = 6)
#par(mfrow =c(5,2),xpd=T)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, arrow.width = 0.2,
                   arrow.size = 0.1, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
dev.off()

################ visualize the cell-cell interaction ligand-receptor pair 
################ get the ligand-receptor pair information
df.net <- subsetCommunication(cellchat)
write.csv(df.net, "net_lr.csv")
df.netp <- subsetCommunication(cellchat, slot.name = "netP")
write.csv(df.netp, "net_pathway.csv")



################ calculate the cell-cell interaction numbers between CNNM2+/CNNM2-
cell_interaction_number<-netVisual_bubble(cellchat, remove.isolate = F,return.data = T)$communication
cell_interaction_number1<-cell_interaction_number[cell_interaction_number$source%in%"CNNM2+_In6"|cell_interaction_number$target%in%"CNNM2+_In6",]
cell_interaction_number2<-cell_interaction_number[cell_interaction_number$source%in%"CNNM2-_In6"|cell_interaction_number$target%in%"CNNM2-_In6",]
cell_interaction_number1$source<-as.character(cell_interaction_number1$source)
cell_interaction_number1$target<-as.character(cell_interaction_number1$target)
cell_interaction_number2$source<-as.character(cell_interaction_number2$source)
cell_interaction_number2$target<-as.character(cell_interaction_number2$target)

#### calculate the cell-cell interaction numbers about CNNM2+/
cell_interaction_number1 <- cell_interaction_number1 %>%
  mutate(temp = ifelse(target == "CNNM2+_In6", source, target),
         source = ifelse(target == "CNNM2+_In6", target, source),
         target = ifelse(target == "CNNM2+_In6", temp, target)) %>%
  select(-temp)
table(cell_interaction_number1$target)

cell_interaction_number2 <- cell_interaction_number2 %>%
  mutate(temp = ifelse(target == "CNNM2-_In6", source, target),
         source = ifelse(target == "CNNM2-_In6", target, source),
         target = ifelse(target == "CNNM2-_In6", temp, target)) %>%
  select(-temp)
table(cell_interaction_number2$target)

#### merge the data 
table1 <- as.data.frame(table(cell_interaction_number1$target))
table2 <- as.data.frame(table(cell_interaction_number2$target))
colnames(table1) <- c("Category", "Counts")
colnames(table2) <- c("Category", "Counts")
merge_data<-rbind(table1,table2)
merge_data$type<-c(rep("CNNM2+",10),rep("CNNM2-",9))
merge_data$type<-factor(merge_data$type,levels = c("CNNM2+","CNNM2-"))


groups<-list(c("CNNM2+","CNNM2-"))
pdf("boxplot_for_the_number_of_interactions.pdf",height = 4.5,width = 5.5)
ggplot(merge_data, aes(x=type, y=Counts, color=type, fill=type)) +
  stat_boxplot(geom = "errorbar", width=0.1) +
  geom_boxplot(alpha = 0.5, size=1.5, width = 0.6,outlier.size = 0,outlier.color = "white",outlier.fill = "white") + 
  geom_jitter(size = 2, width = 0.1) +  
  scale_color_manual(values = c("#e97371","#5ac6e9")) +  
  scale_fill_manual(values = c("#e97371","#5ac6e9")) +
  stat_compare_means(comparisons = groups, method = "t.test",label = "p.signif",size=4.5) +
  theme_bw()+  
  theme(plot.title = element_text(size=12,hjust=0.5), 
        axis.text.x = element_text(size = 11,color = 'black',hjust = 0.5,angle = 0), 
        axis.text.y = element_text(size = 11,color = 'black'),
        legend.text = element_text(size = 11,color = 'black'), 
        legend.title = element_text(size = 11,color = 'black'), 
        panel.grid.major = element_blank(),   
        panel.grid.minor = element_blank())+
  ylab("Number of interactions")
dev.off() 
 


################ plotting default function netVisual_bubble   
levels(cellchat@idents)
pdf("ligand_receptor.pdf")
netVisual_bubble(cellchat, sources.use = c(2,3),targets.use = c(1,4:10), remove.isolate = F)+
  scale_colour_gradient(low="white",high="#B22028")+
  theme(axis.text.x = element_text(size = 10,color = 'black',vjust = 1,angle = 35), 
        axis.text.y = element_text(size = 10,color = 'black'), 
        legend.text = element_text(size = 10,color = 'black'), 
        legend.title = element_text(size = 10,color = 'black'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin=unit(c(1,1,1,1),'lines'),
        legend.position = 'right')
dev.off()


########### get the data for ggplot2
temp<-netVisual_bubble(cellchat, sources.use = c(2,3),targets.use = c(1,4:10), remove.isolate = F,return.data = T)$communication
temp<-na.omit(temp)

########## pval>0.05:1，0.01<pval<=0.05:2，pval<=0.01:3
df.net_positive<-temp[temp$source%in%c("CNNM2+_In6"),]
df.net_negnative<-temp[temp$source%in%c("CNNM2-_In6"),]
df.net_negnative<-df.net_negnative[,c("target","interaction_name_2","prob","pval")]
df.net_positive<-df.net_positive[,c("target","interaction_name_2","prob","pval")]


########## CNNM2+_In_PVALB plot
df.net_positive1<-data.frame(target=rep(names(table(temp$target)),26),
                             interaction_name_2=rep(names(table(temp$interaction_name_2)),each = 8))
########## merge the data
df.net_positive2<-merge(df.net_positive,df.net_positive1,all=TRUE)
df.net_positive2$prob[is.na(df.net_positive2$prob)]<-0
df.net_positive2$pval[is.na(df.net_positive2$pval)]<-1

df.net_positive2$pval[df.net_positive2$pval=="1"]<-"p > 0.05"
df.net_positive2$pval[df.net_positive2$pval=="2"]<-"0.01 < p < 0.05"
df.net_positive2$pval[df.net_positive2$pval=="3"]<-"p < 0.01"
df.net_positive2$pval<-factor(df.net_positive2$pval,levels = c("p > 0.05", "0.01 < p < 0.05", "p < 0.01"))
df.net_positive2$title <- "CNNM2+ In_PVALB"

pdf("CNNM2+_In_PVALB_cellchat.pdf",height = 5.5,width = 6)
ggplot(data = df.net_positive2,aes(x=target,y=interaction_name_2,color = prob,size=pval))+
  geom_point()+
  theme_bw()+
  scale_color_distiller(palette = "Spectral")+
  scale_size_manual(values = c("0.01 < p < 0.05" = 2, "p < 0.01" = 4)) +
  theme(axis.text.x = element_text(size = 10,color = 'black',hjust = 1,vjust = 1,angle = 35), 
        axis.text.y = element_text(size = 10,color = 'black'),
        legend.text = element_text(size = 10,color = 'black'),
        legend.title = element_text(size = 10,color = 'black'), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = 'right')+
  facet_grid(. ~ title)
dev.off()


########## CNNM2-_In_PVALB plot
df.net_negnative1<-data.frame(target=rep(names(table(temp$target)),26),
                              interaction_name_2=rep(names(table(temp$interaction_name_2)),each = 8))
########## merge the data
df.net_negnative2<-merge(df.net_negnative,df.net_negnative1,all=TRUE)
df.net_negnative2$prob[is.na(df.net_negnative2$prob)]<-0
df.net_negnative2$pval[is.na(df.net_negnative2$pval)]<-1

df.net_negnative2$pval[df.net_negnative2$pval=="1"]<-"p > 0.05"
df.net_negnative2$pval[df.net_negnative2$pval=="2"]<-"0.01 < p < 0.05"
df.net_negnative2$pval[df.net_negnative2$pval=="3"]<-"p < 0.01"
df.net_negnative2$pval<-factor(df.net_negnative2$pval,levels = c("p > 0.05", "0.01 < p < 0.05", "p < 0.01"))
df.net_negnative2$title <- "CNNM2- In_PVALB"

pdf("CNNM2-_In_PVALB_cellchat.pdf",height = 5.5,width = 6)
ggplot(data = df.net_negnative2,aes(x=target,y=interaction_name_2,color = prob,size=pval))+
  geom_point()+
  theme_bw()+
  scale_color_distiller(palette = "Spectral")+
  scale_size_manual(values = c("0.01 < p < 0.05" = 2, "p < 0.01" = 4)) +
  theme(axis.text.x = element_text(size = 10,color = 'black',hjust = 1,vjust = 1,angle = 35), 
        axis.text.y = element_text(size = 10,color = 'black'),
        legend.text = element_text(size = 10,color = 'black'),
        legend.title = element_text(size = 10,color = 'black'), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = 'right')+
  facet_grid(. ~ title)
dev.off()

