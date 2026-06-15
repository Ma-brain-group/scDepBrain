#################### Ex-L2/4 cellchat analysis ###################
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
setwd("D:/Project/SingleCell_MDD/SingleCell_analysis/Exitatory_neurons_subtype_analysis/cellchat/")

############ import the singlecell data
single_data_Ex_L2_4<-readRDS("D:/Project/SingleCell_MDD/SingleCell_analysis/Exitatory_neurons_subtype_analysis/Ex_L2_4_NEGR_analysis/NEGR1_EX_L2_4.rds")
single_data<-readRDS("D:/Project/SingleCell_MDD/SingleCell_analysis/annotation/MDD_singlecell_data_reannotation_simple.rds")

########### add another colunm for NEGR1_type
single_data$NEGR1_type<-single_data$anno
Idents(single_data_Ex_L2_4)<-single_data_Ex_L2_4$NEGR1_type
single_data$NEGR1_type[WhichCells(single_data_Ex_L2_4,idents = c("NEGR1+"))] <- "NEGR1+_Ex_L2_4"
single_data$NEGR1_type[WhichCells(single_data_Ex_L2_4,idents = c("NEGR1-"))] <- "NEGR1-_Ex_L2_4"
########### define the Excitatory neurons for other Excitatory neurons
Idents(single_data)<-single_data$NEGR1_type
single_data$NEGR1_type<-as.character(single_data$NEGR1_type)
single_data$NEGR1_type[WhichCells(single_data,idents = c("Excitatory neurons"))] <- "Other Excitatory neurons"


############ create cellchat object
cellchat <- createCellChat(object = single_data, group.by = "NEGR1_type", assay = "RNA")

############# add CellChat ligand-receptor database
CellChatDB <- CellChatDB.human 
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
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
saveRDS(cellchat,file = "cellchat_NEGR1_Ex_L2_4.rds")



############# visualize the result ###################
############# import the cellchat result
cellchat<-readRDS("cellchat_NEGR1_Ex_L2_4.rds")


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


################ calculate the cell-cell interaction numbers between NEGR1+/NEGR1-
cell_interaction_number<-netVisual_bubble(cellchat, remove.isolate = F,return.data = T)$communication
cell_interaction_number1<-cell_interaction_number[cell_interaction_number$source%in%"NEGR1+_Ex_L2_4"|cell_interaction_number$target%in%"NEGR1+_Ex_L2_4",]
cell_interaction_number2<-cell_interaction_number[cell_interaction_number$source%in%"NEGR1-_Ex_L2_4"|cell_interaction_number$target%in%"NEGR1-_Ex_L2_4",]
cell_interaction_number1$source<-as.character(cell_interaction_number1$source)
cell_interaction_number1$target<-as.character(cell_interaction_number1$target)
cell_interaction_number2$source<-as.character(cell_interaction_number2$source)
cell_interaction_number2$target<-as.character(cell_interaction_number2$target)

#### calculate the cell-cell interaction numbers about NEGR1+/NEGR1-
cell_interaction_number1 <- cell_interaction_number1 %>%
  mutate(temp = ifelse(target == "NEGR1+_Ex_L2_4", source, target),
         source = ifelse(target == "NEGR1+_Ex_L2_4", target, source),
         target = ifelse(target == "NEGR1+_Ex_L2_4", temp, target)) %>%
  select(-temp)
table(cell_interaction_number1$target)

cell_interaction_number2 <- cell_interaction_number2 %>%
  mutate(temp = ifelse(target == "NEGR1-_Ex_L2_4", source, target),
         source = ifelse(target == "NEGR1-_Ex_L2_4", target, source),
         target = ifelse(target == "NEGR1-_Ex_L2_4", temp, target)) %>%
  select(-temp)
table(cell_interaction_number2$target)

#### merge the data 
table1 <- as.data.frame(table(cell_interaction_number1$target))
table2 <- as.data.frame(table(cell_interaction_number2$target))
colnames(table1) <- c("Category", "Counts")
colnames(table2) <- c("Category", "Counts")
merge_data<-rbind(table1,table2)
merge_data$type<-c(rep("NEGR1+",10),rep("NEGR1-",9))
merge_data$type<-factor(merge_data$type,levels = c("NEGR1+","NEGR1-"))


groups<-list(c("NEGR1+","NEGR1-"))
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


############ plotting default function netVisual_bubble 
levels(cellchat@idents)
pdf("ligand_receptor.pdf")
netVisual_bubble(cellchat, sources.use = c(5,6),targets.use = c(1:4,7:10), remove.isolate = F)+
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
temp<-netVisual_bubble(cellchat, sources.use = c(5,6),targets.use = c(1:4,7:10), remove.isolate = F,return.data = T)$communication
temp<-na.omit(temp)

########## pval>0.05:1，0.01<pval<=0.05:2，pval<=0.01:3
df.net_positive<-temp[temp$source%in%c("NEGR1+_Ex_L2_4"),]
df.net_negnative<-temp[temp$source%in%c("NEGR1-_Ex_L2_4"),]
df.net_negnative<-df.net_negnative[,c("target","interaction_name_2","prob","pval")]
df.net_positive<-df.net_positive[,c("target","interaction_name_2","prob","pval")]


########## NEGR1+_Ex_L2_4 plot
df.net_positive1<-data.frame(target=rep(names(table(temp$target)),38),
                             interaction_name_2=rep(names(table(temp$interaction_name_2)),each = 8))
########## merge the data
df.net_positive2<-merge(df.net_positive,df.net_positive1,all=TRUE)
df.net_positive2$prob[is.na(df.net_positive2$prob)]<-0
df.net_positive2$pval[is.na(df.net_positive2$pval)]<-1

df.net_positive2$pval[df.net_positive2$pval=="1"]<-"p > 0.05"
df.net_positive2$pval[df.net_positive2$pval=="2"]<-"0.01 < p < 0.05"
df.net_positive2$pval[df.net_positive2$pval=="3"]<-"p < 0.01"
df.net_positive2$pval<-factor(df.net_positive2$pval,levels = c("p > 0.05", "0.01 < p < 0.05", "p < 0.01"))
df.net_positive2$title <- "NEGR1+ Ex-L2/4"

pdf("NEGR1+_Ex_L2_4_cellchat.pdf",height = 6.5,width = 6.5)
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


########## NEGR1-_Ex_L2_4 plot
df.net_negnative1<-data.frame(target=rep(names(table(temp$target)),38),
                              interaction_name_2=rep(names(table(temp$interaction_name_2)),each = 8))
########## merge the data
df.net_negnative2<-merge(df.net_negnative,df.net_negnative1,all=TRUE)
df.net_negnative2$prob[is.na(df.net_negnative2$prob)]<-0
df.net_negnative2$pval[is.na(df.net_negnative2$pval)]<-1

df.net_negnative2$pval[df.net_negnative2$pval=="1"]<-"p > 0.05"
df.net_negnative2$pval[df.net_negnative2$pval=="2"]<-"0.01 < p < 0.05"
df.net_negnative2$pval[df.net_negnative2$pval=="3"]<-"p < 0.01"
df.net_negnative2$pval<-factor(df.net_negnative2$pval,levels = c("p > 0.05", "0.01 < p < 0.05", "p < 0.01"))
df.net_negnative2$title <- "NEGR1- Ex-L2/4"

pdf("NEGR1-_Ex_L2_4_cellchat.pdf",height = 6.5,width = 6.5)
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







