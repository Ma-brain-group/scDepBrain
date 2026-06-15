############# NC Ex-L2/4 cellchat ##################
library(Seurat)
library(dplyr)
library(SeuratData)
library(patchwork)
library(ggplot2)
library(CellChat)
library(ggalluvial)
library(svglite)
options(stringsAsFactors = FALSE)
setwd("/share2/pub/chenchg/chenchg/SingleCell/Brain/SingleCell_data/Exitatory_neurons_subtype_analysis/NC_cellchat/")

############ import the singlecell data
single_data_Ex_L2_4<-readRDS("D:/Project/SingleCell_MDD/SingleCell_analysis/Exitatory_neurons_subtype_analysis/NC_Excitatory_neurons_annotation/NC_NEGR1_EX_L2_4.rds")
single_data<-readRDS("D:/Project/SingleCell_MDD/SingleCell_analysis/16_validation_NC/1_NC_data_process/GSE213982/2023_NC.rds")

########### add another colunm for NEGR1_type
single_data$NEGR1_type<-single_data$broad_cell_type
Idents(single_data)<-single_data$broad_cell_type
single_data$NEGR1_type[WhichCells(single_data,idents = c("Ast"))] <- "Astrocytes"
single_data$NEGR1_type[WhichCells(single_data,idents = c("End"))] <- "Endothelial cells"
single_data$NEGR1_type[WhichCells(single_data,idents = c("ExN"))] <- "Excitatory neurons"
single_data$NEGR1_type[WhichCells(single_data,idents = c("Mic"))] <- "Microglia"
single_data$NEGR1_type[WhichCells(single_data,idents = c("Oli"))] <- "Oligodendrocytes"
single_data$NEGR1_type[WhichCells(single_data,idents = c("OPC"))] <- "OPCs"
single_data$NEGR1_type[WhichCells(single_data,idents = c("InN"))] <- "Inhibitory neurons"
########## remove the Mix
single_data<-subset(single_data,idents="Mix",invert=T)

Idents(single_data_Ex_L2_4)<-single_data_Ex_L2_4$NEGR1_type
single_data$NEGR1_type[WhichCells(single_data_Ex_L2_4,idents = c("NEGR1+"))] <- "NEGR1+_Ex-L2/4"
single_data$NEGR1_type[WhichCells(single_data_Ex_L2_4,idents = c("NEGR1-"))] <- "NEGR1-_Ex-L2/4"

########### define the Excitatory neurons to other Excitatory neurons
Idents(single_data)<-single_data$NEGR1_type
single_data$NEGR1_type[WhichCells(single_data,idents = c("Excitatory neurons"))] <- "Other Excitatory neurons"


########### cellchat for different phenotype
Idents(single_data)<-single_data$phenotype
single_data_Case<-subset(single_data,idents="Case")
single_data_Control<-subset(single_data,idents="Control")

############ create cellchat object for different phenotype
cellchat_Case <- createCellChat(object = single_data_Case, group.by = "NEGR1_type", assay = "RNA")
cellchat_Control<-createCellChat(object = single_data_Control, group.by = "NEGR1_type", assay = "RNA")


############# Case group
CellChatDB <- CellChatDB.human 
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") 
cellchat_Case@DB <- CellChatDB.use
cellchat_Case <- subsetData(cellchat_Case) 
cellchat_Case <- identifyOverExpressedGenes(cellchat_Case)
cellchat_Case <- identifyOverExpressedInteractions(cellchat_Case)
cellchat_Case <- projectData(cellchat_Case, PPI.human)  
cellchat_Case <- computeCommunProb(cellchat_Case) 
cellchat_Case <- filterCommunication(cellchat_Case, min.cells = 10)
cellchat_Case <- computeCommunProbPathway(cellchat_Case)
cellchat_Case <- aggregateNet(cellchat_Case)
cellchat_Case <- netAnalysis_computeCentrality(cellchat_Case, slot.name = "netP")
saveRDS(cellchat_Case,file = "cellchat_Case.rds")


########### Control group
CellChatDB <- CellChatDB.human 
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") 
cellchat_Control@DB <- CellChatDB.use
cellchat_Control <- subsetData(cellchat_Control)
cellchat_Control <- identifyOverExpressedGenes(cellchat_Control)
cellchat_Control <- identifyOverExpressedInteractions(cellchat_Control)
cellchat_Control <- projectData(cellchat_Control, PPI.human)  
cellchat_Control <- computeCommunProb(cellchat_Control) 
cellchat_Control <- filterCommunication(cellchat_Control, min.cells = 10)
cellchat_Control <- computeCommunProbPathway(cellchat_Control)
cellchat_Control <- aggregateNet(cellchat_Control)
cellchat_Control <- netAnalysis_computeCentrality(cellchat_Control, slot.name = "netP")
saveRDS(cellchat_Control,file = "cellchat_Control.rds")


############ import the cellchat result
setwd("D:/Project/SingleCell_MDD/SingleCell_analysis/Exitatory_neurons_subtype_analysis/NC_cellchat/")
cellchat_Case<-readRDS("cellchat_Case.rds")
cellchat_Control<-readRDS("cellchat_Control.rds")

############ merge the case and control cellchat result
list <- list(Case=cellchat_Case, Control=cellchat_Control)
cellchat <- mergeCellChat(list, add.names = names(list), cell.prefix = TRUE)

############ get the cell-cell interaction ligand-receptor pair infoemation
df.net <- subsetCommunication(cellchat_Case)
write.csv(df.net, "Case_net_lr.csv")
df.netp <- subsetCommunication(cellchat_Case, slot.name = "netP")
write.csv(df.netp, "Case_net_pathway.csv")

df.net <- subsetCommunication(cellchat_Control)
write.csv(df.net, "Control_net_lr.csv")
df.netp <- subsetCommunication(cellchat_Control, slot.name = "netP")
write.csv(df.netp, "Control_net_pathway.csv")

############# visualize1: the cell-cell interaction count and weight
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "count")
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
p <- gg1 + gg2

############# visualize2: bubble plot
levels(cellchat@idents$joint)
netVisual_bubble(cellchat, sources.use = c(5,6), targets.use = c(1:4,7:9), comparison = c(1, 2), angle.x = 45)



############# get the data for ggplot2
temp<-netVisual_bubble(cellchat_Case, sources.use = c(5,6), targets.use = c(1:4,7:9), remove.isolate = F,return.data = T)$communication
temp<-na.omit(temp)
df.net_positive<-temp[temp$source%in%c("NEGR1+_Ex-L2/4"),]
df.net_negnative<-temp[temp$source%in%c("NEGR1-_Ex-L2/4"),]

############ case NEGR1-_Ex-L2/4
df.net_negnative<-df.net_negnative[,c("target","interaction_name_2","prob","pval")]

df.net_negnative1<-data.frame(target=rep(names(table(temp$target)),18),
                              interaction_name_2=rep(names(table(temp$interaction_name_2)),each = 7))
########### merge
df.net_negnative2<-merge(df.net_negnative,df.net_negnative1,all=TRUE)
df.net_negnative2$prob[is.na(df.net_negnative2$prob)]<-0
df.net_negnative2$pval[is.na(df.net_negnative2$pval)]<-1

df.net_negnative2$pval[df.net_negnative2$pval=="1"]<-"p > 0.05"
df.net_negnative2$pval[df.net_negnative2$pval=="2"]<-"0.01 < p < 0.05"
df.net_negnative2$pval[df.net_negnative2$pval=="3"]<-"p < 0.01"
df.net_negnative2$pval<-factor(df.net_negnative2$pval,levels = c("p > 0.05", "0.01 < p < 0.05", "p < 0.01"))
df.net_negnative2$title <- "Case NEGR1- Ex-L2/4"

pdf("NC_Case_NEGR1-_Ex-L2_4_cellchat.pdf",height = 5,width = 5.5)
ggplot(df.net_negnative2,aes(x=target,y=interaction_name_2,color = prob,size = pval))+
  geom_point()+
  theme_bw()+
  scale_color_distiller(palette = "Spectral")+
  scale_size_manual(values = c("0.01 < p < 0.05" = 2, "p < 0.01" = 4)) +
  theme(axis.text.x = element_text(size = 10,color = 'black',vjust = 1,hjust = 1,angle = 35), 
        axis.text.y = element_text(size = 10,color = 'black'), 
        legend.text = element_text(size = 10,color = 'black'), 
        legend.title = element_text(size = 10,color = 'black'), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin=unit(c(1,1,1,1),'lines'),
        legend.position = 'right')+
  facet_grid(. ~ title)
dev.off()  


########### case NEGR1+_Ex-L2/4
df.net_positive<-df.net_positive[,c("target","interaction_name_2","prob","pval")]
df.net_positive1<-data.frame(target=rep(names(table(temp$target)),18),
                             interaction_name_2=rep(names(table(temp$interaction_name_2)),each = 7))
df.net_positive2<-merge(df.net_positive,df.net_positive1,all=TRUE)
df.net_positive2$prob[is.na(df.net_positive2$prob)]<-0
df.net_positive2$pval[is.na(df.net_positive2$pval)]<-1

df.net_positive2$pval[df.net_positive2$pval=="1"]<-"p > 0.05"
df.net_positive2$pval[df.net_positive2$pval=="2"]<-"0.01 < p < 0.05"
df.net_positive2$pval[df.net_positive2$pval=="3"]<-"p < 0.01"
df.net_positive2$pval<-factor(df.net_positive2$pval,levels = c("p > 0.05", "0.01 < p < 0.05", "p < 0.01"))
df.net_positive2$title<-"Case NEGR1+ Ex-L2/4"

pdf("NC_Case_NEGR1+_Ex-L2_4_cellchat.pdf",height = 5,width = 6)
ggplot(df.net_positive2,aes(x=target,y=interaction_name_2,color = prob, size = pval))+
  geom_point()+
  theme_bw()+
  scale_color_distiller(palette = "Spectral")+
  scale_size_manual(values = c("0.01 < p < 0.05" = 2, "p < 0.01" = 4)) +
  theme(axis.text.x = element_text(size = 10,color = 'black',vjust = 1,hjust = 1,angle = 35), 
        axis.text.y = element_text(size = 10,color = 'black'), 
        legend.text = element_text(size = 10,color = 'black'), 
        legend.title = element_text(size = 10,color = 'black'), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin=unit(c(1,1,1,1),'lines'),
        legend.position = 'right')+
  facet_grid(. ~ title)
dev.off()  




############ control NEGR1-_Ex-L2/4
temp<-netVisual_bubble(cellchat_Control, sources.use = c(5,6), targets.use = c(1:4,7:9), remove.isolate = F,return.data = T)$communication
temp<-na.omit(temp)
df.net_negnative<-temp[temp$source%in%c("NEGR1-_Ex-L2/4"),]
df.net_negnative<-df.net_negnative[,c("target","interaction_name_2","prob","pval")]

df.net_negnative1<-data.frame(target=rep(names(table(temp$target)),22),
                              interaction_name_2=rep(names(table(temp$interaction_name_2)),each = 7))
df.net_negnative2<-merge(df.net_negnative,df.net_negnative1,all=TRUE)
df.net_negnative2$prob[is.na(df.net_negnative2$prob)]<-0
df.net_negnative2$pval[is.na(df.net_negnative2$pval)]<-1

df.net_negnative2$pval[df.net_negnative2$pval=="1"]<-"p > 0.05"
df.net_negnative2$pval[df.net_negnative2$pval=="2"]<-"0.01 < p < 0.05"
df.net_negnative2$pval[df.net_negnative2$pval=="3"]<-"p < 0.01"
df.net_negnative2$pval<-factor(df.net_negnative2$pval,levels = c("p > 0.05", "0.01 < p < 0.05", "p < 0.01"))
df.net_negnative2$title<-"Control NEGR1- Ex-L2/4"

pdf("NC_Control_NEGR1-_Ex-L2_4_cellchat.pdf",height = 5,width = 5.5)
ggplot(df.net_negnative2,aes(x=target,y=interaction_name_2,color = prob,size = pval))+
  geom_point()+
  theme_bw()+
  scale_color_distiller(palette = "Spectral")+
  scale_size_manual(values = c("0.01 < p < 0.05" = 2, "p < 0.01" = 4)) +
  theme(axis.text.x = element_text(size = 10,color = 'black',vjust = 1,hjust = 1,angle = 35),
        axis.text.y = element_text(size = 10,color = 'black'),
        legend.text = element_text(size = 10,color = 'black'), 
        legend.title = element_text(size = 10,color = 'black'), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin=unit(c(1,1,1,1),'lines'),
        legend.position = 'right')+
  facet_grid(. ~ title)
dev.off()  


########### control NEGR1+_Ex-L2/4
df.net_positive<-temp[temp$source%in%c("NEGR1+_Ex-L2/4"),]
df.net_positive<-df.net_positive[,c("target","interaction_name_2","prob","pval")]
df.net_positive1<-data.frame(target=rep(names(table(temp$target)),22),
                             interaction_name_2=rep(names(table(temp$interaction_name_2)),each = 7))
df.net_positive2<-merge(df.net_positive,df.net_positive1,all=TRUE)
df.net_positive2$prob[is.na(df.net_positive2$prob)]<-0
df.net_positive2$pval[is.na(df.net_positive2$pval)]<-1

df.net_positive2$pval[df.net_positive2$pval=="1"]<-"p > 0.05"
df.net_positive2$pval[df.net_positive2$pval=="2"]<-"0.01 < p < 0.05"
df.net_positive2$pval[df.net_positive2$pval=="3"]<-"p < 0.01"
df.net_positive2$pval<-factor(df.net_positive2$pval,levels = c("p > 0.05", "0.01 < p < 0.05", "p < 0.01"))
df.net_positive2$title<-"Control NEGR1+ Ex-L2/4"

pdf("NC_Control_NEGR1+_Ex-L2_4_cellchat.pdf",height = 5,width = 5.5)
ggplot(df.net_positive2,aes(x=target,y=interaction_name_2,color = prob, size = pval))+
  geom_point()+
  theme_bw()+
  scale_color_distiller(palette = "Spectral")+
  scale_size_manual(values = c("0.01 < p < 0.05" = 2, "p < 0.01" = 4)) +
  theme(axis.text.x = element_text(size = 10,color = 'black',vjust = 1,hjust = 1,angle = 35), 
        axis.text.y = element_text(size = 10,color = 'black'), 
        legend.text = element_text(size = 10,color = 'black'), 
        legend.title = element_text(size = 10,color = 'black'), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin=unit(c(1,1,1,1),'lines'),
        legend.position = 'right')+
  facet_grid(. ~ title)
dev.off()  





