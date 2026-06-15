#################### Addmodule function validate the risk pathways in NC In_PVALB ####################
library(Seurat)
library(scPagwas)
library(ggplot2)
library(rstatix)
library(ggpubr)
library(DOSE)
library(GOSemSim)
library(clusterProfiler)
library(org.Hs.eg.db)
setwd("D:/Project/SingleCell_MDD/Figure/NC_In_PVALB/AddmoduleScore/")

#################### import the singlecell data
single_data<-readRDS("D:/Project/SingleCell_MDD/SingleCell_analysis/16_validation_NC/2_INs_subset/NC_In_PVALB.rds")
table(single_data$CNNM2_type)


############### get the KEGG pathway genes
kegg_geneset<-Genes_by_pathway_kegg
############## Glutamatergic synapse
Glutamatergic_synapse<-kegg_geneset["hsa04724"]
single_data <- AddModuleScore(object = single_data,features =Glutamatergic_synapse ,ctrl = 100, name = 'Glutamatergic_synapse_score')
data_score<-FetchData(single_data,vars = c("CNNM2_type","Glutamatergic_synapse_score1"))
data_score$CNNM2_type<-factor(data_score$CNNM2_type,levels = c("CNNM2+_Case","CNNM2-_Case","CNNM2+_Control","CNNM2-_Control"))
data_score1<-aggregate(data_score$Glutamatergic_synapse_score1,by=list(data_score$CNNM2_type),median)
groups <- list(c("CNNM2+_Case","CNNM2-_Case"),c("CNNM2+_Control","CNNM2-_Control"),
               c("CNNM2+_Case","CNNM2+_Control"),c("CNNM2-_Case","CNNM2-_Control"))

pdf("Glutamatergic_synapse_score.pdf",height = 4.5,width = 5.5)
ggplot(data_score, aes(x=CNNM2_type, y=Glutamatergic_synapse_score1)) + 
  stat_boxplot(geom = "errorbar",width=0.05, size=0.5,position=position_dodge(0.6),color= "black")+
  theme_bw()+
  geom_boxplot(position = position_dodge(0.6),
               size = 0.5,
               width = 0.8,
               fill = c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0'),
               color = "black",
               outlier.color = "black",
               outlier.fill = "black",
               outlier.shape = 19,
               outlier.size = 1,
               outlier.stroke = 0.5,
               outlier.alpha = 45,
               notch = F,
               notchwidth = 0.5)+
  xlab("")+
  ylab("Glutamatergic synapse score")+
  theme(axis.text.x = element_text(size = 10,color = 'black',vjust = 1,hjust = 0.5,angle = 0),
        axis.text.y = element_text(size = 10,color = 'black'), 
        legend.text = element_text(size = 10,color = 'black'), 
        legend.title = element_text(size = 10,color = 'black'), 
        panel.grid.major =  element_blank(),
        panel.grid.minor = element_blank())+
  stat_compare_means(comparisons=groups,
                     method = "wilcox.test",
                     label = "p.format",
                     size=4.5)
dev.off()





################# get the GO BP genes
#### modulation of chemical synaptic transmission
modulation_of_chemical_synaptic_transmission<- get("GO:0050804", org.Hs.egGO2ALLEGS) %>% mget(org.Hs.egSYMBOL) %>% unlist()
modulation_of_chemical_synaptic_transmission<-list(modulation_of_chemical_synaptic_transmission)
single_data <- AddModuleScore(object = single_data,features =modulation_of_chemical_synaptic_transmission ,ctrl = 500, name = 'modulation_of_chemical_synaptic_transmission_score')
data_score<-FetchData(single_data,vars = c("CNNM2_type","modulation_of_chemical_synaptic_transmission_score1"))
data_score$CNNM2_type<-factor(data_score$CNNM2_type,levels = c("CNNM2+_Case","CNNM2-_Case","CNNM2+_Control","CNNM2-_Control"))
data_score1<-aggregate(data_score$modulation_of_chemical_synaptic_transmission_score1,by=list(data_score$CNNM2_type),median)
colnames(data_score1)<-c("CNNM2_type","median_score")
write.csv(data_score1,file = "modulation_of_chemical_synaptic_transmission.csv",quote = F)


pdf("GO_modulation_of_chemical_synaptic_transmission.pdf",height = 4.5,width = 5.5)
ggplot(data_score, aes(x=CNNM2_type, y=modulation_of_chemical_synaptic_transmission_score1)) + 
  stat_boxplot(geom = "errorbar",width=0.05, size=0.5,position=position_dodge(0.6),color= "black")+
  theme_bw()+
  geom_boxplot(position = position_dodge(0.6),
               size = 0.5,
               width = 0.65,
               fill = c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0'),
               color = "black",
               outlier.color = "black",
               outlier.fill = "black",
               outlier.shape = 19,
               outlier.size = 1,
               outlier.stroke = 0.5,
               outlier.alpha = 45,
               notch = F,
               notchwidth = 0.5)+
  xlab("")+
  ylab("Modulation of chemical synaptic transmission score")+
  theme(axis.text.x = element_text(size = 10,color = 'black',vjust = 1,hjust = 0.5,angle = 0),
        axis.text.y = element_text(size = 10,color = 'black'), 
        legend.text = element_text(size = 10,color = 'black'), 
        legend.title = element_text(size = 10,color = 'black'), 
        panel.grid.major =  element_blank(),
        panel.grid.minor = element_blank())+
  stat_compare_means(comparisons=groups,
                     method = "wilcox.test",
                     label = "p.signif",
                     size=4.5)
dev.off()


#### regulation of nervous system process
regulation_of_nervous_system_process<- get("GO:0031644", org.Hs.egGO2ALLEGS) %>% mget(org.Hs.egSYMBOL) %>% unlist()
regulation_of_nervous_system_process<-list(regulation_of_nervous_system_process)
single_data <- AddModuleScore(object = single_data,features =regulation_of_nervous_system_process ,ctrl = 100, name = 'regulation_of_nervous_system_process')
data_score<-FetchData(single_data,vars = c("CNNM2_type","regulation_of_nervous_system_process1"))
data_score$CNNM2_type<-factor(data_score$CNNM2_type,levels = c("CNNM2+_Case","CNNM2-_Case","CNNM2+_Control","CNNM2-_Control"))
data_score1<-aggregate(data_score$regulation_of_nervous_system_process1,by=list(data_score$CNNM2_type),median)
colnames(data_score1)<-c("CNNM2_type","median_score")
write.csv(data_score1,file = "regulation_of_nervous_system_process.csv",quote = F)

pdf("GO_regulation_of_nervous_system_process.pdf",height = 4.5,width = 5.5)
ggplot(data_score, aes(x=CNNM2_type, y=regulation_of_nervous_system_process1)) + 
  stat_boxplot(geom = "errorbar",width=0.05, size=0.5,position=position_dodge(0.6),color= "black")+
  theme_bw()+
  geom_boxplot(position = position_dodge(0.6),
               size = 0.5,
               width = 0.65,
               fill = c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0'),
               color = "black",
               outlier.color = "black",
               outlier.fill = "black",
               outlier.shape = 19,
               outlier.size = 1,
               outlier.stroke = 0.5,
               outlier.alpha = 45,
               notch = F,
               notchwidth = 0.5)+
  xlab("")+
  ylab("regulation of nervous system process")+
  theme(axis.text.x = element_text(size = 10,color = 'black',vjust = 1,hjust = 0.5,angle = 0),
        axis.text.y = element_text(size = 10,color = 'black'), 
        legend.text = element_text(size = 10,color = 'black'), 
        legend.title = element_text(size = 10,color = 'black'), 
        panel.grid.major =  element_blank(),
        panel.grid.minor = element_blank())+
  stat_compare_means(comparisons=groups,
                     method = "wilcox.test",
                     label = "p.signif",
                     size=4.5)
dev.off()



################# cell-cell adhesion via plasma-membrane adhesion molecules
cell_cell_adhesion_via_plasma_membrane_adhesion_molecules<- get("GO:0098742", org.Hs.egGO2ALLEGS) %>% mget(org.Hs.egSYMBOL) %>% unlist()
cell_cell_adhesion_via_plasma_membrane_adhesion_molecules<-list(cell_cell_adhesion_via_plasma_membrane_adhesion_molecules)
single_data <- AddModuleScore(object = single_data,features =cell_cell_adhesion_via_plasma_membrane_adhesion_molecules ,ctrl = 500, name = 'cell_cell_adhesion_via_plasma_membrane_adhesion_molecules_score')
data_score<-FetchData(single_data,vars = c("CNNM2_type","cell_cell_adhesion_via_plasma_membrane_adhesion_molecules_score1"))
data_score$CNNM2_type<-factor(data_score$CNNM2_type,levels = c("CNNM2+_Case","CNNM2-_Case","CNNM2+_Control","CNNM2-_Control"))
data_score1<-aggregate(data_score$cell_cell_adhesion_via_plasma_membrane_adhesion_molecules_score1,by=list(data_score$CNNM2_type),median)
colnames(data_score1)<-c("CNNM2_type","median_score")
write.csv(data_score1,file = "cell_cell_adhesion_via_plasma_membrane_adhesion_molecules.csv",quote = F)

pdf("GO_cell_cell_adhesion_via_plasma_membrane_adhesion_molecules_score.pdf",height = 4.5,width = 5.5)
ggplot(data_score, aes(x=CNNM2_type, y=cell_cell_adhesion_via_plasma_membrane_adhesion_molecules_score1)) + 
  stat_boxplot(geom = "errorbar",width=0.05, size=0.5,position=position_dodge(0.6),color= "black")+
  theme_bw()+
  geom_boxplot(position = position_dodge(0.6),
               size = 0.5,
               width = 0.65,
               fill = c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0'),
               color = "black",
               outlier.color = "black",
               outlier.fill = "black",
               outlier.shape = 19,
               outlier.size = 1,
               outlier.stroke = 0.5,
               outlier.alpha = 45,
               notch = F,
               notchwidth = 0.5)+
  xlab("")+
  ylab("cell-cell adhesion via plasma-membrane adhesion molecules")+
  theme(axis.text.x = element_text(size = 10,color = 'black',vjust = 1,hjust = 0.5,angle = 0),
        axis.text.y = element_text(size = 10,color = 'black'), 
        legend.text = element_text(size = 10,color = 'black'), 
        legend.title = element_text(size = 10,color = 'black'), 
        panel.grid.major =  element_blank(),
        panel.grid.minor = element_blank())+
  stat_compare_means(comparisons=groups,
                     method = "wilcox.test",
                     label = "p.signif",
                     size=4.5)
dev.off()








