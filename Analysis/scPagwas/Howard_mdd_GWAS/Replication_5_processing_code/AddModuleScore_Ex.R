#################### Addmodule function validate the risk pathways in NC Ex-L2/4 ####################
library(Seurat)
library(scPagwas)
library(ggplot2)
library(rstatix)
library(ggpubr)
library(DOSE)
library(GOSemSim)
library(clusterProfiler)
library(org.Hs.eg.db)
setwd("D:/Project/SingleCell_MDD/Figure/NC_Ex-L2-4/AddmoduleScore/")

#################### import the singlecell data
single_data<-readRDS("D:/Project/SingleCell_MDD/SingleCell_analysis/Exitatory_neurons_subtype_analysis/NC_Excitatory_neurons_annotation/NC_NEGR1_EX_L2_4.rds")
Idents(single_data)<-single_data$NEGR1_type
single_data<-subset(single_data,idents=c("NEGR1+","NEGR1-"))
table(single_data$NEGR1_type)
single_data$NEGR1_type<-paste(single_data$NEGR1_type,single_data$phenotype,sep = "_")

#################### get the KEGG pathway genes
kegg_geneset<-Genes_by_pathway_kegg
#################### Glutamatergic synapse
Glutamatergic_synapse<-kegg_geneset["hsa04724"]
single_data <- AddModuleScore(object = single_data,features =Glutamatergic_synapse ,ctrl = 100, name = 'Glutamatergic_synapse_score')
data_score<-FetchData(single_data,vars = c("NEGR1_type","Glutamatergic_synapse_score1"))
data_score$NEGR1_type<-factor(data_score$NEGR1_type,levels = c("NEGR1+_Case","NEGR1-_Case","NEGR1+_Control","NEGR1-_Control"))
data_score1<-aggregate(data_score$Glutamatergic_synapse_score1,by=list(data_score$NEGR1_type),median)
groups <- list(c("NEGR1+_Case","NEGR1-_Case"),c("NEGR1+_Control","NEGR1-_Control"),
               c("NEGR1+_Case","NEGR1+_Control"),c("NEGR1-_Case","NEGR1-_Control"))

pdf("Glutamatergic_synapse_score.pdf",height = 4.5,width = 5.5)
ggplot(data_score, aes(x=NEGR1_type, y=Glutamatergic_synapse_score1)) + 
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
  ylab("Glutamatergic synapse score")+
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




geneset<-read.csv("D:/Project/SingleCell_MDD/Figure/Figure4/NEGR1+ vs NEGR1- GO BP.csv",row.names = 1)
G0_id<-geneset$ID
G0_id<-G0_id[G0_id!=""]

for (i in G0_id) {
  modulation_of_chemical_synaptic_transmission<- get(i, org.Hs.egGO2ALLEGS) %>% mget(org.Hs.egSYMBOL) %>% unlist()
  modulation_of_chemical_synaptic_transmission<-list(modulation_of_chemical_synaptic_transmission)
  single_data <- AddModuleScore(object = single_data,features =modulation_of_chemical_synaptic_transmission ,ctrl = 100, name = 'modulation_of_chemical_synaptic_transmission_score')
  data_score<-FetchData(single_data,vars = c("NEGR1_type","modulation_of_chemical_synaptic_transmission_score1"))
  data_score$NEGR1_type<-factor(data_score$NEGR1_type,levels = c("NEGR1+_Case","NEGR1-_Case","NEGR1+_Control","NEGR1-_Control"))
  data_score1<-aggregate(data_score$modulation_of_chemical_synaptic_transmission_score1,by=list(data_score$NEGR1_type),median)
  colnames(data_score1)<-c("NEGR1_type","median_score")
  if(data_score1$median_score[1]>data_score1$median_score[3]){
    i<-gsub(":","_",i)
    print(i)
    print(data_score1)
    write.csv(data_score1,file = paste0(i,"_median_score.csv"),quote = F)
  }
}



################# get the GO BP genes
#### vocalization behavior
vocalization_behavior<- get("GO:0071625", org.Hs.egGO2ALLEGS) %>% mget(org.Hs.egSYMBOL) %>% unlist()
vocalization_behavior<-list(vocalization_behavior)
single_data <- AddModuleScore(object = single_data,features =vocalization_behavior ,ctrl = 100, name = 'vocalization_behavior_score')
data_score<-FetchData(single_data,vars = c("NEGR1_type","vocalization_behavior_score1"))
data_score$NEGR1_type<-factor(data_score$NEGR1_type,levels = c("NEGR1+_Case","NEGR1-_Case","NEGR1+_Control","NEGR1-_Control"))
data_score1<-aggregate(data_score$vocalization_behavior_score1,by=list(data_score$NEGR1_type),median)
colnames(data_score1)<-c("NEGR1_type","median_score")
write.csv(data_score1,file = "vocalization_behavior.csv",quote = F)


pdf("GO_vocalization_behavior.pdf",height = 4.5,width = 5.5)
ggplot(data_score, aes(x=NEGR1_type, y=vocalization_behavior_score1)) + 
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
  ylab("Vocalization behavior score")+
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


#### regulation of basement membrane organization
regulation_of_basement_membrane_organization<- get("GO:0110011", org.Hs.egGO2ALLEGS) %>% mget(org.Hs.egSYMBOL) %>% unlist()
regulation_of_basement_membrane_organization<-list(regulation_of_basement_membrane_organization)
single_data <- AddModuleScore(object = single_data,features =regulation_of_basement_membrane_organization ,ctrl = 100, name = 'regulation_of_basement_membrane_organization')
data_score<-FetchData(single_data,vars = c("NEGR1_type","regulation_of_basement_membrane_organization1"))
data_score$NEGR1_type<-factor(data_score$NEGR1_type,levels = c("NEGR1+_Case","NEGR1-_Case","NEGR1+_Control","NEGR1-_Control"))
data_score1<-aggregate(data_score$regulation_of_basement_membrane_organization1,by=list(data_score$NEGR1_type),median)
colnames(data_score1)<-c("NEGR1_type","median_score")
write.csv(data_score1,file = "regulation_of_basement_membrane_organization.csv",quote = F)

pdf("GO_regulation_of_basement_membrane_organization.pdf",height = 4.5,width = 5.5)
ggplot(data_score, aes(x=NEGR1_type, y=regulation_of_basement_membrane_organization1)) + 
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
  ylab("Regulation of basement membrane organization score")+
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
