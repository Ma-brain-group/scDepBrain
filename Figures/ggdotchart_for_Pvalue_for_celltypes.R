#################### visualize result of the celltypes in four software ##################
rm(list = ls())
library(Seurat)
library(ggplot2)
library(ggvenn) 
library(cowplot)
library(ggpubr)

#################### set the file path
setwd("D:/Project/SingleCell_MDD/Figure/Figure2/")

#################### import the result
##magma
magma<-read.table("D:/Project/SingleCell_MDD/GWAS_combine_singlecell_analysis/1_GWAS_analysis_for_broad_celltype/magma/ieu_b_102.gsa.out",header = T)
##ldsc
ldsc<-read.table("D:/Project/SingleCell_MDD/GWAS_combine_singlecell_analysis/1_GWAS_analysis_for_broad_celltype/LDSC_SEG/ieu_b_102.cell_type_results.txt",header = T)
##scPagwas
scPagwas<-read.csv("D:/Project/SingleCell_MDD/GWAS_combine_singlecell_analysis/1_GWAS_analysis_for_broad_celltype/scPagwas/scPagwas_singlecell_result_block_annotation_hg37_1/ieu_b_102/Clinically-derived MDD_celltype_result.csv",header = T,sep = ",")
scPagwas[scPagwas$celltype=="Excitatory neurons",1]<-"Excitatory.neurons"
scPagwas[scPagwas$celltype=="Inhibitory neurons",1]<-"Inhibitory.neurons"
scPagwas[scPagwas$celltype=="Purkinje neurons",1]<-"Purkinje.neurons"
scPagwas[scPagwas$celltype=="Endothelial cells",1]<-"Endothelial.cells"
##scDRS
scDRS<-read.table("D:/Project/SingleCell_MDD/GWAS_combine_singlecell_analysis/1_GWAS_analysis_for_broad_celltype/scDRS/MDD.scdrs_ct.anno",header = T,sep = "\t")
scDRS<-scDRS[,c("X","assoc_mcp")]

#################### get the celltype and P value
magma<-magma[,c(1,7)]
ldsc<-ldsc[,c(1,4)]
colnames(magma)<-c("cell_type","P")
colnames(ldsc)<-c("cell_type","P")
colnames(scDRS)<-c("cell_type","P")
colnames(scPagwas)<-c("cell_type","P")

#################### calculate the number of Significant celltypes
magma[magma$P<0.05,3]<-1
magma[magma$P>0.05,3]<-0
ldsc[ldsc$P<0.05,3]<-1
ldsc[ldsc$P>0.05,3]<-0
scDRS[scDRS$P<0.05,3]<-1
scDRS[scDRS$P>0.05,3]<-0
scPagwas[scPagwas$P<0.05,3]<-1
scPagwas[scPagwas$P>0.05,3]<-0

#################### change the character to factor
colnames(magma)[3]<-"Significance"
colnames(ldsc)[3]<-"Significance"
colnames(scDRS)[3]<-"Significance"
colnames(scPagwas)[3]<-"Significance"
magma$Significance<-factor(magma$Significance)
ldsc$Significance<-factor(ldsc$Significance)
scDRS$Significance<-factor(scDRS$Significance)
scPagwas$Significance<-factor(scPagwas$Significance)

#################### -log10(P value)
magma$`-log10(P)`<--log10(magma$P)
ldsc$`-log10(P)`<--log10(ldsc$P)
scDRS$`-log10(P)`<--log10(scDRS$P)
scPagwas$`-log10(P)`<--log10(scPagwas$P)


#################### magma plot
magma$cell_type<-factor(magma$cell_type,levels=rev(c("Excitatory.neurons","Inhibitory.neurons","Purkinje.neurons","Endothelial.cells",
                                                     "Oligodendrocytes","OPCs","Microglia","Astrocytes")))

pdf("ggdotchart_of_magma.pdf",width = 4.5,height = 3)
ggdotchart(magma, x = "cell_type", y = "-log10(P)",
           color = "cell_type",                               ## 指定上色组别
           palette = c("#d5231d","#e88f18","#b698c5","#e47faf","#3777ac","#a05528","#8e4c99","#4ea64a"),
           sorting = "descending",                       # Sort value in descending order
           add = "segments",                             # Add segments from y = 0 to dots
           add.params = list(color = "lightgray", size = 2), # Change segment color and size
           group = "cell_type",                                # Order by groups
           dot.size = 7,                                 # Large dot size
           label = round(magma$`-log10(P)`,3), # Add mpg values as dot labels
           font.label = list(color = "black", size = 8, vjust = 0.5),               # Adjust label parameters
           ggtheme = theme_pubr())+
  coord_flip()+
  geom_hline(yintercept = -log10(0.05), linetype = 2, color = "lightgray")+
  theme(axis.text.x = element_text (size = 10,color = 'black',angle = 0,hjust = 0.5),
        axis.text.y = element_text(size = 10,color = 'black'),
        legend.text = element_text(size = 10,color = 'black'),
        legend.title = element_text(size = 10,color = 'black'),
        legend.position = 'none'
  )+
  scale_y_continuous(limits = c (0, 10), breaks = seq (0, 10, 2))
dev.off()



#################### ldsc plot
ldsc$cell_type<-factor(ldsc$cell_type,levels=rev(c("Excitatory.neurons","Inhibitory.neurons","Purkinje.neurons","Endothelial.cells",
                                                     "Oligodendrocytes","OPCs","Microglia","Astrocytes")))

pdf("ggdotchart_of_ldsc.pdf",width = 4.5,height = 3)
ggdotchart(ldsc, x = "cell_type", y = "-log10(P)",
           color = "cell_type",                               ## 指定上色组别
           palette = c("#d5231d","#e88f18","#b698c5","#e47faf","#3777ac","#a05528","#8e4c99","#4ea64a"),
           sorting = "descending",                       # Sort value in descending order
           add = "segments",                             # Add segments from y = 0 to dots
           add.params = list(color = "lightgray", size = 2), # Change segment color and size
           group = "cell_type",                                # Order by groups
           dot.size = 7,                                 # Large dot size
           label = round(ldsc$`-log10(P)`,3), # Add mpg values as dot labels
           font.label = list(color = "black", size = 8, vjust = 0.5),               # Adjust label parameters
           ggtheme = theme_pubr())+
  coord_flip()+
  geom_hline(yintercept = -log10(0.05), linetype = 2, color = "lightgray")+
  theme(axis.text.x = element_text (size = 10,color = 'black',angle = 0,hjust = 0.5),
        axis.text.y = element_text(size = 10,color = 'black'),
        legend.text = element_text(size = 10,color = 'black'),
        legend.title = element_text(size = 10,color = 'black'),
        legend.position = 'none'
  )+
  scale_y_continuous(limits = c (0, 10), breaks = seq (0, 10, 2))
dev.off()



#################### scDRS plot
scDRS$cell_type<-factor(scDRS$cell_type,levels=rev(c("Excitatory neurons","Inhibitory neurons","Purkinje neurons","Endothelial cells",
                                                   "Oligodendrocytes","OPCs","Microglia","Astrocytes")))

pdf("ggdotchart_of_scDRS.pdf",width = 4.5,height = 3)
ggdotchart(scDRS, x = "cell_type", y = "-log10(P)",
           color = "cell_type",                               ## 指定上色组别
           palette = c("#d5231d","#e88f18","#b698c5","#e47faf","#3777ac","#a05528","#8e4c99","#4ea64a"),
           sorting = "descending",                       # Sort value in descending order
           add = "segments",                             # Add segments from y = 0 to dots
           add.params = list(color = "lightgray", size = 2), # Change segment color and size
           group = "cell_type",                                # Order by groups
           dot.size = 7,                                 # Large dot size
           label = round(scDRS$`-log10(P)`,3), # Add mpg values as dot labels
           font.label = list(color = "black", size = 8, vjust = 0.5),               # Adjust label parameters
           ggtheme = theme_pubr())+
  coord_flip()+
  geom_hline(yintercept = -log10(0.05), linetype = 2, color = "lightgray")+
  theme(axis.text.x = element_text (size = 10,color = 'black',angle = 0,hjust = 0.5),
        axis.text.y = element_text(size = 10,color = 'black'),
        legend.text = element_text(size = 10,color = 'black'),
        legend.title = element_text(size = 10,color = 'black'),
        legend.position = 'none'
  )+
  scale_y_continuous(limits = c (0, 10), breaks = seq (0, 10, 2))
dev.off()




#################### scPagwas plot
scPagwas$cell_type<-factor(scPagwas$cell_type,levels=rev(c("Excitatory.neurons","Inhibitory.neurons","Purkinje.neurons","Endothelial.cells",
                                                     "Oligodendrocytes","OPCs","Microglia","Astrocytes")))

pdf("ggdotchart_of_scPagwas.pdf",width = 4.5,height = 3)
ggdotchart(scPagwas, x = "cell_type", y = "-log10(P)",
           color = "cell_type",                               ## 指定上色组别
           palette = c("#d5231d","#e88f18","#b698c5","#e47faf","#3777ac","#a05528","#8e4c99","#4ea64a"),
           sorting = "descending",                       # Sort value in descending order
           add = "segments",                             # Add segments from y = 0 to dots
           add.params = list(color = "lightgray", size = 2), # Change segment color and size
           group = "cell_type",                                # Order by groups
           dot.size = 7,                                 # Large dot size
           label = round(scPagwas$`-log10(P)`,3), # Add mpg values as dot labels
           font.label = list(color = "black", size = 8, vjust = 0.5),               # Adjust label parameters
           ggtheme = theme_pubr())+
  coord_flip()+
  geom_hline(yintercept = -log10(0.05), linetype = 2, color = "lightgray")+
  theme(axis.text.x = element_text (size = 10,color = 'black',angle = 0,hjust = 0.5),
        axis.text.y = element_text(size = 10,color = 'black'),
        legend.text = element_text(size = 10,color = 'black'),
        legend.title = element_text(size = 10,color = 'black'),
        legend.position = 'none'
  )+
  scale_y_continuous(limits = c (0, 10), breaks = seq (0, 10, 2))
dev.off()
