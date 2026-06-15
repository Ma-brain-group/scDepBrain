################ Stratified analysis visualize the result of Magma #####################
rm(list = ls())
library(Seurat)
library(ggplot2)
library(ggvenn) 
library(cowplot)
library(ggpubr)

############ change the file path
setwd("D:/Project/SingleCell_MDD/GWAS_combine_singlecell_analysis/6_分层分析/Sex/Magma/") 
list.files(getwd())

############ plot the result
file=list.files(getwd())
for (i in 1:length(file)) {
  print(file[i])
  magma<-read.table(file[i],header = T)
  #### get the information 
  magma<-magma[,c(1,7)]
  colnames(magma)<-c("cell_type","P")
  #### -log10(P)
  magma$`-log10(P)`<--log10(magma$P)
  #### change the character to factor
  magma$cell_type<-factor(magma$cell_type,levels=rev(c("Excitatory.neurons","Inhibitory.neurons","Purkinje.neurons","Endothelial.cells",
                                                     "Oligodendrocytes","OPCs","Microglia","Astrocytes")))
  #### plot
  p<-ggdotchart(magma, x = "cell_type", y = "-log10(P)",
                   color = "cell_type",                               
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
  ggsave(paste0(file[i],"_magma.pdf"),p,width = 4.5,height = 3)
  
  magma<-NULL
}
