################# TRS plot of brain region ################
setwd("D:/Project/SingleCell_MDD/SingleCell_analysis/4_TRS_brain_region_plot/")
library(Seurat)
library(ggplot2)
library(cowplot)

################ import the singlecell data
single_data<-readRDS("D:/Project/SingleCell_MDD/SingleCell_analysis/1_Singlecell_data/simple.rds")


################ boxplot for all cells
metadata<-single_data@meta.data
metadata_brain_region<-metadata[,c("tissue.region","scPagwas.TRS.Score")]
metadata_brain_region1<-aggregate(metadata_brain_region$scPagwas.TRS.Score,by=list(metadata_brain_region$tissue.region),FUN=median)
metadata_brain_region1$color<-c("#ED4437","#E1884A","#8ACC72","#1F78B4","#89C8E8","#B3446C","#EBD57C","#E68FAC","#CAA2F4",
                                "#96873B","#B49D99","#B37557","#FC9A9A","#6A3D9A")
metadata_brain_region1<-metadata_brain_region1[order(metadata_brain_region1$x,decreasing = T),]

######top的脑区为M1C，MTG，S1C，绘图
metadata_brain_region$tissue.region<-factor(metadata_brain_region$tissue.region,levels=metadata_brain_region1$Group.1)

pdf("TRS_brain_region_plot_in_all_cells.pdf",height = 5,width = 6)
ggplot(metadata_brain_region, aes(x=tissue.region, y=scPagwas.TRS.Score)) + 
  stat_boxplot(geom = "errorbar",width=0.05, size=0.5,position=position_dodge(0.6),color= "black")+
  theme_bw()+
  geom_boxplot(position = position_dodge(0.6),
               size = 0.5,
               width = 0.85,
               fill = metadata_brain_region1$color,
               color = "black",
               outlier.color = "black",
               outlier.fill = "black",
               outlier.shape = 19,
               outlier.size = 1,
               outlier.stroke = 0.5,
               outlier.alpha = 45,
               notch = F,
               notchwidth = 0.5)+
  xlab("Brain region")+
  ylab("scPagwas TRS Score")+
  ggtitle("All cells")+
  theme(axis.text.x = element_text(size = 10,color = 'black',vjust = 1,hjust = 1,angle = 35),
        axis.text.y = element_text(size = 10,color = 'black'),
        legend.text = element_text(size = 10,color = 'black'), 
        legend.title = element_text(size = 10,color = 'black'), 
        plot.title = element_text(size=12,hjust=0.5), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
dev.off()


############# subset the Inhibitory neurons
Idents(single_data)<-single_data$anno
single_data<-subset(single_data,idents = "Inhibitory neurons")

metadata<-single_data@meta.data
metadata_brain_region<-metadata[,c("tissue.region","scPagwas.TRS.Score")]

metadata_brain_region1<-aggregate(metadata_brain_region$scPagwas.TRS.Score,by=list(metadata_brain_region$tissue.region),FUN=median)
metadata_brain_region1$color<-c("#ED4437","#E1884A","#8ACC72","#1F78B4","#89C8E8","#B3446C","#EBD57C","#E68FAC","#CAA2F4",
                                "#96873B","#B49D99","#B37557","#FC9A9A","#6A3D9A")
metadata_brain_region1<-metadata_brain_region1[order(metadata_brain_region1$x,decreasing = T),]
metadata_brain_region$tissue.region<-factor(metadata_brain_region$tissue.region,levels=metadata_brain_region1$Group.1)

############# boxplot for Inhibitory neurons
pdf("TRS_brain_region_plot_in_Inhibitory_neurons.pdf",height = 5,width = 6)
ggplot(metadata_brain_region, aes(x=tissue.region, y=scPagwas.TRS.Score)) + 
  stat_boxplot(geom = "errorbar",width=0.05, size=0.5,position=position_dodge(0.6),color= "black")+
  theme_bw()+
  geom_boxplot(position = position_dodge(0.6),
               size = 0.5,
               width = 0.85,
               fill = metadata_brain_region1$color,
               color = "black",
               outlier.color = "black",
               outlier.fill = "black",
               outlier.shape = 19,
               outlier.size = 1,
               outlier.stroke = 0.5,
               outlier.alpha = 45,
               notch = F,
               notchwidth = 0.5)+
  xlab("Brain region")+
  ylab("scPagwas TRS Score")+
  ggtitle("Inhibitory neurons")+
  theme(axis.text.x = element_text(size = 10,color = 'black',vjust = 1,hjust = 1,angle = 35), 
        axis.text.y = element_text(size = 10,color = 'black'), 
        legend.text = element_text(size = 10,color = 'black'), 
        legend.title = element_text(size = 10,color = 'black'), 
        plot.title = element_text(size=12,hjust=0.5), 
        panel.grid.major =  element_blank(),
        panel.grid.minor = element_blank())
dev.off()


############## boxplot for In_PVALB
single_data<-readRDS("D:/Project/SingleCell_MDD/SingleCell_analysis/1_Singlecell_data/IN.rds")

Idents(single_data)<-single_data$final_anno
single_data<-subset(single_data,idents = "In6")

metadata<-single_data@meta.data
metadata_brain_region<-metadata[,c("tissue.region","scPagwas.TRS.Score")]

metadata_brain_region1<-aggregate(metadata_brain_region$scPagwas.TRS.Score,by=list(metadata_brain_region$tissue.region),FUN=median)
metadata_brain_region1$color<-c("#ED4437","#E1884A","#8ACC72","#1F78B4","#89C8E8","#B3446C","#EBD57C","#E68FAC","#CAA2F4",
                                "#96873B","#B49D99","#B37557","#FC9A9A","#6A3D9A")
metadata_brain_region1<-metadata_brain_region1[order(metadata_brain_region1$x,decreasing = T),]
metadata_brain_region$tissue.region<-factor(metadata_brain_region$tissue.region,levels=metadata_brain_region1$Group.1)

pdf("TRS_brain_region_plot_in_In_PVALB.pdf",height = 5,width = 6)
ggplot(metadata_brain_region, aes(x=tissue.region, y=scPagwas.TRS.Score)) + 
  stat_boxplot(geom = "errorbar",width=0.05, size=0.5,position=position_dodge(0.6),color= "black")+
  theme_bw()+
  geom_boxplot(position = position_dodge(0.6),
               size = 0.5,
               width = 0.85,
               fill = metadata_brain_region1$color,
               color = "black",
               outlier.color = "black",
               outlier.fill = "black",
               outlier.shape = 19,
               outlier.size = 1,
               outlier.stroke = 0.5,
               outlier.alpha = 45,
               notch = F,
               notchwidth = 0.5)+
  xlab("Brain region")+
  ylab("scPagwas TRS Score")+
  ggtitle("In_PVALB")+
  theme(axis.text.x = element_text(size = 10,color = 'black',vjust = 1,hjust=1,angle = 35), 
        axis.text.y = element_text(size = 10,color = 'black'), 
        legend.text = element_text(size = 10,color = 'black'), 
        legend.title = element_text(size = 10,color = 'black'), 
        plot.title = element_text(size=12,hjust=0.5), 
        panel.grid.major =  element_blank(),
        panel.grid.minor = element_blank())
dev.off()




################### import the Excitatory neurons singlecell data
single_data<-readRDS("D:/Project/SingleCell_MDD/SingleCell_analysis/annotation/Excitatory_neuron_subtypes_reannotation.rds")
metadata<-single_data@meta.data
metadata_brain_region<-metadata[,c("tissue.region","scPagwas.TRS.Score")]

metadata_brain_region1<-aggregate(metadata_brain_region$scPagwas.TRS.Score,by=list(metadata_brain_region$tissue.region),FUN=median)
metadata_brain_region1$color<-c("#ED4437","#E1884A","#8ACC72","#1F78B4","#89C8E8","#B3446C","#EBD57C","#E68FAC","#CAA2F4",
                                "#96873B","#B49D99","#B37557","#FC9A9A","#6A3D9A")
metadata_brain_region1<-metadata_brain_region1[order(metadata_brain_region1$x,decreasing = T),]
metadata_brain_region$tissue.region<-factor(metadata_brain_region$tissue.region,levels=metadata_brain_region1$Group.1)

pdf("TRS_brain_region_plot_in_Excitatory_neurons.pdf",height = 5,width = 6)
ggplot(metadata_brain_region, aes(x=tissue.region, y=scPagwas.TRS.Score)) + 
  stat_boxplot(geom = "errorbar",width=0.05, size=0.5,position=position_dodge(0.6),color= "black")+
  theme_bw()+
  geom_boxplot(position = position_dodge(0.6),
               size = 0.5,
               width = 0.85,
               fill = metadata_brain_region1$color,
               color = "black",
               outlier.color = "black",
               outlier.fill = "black",
               outlier.shape = 19,
               outlier.size = 1,
               outlier.stroke = 0.5,
               outlier.alpha = 45,
               notch = F,
               notchwidth = 0.5)+
  xlab("Brain region")+
  ylab("scPagwas TRS Score")+
  ggtitle("Excitatory neurons")+
  theme(axis.text.x = element_text(size = 10,color = 'black',vjust = 1,hjust=1,angle = 35), 
        axis.text.y = element_text(size = 10,color = 'black'), 
        legend.text = element_text(size = 10,color = 'black'), 
        legend.title = element_text(size = 10,color = 'black'), 
        plot.title = element_text(size=12,hjust=0.5), 
        panel.grid.major =  element_blank(),
        panel.grid.minor = element_blank())
dev.off()



################### import the Ex-L2/4 singlecell data
single_data_Ex_L2_4<-readRDS("D:/Project/SingleCell_MDD/SingleCell_analysis/Exitatory_neurons_subtype_analysis/Ex_L2_4_NEGR_analysis/NEGR1_EX_L2_4.rds")
metadata<-single_data_Ex_L2_4@meta.data
metadata_brain_region<-metadata[,c("tissue.region","scPagwas.TRS.Score")]

metadata_brain_region1<-aggregate(metadata_brain_region$scPagwas.TRS.Score,by=list(metadata_brain_region$tissue.region),FUN=median)
metadata_brain_region1$color<-c("#ED4437","#E1884A","#8ACC72","#1F78B4","#89C8E8","#B3446C","#EBD57C","#E68FAC","#CAA2F4",
                                "#96873B","#B49D99","#B37557","#FC9A9A","#6A3D9A")
metadata_brain_region1<-metadata_brain_region1[order(metadata_brain_region1$x,decreasing = T),]
metadata_brain_region$tissue.region<-factor(metadata_brain_region$tissue.region,levels=metadata_brain_region1$Group.1)

pdf("TRS_brain_region_plot_in_Ex_L2_4.pdf",height = 5,width = 6)
ggplot(metadata_brain_region, aes(x=tissue.region, y=scPagwas.TRS.Score)) + 
  stat_boxplot(geom = "errorbar",width=0.05, size=0.5,position=position_dodge(0.6),color= "black")+
  theme_bw()+
  geom_boxplot(position = position_dodge(0.6),
               size = 0.5,
               width = 0.85,
               fill = metadata_brain_region1$color,
               color = "black",
               outlier.color = "black",
               outlier.fill = "black",
               outlier.shape = 19,
               outlier.size = 1,
               outlier.stroke = 0.5,
               outlier.alpha = 45,
               notch = F,
               notchwidth = 0.5)+
  xlab("Brain region")+
  ylab("scPagwas TRS Score")+
  ggtitle("Ex-L2/4")+
  theme(axis.text.x = element_text(size = 10,color = 'black',vjust = 1,hjust=1,angle = 35), 
        axis.text.y = element_text(size = 10,color = 'black'), 
        legend.text = element_text(size = 10,color = 'black'), 
        legend.title = element_text(size = 10,color = 'black'), 
        plot.title = element_text(size=12,hjust=0.5), 
        panel.grid.major =  element_blank(),
        panel.grid.minor = element_blank())
dev.off()
