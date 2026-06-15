################## In_PVALB CNNM2+ vs CNNM2- enrichment analysis ########################
library(Seurat)
library(DOSE)
library(org.Hs.eg.db)
library(topGO)
library(clusterProfiler)
library(pathview)
library(enrichplot)
library(ggplot2)
setwd("D:/Project/SingleCell_MDD/Figure/Figure3/")

################## import the In_PVALB singlecell data
single_data<-readRDS("D:/Project/SingleCell_MDD/SingleCell_analysis/11_CNNM2_positive_vs_CNNM2_negnative_IN6/CNNM2_IN6.rds")
table(single_data$CNNM2_type)

################## calculate the DEGs
markers<-FindMarkers(single_data,assay = "RNA",group.by = "CNNM2_type",ident.1 = "CNNM2+",ident.2 = "CNNM2-",only.pos = T)
markers<-markers[markers$p_val_adj<0.05,]
write.csv(markers,file = "CNNM2_positive_vs_CNNM2_negnative_In_PVALB_DEG.csv",quote = F)

################## GO BP analysis
markers <- bitr(row.names(markers), fromType="SYMBOL",toType="ENTREZID", OrgDb="org.Hs.eg.db")
ego_up <- enrichGO(markers$ENTREZID,
                   OrgDb = "org.Hs.eg.db",
                   keyType="ENTREZID",
                   pAdjustMethod="BH",
                   ont="BP",
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05,
                   readable = T)
### simplify the GO terms
ego_up <- clusterProfiler::simplify(ego_up, cutoff=0.7, by="p.adjust", select_fun=min)
temp<-ego_up@result
temp<-temp[temp$p.adjust<0.05,]
write.csv(temp,file = "CNNM2+ vs CNNM2- GO BP.csv",quote = F)

### plot top 15 terms
temp$`-log10(fdr)`<- -log10(temp$p.adjust)
p<-ggplot(data = temp[1:10,], 
          aes(`-log10(fdr)`, reorder(Description,`-log10(fdr)`))) +
  geom_bar(stat="identity",
           alpha=0.5,
           fill="#E1884A",
           color="black",
           width = 0.8) + 
  theme_classic(base_size = 12)+
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
  )+
  geom_vline(aes(xintercept=-log10(0.05)),color="grey",linetype="dashed")
ggsave("CNNM2+ vs CNNM2- GO BP.pdf",p,height = 4,width = 7)


################# KEGG analysis
kk <- enrichKEGG(markers$ENTREZID,
                 keyType = "kegg",
                 organism="hsa",
                 pAdjustMethod = "BH",
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 0.05)
kk<-setReadable(kk, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
temp1<-kk@result
temp1<-temp1[temp1$p.adjust<0.05,]
write.csv(temp1,file = "CNNM2+ vs CNNM2- KEGG.csv",quote = F)

## plot top 10 terms
temp1$`-log10(fdr)`<- -log10(temp1$p.adjust)
p1<-ggplot(data = temp1[1:10,], 
           aes(`-log10(fdr)`, reorder(Description,`-log10(fdr)`))) +
  geom_bar(stat="identity",
           alpha=0.5,
           fill="#E1884A",
           color="black",
           width = 0.8) + 
  theme_classic(base_size = 12)+
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
  )+
  geom_vline(aes(xintercept=-log10(0.05)),color="grey",linetype="dashed")
ggsave("CNNM2+ vs CNNM2- KEGG.pdf",p1,height = 4,width = 7)











