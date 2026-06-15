##################### Ex-L2/4 regulon plot ####################
rm(list = ls())
library(ggplot2)
library(dplyr)
library(tidyr)
setwd("D:/Project/SingleCell_MDD/Figure/Figure5/pyscenic_Ex-L2_4/")

##################### import the regulon result
regulon<-read.csv("SingleCell.out.regulons.csv",header = T)
regulon<-regulon[c(-1,-2),]

##################### get the each TF target genes
extract_genes<-function(info) {
  genes<-regmatches(info, gregexpr("'\\w+'", info))[[1]]
  genes<-gsub("'", "", genes)
  return(genes)
}
regulon$target_gene<-lapply(regulon$Enrichment.6, extract_genes)
regulon$target_gene<- sapply(regulon$target_gene, paste, collapse = ", ")

##################### for each TF, get the target genes and remove duplicates
regulon1 <- regulon %>%
  mutate(target_gene = strsplit(as.character(target_gene), ",")) %>%
  unnest(target_gene)
regulon1$target_gene<-gsub(" ", "", regulon1$target_gene)

regulon1 <- regulon1 %>%
  group_by(X) %>%
  summarise(all_target_gene = paste(unique(target_gene), collapse = ","))

regulon1 <- regulon1 %>%
  mutate(target_gene_count = lengths(unique(strsplit(all_target_gene, ","))))
###################### save the result
write.csv(regulon1,file = "Ex-L2_4_TF_target_genes.csv",quote = F)


###################### plot
regulon1$color<-c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87')
regulon1<-regulon1[order(regulon1$target_gene_count,decreasing = T),]
regulon1$X<-factor(regulon1$X,levels = regulon1$X)
pdf("Ex-L2_4_TF_target_genes.pdf",height = 4.5,width = 5)
ggplot(regulon1,aes(X,target_gene_count))+
  geom_col(aes(fill=X),width = 0.6)+
  scale_fill_manual(values = regulon1$color)+  
  geom_text(aes(label=target_gene_count),size = 3,vjust=-0.5) + 
  theme_bw()+
  xlab(" ")+
  ylab("Target gene numbers")+
  theme(axis.text.x = element_text(size = 10,color = 'black',vjust = 0.5,hjust = 1,angle = 90), 
        axis.text.y = element_text(size = 10,color = 'black'), 
        legend.text = element_text(size = 10,color = 'black'), 
        legend.title = element_text(size = 10,color = 'black'), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
dev.off()


##################### calculate the average target genes for each TF
regulon2<- regulon %>%
  mutate(gene_count = sapply(strsplit(target_gene, ","), length)) %>%
  group_by(X) %>%
  summarise(Average_Gene_Count = mean(gene_count))
write.csv(regulon2,file = "Ex-L2_4_average_TF_target_genes.csv",quote = F)

regulon2$color<-c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87')
regulon2<-regulon2[order(regulon2$Average_Gene_Count,decreasing = T),]
regulon2$X<-factor(regulon2$X,levels = regulon2$X)

pdf("Ex-L2_4_average_TF_target_genes.pdf",height = 4.5,width = 5)
ggplot(regulon2,aes(X,Average_Gene_Count))+
  geom_col(aes(fill=X),width = 0.6)+
  scale_fill_manual(values = regulon2$color)+  
  geom_text(aes(label=Average_Gene_Count),size = 3,vjust=-0.5) + 
  theme_bw()+
  xlab(" ")+
  ylab("Average target gene numbers")+
  theme(axis.text.x = element_text(size = 10,color = 'black',vjust = 0.5,hjust = 1,angle = 90), 
        axis.text.y = element_text(size = 10,color = 'black'), 
        legend.text = element_text(size = 10,color = 'black'), 
        legend.title = element_text(size = 10,color = 'black'), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
dev.off()








