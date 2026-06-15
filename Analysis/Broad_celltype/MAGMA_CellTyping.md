
# Step 1 magma.sh
```sh
#PBS -N magma
#PBS -q workq
#PBS -l mem=150gb
#PBS -l ncpus=2

##magma��������
cd /share2/pub/chenchg/chenchg/SingleCell/Brain/GWAS_data/magma/

##step1��ע�ͣ���GWAS�ϵ�snpע�͵�������
##--snp-loc��magma�����ļ�2������ע�͵ģ��������У�SNP,CHR,POS��Ϣ
##--gene-loc��������magma���صĻ���ע����Ϣ��ע�������汾
mkdir gene_annotation
cd gene_annotation
for i in 2013_NG_15_risk_loci_MDD 2020_MP_GWAS_MDD 2022_uk_Depressive_symptoms broad_depression ICD_10-coded_MDD probable_MDD uk_self_depression ieu_b_102
do
magma --annotate window=50,50 \
      --snp-loc /share2/pub/chenchg/chenchg/SingleCell/Brain/GWAS_data/magma/magma_input2/"$i"_magma_input2.txt \  
      --gene-loc /share2/pub/chenchg/chenchg/SingleCell/Brain/GWAS_data/magma/magma_data/NCBI/NCBI37.3.gene.loc.extendedMHCexcluded \
      --out $i.int.annotated_50kbup_50_down
done

##step2��gene-based��������
##--pval��ÿһ��SNP��Pֵ��ncol/N��������������Ϣ
##Magma�ҵ����Ŵ���صĻ��򣬼�����һ������
for i in 2013_NG_15_risk_loci_MDD 2020_MP_GWAS_MDD 2022_uk_Depressive_symptoms broad_depression ICD_10-coded_MDD probable_MDD uk_self_depression ieu_b_102
do
magma --bfile /share/pub/dengcy/Singlecell/COVID19/MAGMA/g1000_eur/g1000_eur \  
      --pval /share2/pub/chenchg/chenchg/SingleCell/Brain/GWAS_data/magma/magma_input1/"$i"_magma_input1.txt ncol=3 \  
      --gene-annot /share2/pub/chenchg/chenchg/SingleCell/Brain/GWAS_data/magma/gene_annotation/"$i".int.annotated_50kbup_50_down.genes.annot \  
      --out /share2/pub/chenchg/chenchg/SingleCell/Brain/GWAS_data/magma/gene_based/"$i"
done

##step3��geneset���������Կ���ͨ·����������
##--gene-results����һ�����ɵ�raw�ļ�
##--set-annot��������MsigDB���ݿ����ص�ͨ·��Ϣ��Ҳ�����Զ������
##��MAGMA�ҵ��Ļ��򸻼���ϸ����������Ļ�����
for i in 2013_NG_15_risk_loci_MDD 2020_MP_GWAS_MDD 2022_uk_Depressive_symptoms broad_depression ICD_10-coded_MDD probable_MDD uk_self_depression ieu_b_102
do
magma --gene-results /share2/pub/chenchg/chenchg/SingleCell/Brain/GWAS_data/magma/gene_based/"$i".genes.raw --set-annot /share/pub/qiuf/brain/01-data/GWAS/MAGMA_anno1/top10.txt --out /share2/pub/chenchg/chenchg/SingleCell/Brain/GWAS_data/magma/geneset/"$i"
done


```

# Step 2 MAGMA based gene results
```R
############## magma gene #################
setwd("D:/Project/SingleCell_MDD/GWAS_data/magma/")
magma<-read.table("ieu_b_102.genes.out",sep = "",header = T)

############## P值做一下fdr矫正 ##############
magma$fdr<-p.adjust(magma$P,method = "fdr")
############# 905个fdr显著的基因
magma<-magma[magma$fdr<0.05,]

############# 将基因ID转为symbol ###############
library(clusterProfiler)
library(org.Hs.eg.db)

gene<-bitr(magma$GENE,
           fromType="ENTREZID",
           toType="SYMBOL",
           OrgDb = "org.Hs.eg.db")
colnames(gene)[1]<-"GENE"
gene1<-merge(gene,magma,by="GENE")
gene1<-gene1[order(gene1$fdr,decreasing = F),]
################# 保存结果 ######################
write.csv(gene1,file="magma_gene_fdr_0.05.csv",row.names = T,quote=F)


############### Magma CNNM2 top SNP #################
library(data.table)
SNP<-read.table("ieu_b_102.int.annotated_50kbup_50_down.genes.annot",fill = T)

############### 提取出CNNM2对应的snp ##################
gene1[gene1$SYMBOL=="CNNM2",]
SNP1<-SNP[SNP[,1]=="54805",]
SNP1<-na.omit(SNP1)
SNP1<-SNP1[,3:32]

############## 在GWAS SUMMARY 数据中提取出对应的SNP #################
GWAS_summary<-read.table("ieu_b_102_magma_input1.txt",sep = "\t",header = T)
GWAS_summary<-GWAS_summary[order(GWAS_summary$p),]
GWAS_summary1<-GWAS_summary[match(SNP1[1,],GWAS_summary$SNP),]
GWAS_summary1<-GWAS_summary1[order(GWAS_summary1$p,decreasing = F),]

############# 保存
write.csv(GWAS_summary1,file = "MAGMA_CNNM2_to_SNP.csv",quote = F)

```















