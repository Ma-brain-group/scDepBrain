################### scPagwas input GWAS file ###################
library(readr)
library(vcfR)
library(dplyr)
library(stringr)

################### Linux passing parameters
#args <- commandArgs(T)
#path <- print(args[1])
#memo <- print(args[2])
#vcf <- print(args[3])
#output <- print(args[4])
#maf<- print(args[5])
#path<-"C:/Users/g/Desktop/UK_biobank_data_for_depression"
#memo<-"2013_NG_15_risk_loci_MDD"
#vcf<-"ieu-a-805.vcf"
setwd("D:/Project/SingleCell_MDD/GWAS_combine_singlecell_analysis/scPagwas_Height/")

#################### import the GWAS summary data
GWAS_raw <-read.vcfR("ukb-b-10787.vcf")
GWAS <- cbind(GWAS_raw@fix,GWAS_raw@gt) %>% as.data.frame()
colnames(GWAS )<-c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","ieu")

#################### subset the column for FORMAT
####ES:SE:LP:SS:ID
a<-GWAS[1,9]
####-0.0408012:0.019:1.49377:55374:rs3094315
b<-data.frame(GWAS[,10])
a<-t(apply(data.frame(a),1,function(x){unlist(strsplit(x,":"))}))
n<-length(unlist(strsplit(b[1,],":")))
####get the information for FORMAT
b<-t(apply(b,1,function(x){str_split_fixed(x,":",n=n)}))
b<-as.data.frame(b)
colnames(b)<-a

#################### merge the data
sumstats<-cbind(GWAS,b)
#################### get the column for chrom，rsid，pos，beta，se，maf
#################### if the AF dosen't exsit, maf=0.1
maf<-0.1
if ("AF"%in%colnames(sumstats)) {
  sumstats1<-data.frame(chrom=sumstats$CHROM,rsid=sumstats$ID,pos=sumstats$POS,beta=sumstats$ES,se=sumstats$SE,maf=sumstats$AF)
}else{
  sumstats1<-data.frame(chrom=sumstats$CHROM,rsid=sumstats$ID,pos=sumstats$POS,beta=sumstats$ES,se=sumstats$SE,maf=maf)
}

#################### remove the MHC information
sumstats1<- sumstats1[!( sumstats1$chrom=="6"& sumstats1$pos > 25000000 &  sumstats1$pos < 34000000),] 

#################### save the result
write.table(sumstats1,file = "Height_ukb-b-10787.pagwas.txt",row.names = F,sep = "\t",quote = F)













