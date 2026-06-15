
# Step 1 S-PrediXcan_SMultiXcan.sh
```sh
#!/usr/bin/env bash
#PBS -N Spredixcan_Smultixcan
#PBS -q workq
#PBS -l ncpus=10
#PBS -l mem=150gb

export DATA=/share2/pub/chenchg/chenchg/SingleCell/Brain/GWAS_data/Spredixcan_Smultixcan/Reference_data/
export GWAS_TOOLS=/share2/pub/chenchg/chenchg/software/summary-gwas-imputation-master/src
export METAXCAN=/share2/pub/chenchg/chenchg/software/MetaXcan-master/software
#export DIS=AD
export GWAS_DATA=/share2/pub/chenchg/chenchg/SingleCell/Brain/GWAS_data/Spredixcan_Smultixcan/Spredixcan_Smultixcan_input
export OUTPUT=/share2/pub/chenchg/chenchg/SingleCell/Brain/GWAS_data/Spredixcan_Smultixcan/Spredixcan_Smultixcan_results

mkdir $OUTPUT
mkdir $OUTPUT/harmonized_gwas
mkdir $OUTPUT/summary_imputation_1000G
mkdir $OUTPUT/processed_summary_imputation_1000G
mkdir $OUTPUT/spredixcan
mkdir $OUTPUT/spredixcan/eqtl
mkdir $OUTPUT/spredixcan/eqtl/mashr
mkdir $OUTPUT/smultixcan
mkdir $OUTPUT/smultixcan/eqtl

####### 第一步，将GWAS版本从HG37转为HG38
####### GWAS结果文件的染色体用数值表示，不要用chr+数值的方式
####### GWAS文件格式必须是gz压缩格式，分隔符为tab
source /share2/pub/chenchg/chenchg/miniconda3/bin/activate gwasimputation

python $GWAS_TOOLS/gwas_parsing.py \
-gwas_file $GWAS_DATA/Spredixcan_Smultixcan_input.txt.gz \
-liftover $DATA/data/liftover/hg19ToHg38.over.chain.gz \
-snp_reference_metadata $DATA/data/reference_panel_1000G/variant_metadata.txt.gz METADATA \
-output_column_map ID variant_id \
-output_column_map REF non_effect_allele \
-output_column_map ALT effect_allele \
-output_column_map ES effect_size \
-output_column_map  rawP pvalue \
-output_column_map CHROM chromosome \
--chromosome_format \
-output_column_map POS position \
--insert_value sample_size 500199 --insert_value n_cases 170756 \
-output_order variant_id panel_variant_id chromosome position effect_allele non_effect_allele frequency pvalue zscore effect_size standard_error sample_size n_cases \
-output $OUTPUT/harmonized_gwas/meta_GWAS_harmonized.txt.gz


######## 第二步，对gwas结果文件进行imputation
######## 循环对22条染色体，10个批次进行imputation
for chr in {1..22}; do
for batch in {0..9}; do
python $GWAS_TOOLS/gwas_summary_imputation.py \
-by_region_file $DATA/data/eur_ld.bed.gz \
-gwas_file $OUTPUT/harmonized_gwas/meta_GWAS_harmonized.txt.gz \
-parquet_genotype $DATA/data/reference_panel_1000G/chr${chr}.variants.parquet \
-parquet_genotype_metadata $DATA/data/reference_panel_1000G/variant_metadata.parquet \
-window 100000 \
-parsimony 7 \
-chromosome ${chr} \
-regularization 0.1 \
-frequency_filter 0.01 \
-sub_batches 10 \
-sub_batch ${batch} \
--standardise_dosages \
-output $OUTPUT/summary_imputation_1000G/GWAS_chr${chr}_sb${batch}_reg0.1_ff0.01_by_region.txt.gz
done
done


######## 第三步，合并220个批次的数据
python $GWAS_TOOLS/gwas_summary_imputation_postprocess.py \
-gwas_file $OUTPUT/harmonized_gwas/meta_GWAS_harmonized.txt.gz \
-folder $OUTPUT/summary_imputation_1000G \
-pattern "GWAS_*" \
-parsimony 7 \
-output $OUTPUT/processed_summary_imputation_1000G/imputed_GWAS.txt.gz
conda deactivate


######## 第四步，使用单个组织预测基因表达（这里只用13个脑区的模型数据进行预测）
######## S-PrediXcan mashr eqtl
source /share2/pub/chenchg/chenchg/miniconda3/bin/activate MetaXcan

for db in `ls $DATA/eqtl/mashr_brain/*db`; do
python $METAXCAN/SPrediXcan.py \
--gwas_file  $OUTPUT/processed_summary_imputation_1000G/imputed_GWAS.txt.gz \
--snp_column panel_variant_id --effect_allele_column effect_allele --non_effect_allele_column non_effect_allele --zscore_column zscore \
--model_db_path ${db} \
--covariance ${db%.db}.txt.gz \
--keep_non_rsid --additional_output --model_db_snp_key varID \
--throw \
--output_file $OUTPUT/spredixcan/eqtl/mashr/GWAS_${db##*/}.csv
done
  

######## 第四步，对单个组织预测的结果进行meta分析（用13个脑区的结果进行meta分析）
######## S-MultiXcan eqtl
python $METAXCAN/SMulTiXcan.py \
--models_folder $DATA/eqtl/mashr_brain \
--models_name_pattern "mashr_(.*).db" \
--snp_covariance $DATA/gtex_v8_expression_mashr_snp_smultixcan_covariance.txt.gz \
--metaxcan_folder $OUTPUT/spredixcan/eqtl/mashr/ \
--metaxcan_filter "GWAS_mashr_(.*).db.csv" \
--metaxcan_file_name_parse_pattern "(.*)_mashr_(.*).db.csv" \
--gwas_file  $OUTPUT/processed_summary_imputation_1000G/imputed_GWAS.txt.gz \
--snp_column panel_variant_id --effect_allele_column effect_allele --non_effect_allele_column non_effect_allele --zscore_column zscore --keep_non_rsid --model_db_snp_key varID \
--cutoff_condition_number 30 \
--verbosity 7 \
--throw \
--output $OUTPUT/smultixcan/eqtl/GWAS_mashr_smultixcan_eqtl.txt

```

# Step 2 S-PrediXcan_SMultiXcan.R
```R
############# Spredixcan_Smultixcan input ################
library(readr)
library(vcfR)
library(dplyr)
library(stringr)

setwd("C:/Users/Lenovo/Desktop/")

############# GWAS数据读入
GWAS_raw <-read.vcfR("ieu-b-102.vcf")
GWAS <- cbind(GWAS_raw@fix,GWAS_raw@gt) %>% as.data.frame()
colnames(GWAS )<-c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","ieu")

############ 输入文件应有的列
############ CHROM POS ID REF ALT ES SE rawP AF
a<-GWAS[1,9]##ES:SE:LP:SS:ID
b<-data.frame(GWAS[,10])##-0.0408012:0.019:1.49377:55374:rs3094315
a<-t(apply(data.frame(a),1,function(x){unlist(strsplit(x,":"))}))
n<-length(unlist(strsplit(b[1,],":")))

########### str_split_fixed按照指定分隔符拆分，返回矩阵格式，n为返回的列数,拆分ES:SE:LP:SS:ID列
b<-t(apply(b,1,function(x){str_split_fixed(x,":",n=n)}))
b<-as.data.frame(b)
colnames(b)<-a

########### 合并一下文件
sumstats<-cbind(GWAS,b)

########### 提取出必须列CHROM POS ID REF ALT ES SE rawP AF（maf）
########### CHROM只需数字即可
########### 若GWAS本身就不存在maf，给一个固定的maf值0.1
LP<-sumstats$LP
LP<-as.numeric(LP)
if ("AF"%in%colnames(sumstats)) {
  sumstats1<-data.frame(CHROM=sumstats$CHROM,POS=sumstats$POS,ID=sumstats$ID,REF=sumstats$REF,ALT=sumstats$ALT,
                        ES=sumstats$ES,SE=sumstats$SE,rawP=10^-LP,AF=sumstats$AF)
}else{
  sumstats1<-data.frame(CHROM=sumstats$CHROM,POS=sumstats$POS,ID=sumstats$ID,REF=sumstats$REF,ALT=sumstats$ALT,
                        ES=sumstats$ES,SE=sumstats$SE,rawP=10^LP,AF=0.1)
}

########### 去除MHC区域的位点
sumstats1<- sumstats1[!( sumstats1$CHROM=="6"& sumstats1$POS > 25000000 &  sumstats1$POS < 34000000),] 

########### 输出
sumstats1<-na.omit(sumstats1)
write.table(sumstats1,file = "Spredixcan_Smultixcan_input.txt",row.names = F,sep = "\t",quote = F)
```


