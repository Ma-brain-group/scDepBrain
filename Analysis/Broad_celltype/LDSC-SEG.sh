#PBS -N ldsc
#PBS -q workq
#PBS -l mem=150gb
#PBS -l ncpus=3

##ldsc���д���

##����ldsc����
source activate ldsc

##������ϴ
for i in 2013_NG_15_risk_loci_MDD 2020_MP_GWAS_MDD 2022_uk_Depressive_symptoms broad_depression ICD_10-coded_MDD probable_MDD uk_self_depression ieu_b_102
do
python /share2/pub/chenchg/chenchg/software/ldsc/munge_sumstats.py \
--sumstats /share2/pub/chenchg/chenchg/SingleCell/Brain/GWAS_data/LDSC_SEG/ldsc_input_summary/$i.sumstats \
--merge-alleles /share2/pub/chenchg/chenchg/SingleCell/Brain/GWAS_data/LDSC_SEG/w_hm3.snplist \
--chunksize 500000 \
--out /share2/pub/chenchg/chenchg/SingleCell/Brain/GWAS_data/LDSC_SEG/ldsc_sumstats/$i 
done
##--N 174519 �����ļ�������������Ϣ


##scRNA������ȡ������������LD����������partitioned LD �����ع��Լ��ϸ��������
##���ⲽ֮ǰ��������ÿ��ϸ�������������top10%�Ļ�����Ϊ��ϸ����������������

##��һ��������ע���ļ�annot file
for i in Astrocytes Endothelial.cells Excitatory.neurons Inhibitory.neurons Microglia Oligodendrocytes OPCs Purkinje.neurons control
do
  for j in $(seq 1 22)
    do
      python /share2/pub/chenchg/chenchg/software/ldsc/make_annot.py \
       --gene-set-file /share/pub/qiuf/brain/01-data/GWAS/LDSC_anno1/Bed/$i.bed \
       --gene-coord-file /share/pub/qiuf/brain/01-data/GWAS/LDSC_anno2/test/ENSG_coord.txt \
       --windowsize 100000 \
       --bimfile /share2/pub/chenchg/chenchg/SingleCell/Brain/GWAS_data/LDSC_SEG/1000G_EUR_Phase3_plink/1000G.EUR.QC.$j.bim \
       --annot-file /share2/pub/chenchg/chenchg/SingleCell/Brain/GWAS_data/LDSC_SEG/annot/$i.$j.annot.gz
    done
done

##�ڶ�����ʹ��ע���ļ�����LD����������.M�ļ�
for i in Astrocytes Endothelial.cells Excitatory.neurons Inhibitory.neurons Microglia Oligodendrocytes OPCs Purkinje.neurons control
do
  for j in $(seq 1 22)
    do
      python /share2/pub/chenchg/chenchg/software/ldsc/ldsc.py --l2 --bfile /share2/pub/chenchg/chenchg/SingleCell/Brain/GWAS_data/LDSC_SEG/1000G_EUR_Phase3_plink/1000G.EUR.QC.$j --ld-wind-cm 1 --annot /share2/pub/chenchg/chenchg/SingleCell/Brain/GWAS_data/LDSC_SEG/annot/$i.$j.annot.gz --thin-annot --out /share2/pub/chenchg/chenchg/SingleCell/Brain/GWAS_data/LDSC_SEG/annot/$i.$j --print-snps /share2/pub/chenchg/chenchg/SingleCell/Brain/GWAS_data/LDSC_SEG/hapmap3_snps/hm.$j.snp
    done
done


##LDSC-SEG
for i in 2013_NG_15_risk_loci_MDD 2020_MP_GWAS_MDD 2022_uk_Depressive_symptoms broad_depression ICD_10-coded_MDD probable_MDD uk_self_depression ieu_b_102
do
   python /share2/pub/chenchg/chenchg/software/ldsc/ldsc.py \
    --h2-cts /share2/pub/chenchg/chenchg/SingleCell/Brain/GWAS_data/LDSC_SEG/ldsc_sumstats/$i.sumstats.gz \
    --ref-ld-chr /share2/pub/chenchg/chenchg/SingleCell/Brain/GWAS_data/LDSC_SEG/1000G_EUR_Phase3_baseline/baseline. \
    --out /share2/pub/chenchg/chenchg/SingleCell/Brain/GWAS_data/LDSC_SEG/ldsc_result/$i \
    --ref-ld-chr-cts /share2/pub/chenchg/chenchg/SingleCell/Brain/GWAS_data/LDSC_SEG/ldcts.ldcts \
    --w-ld-chr /share2/pub/chenchg/chenchg/SingleCell/Brain/GWAS_data/LDSC_SEG/weights_hm3_no_hla/weights.
done























