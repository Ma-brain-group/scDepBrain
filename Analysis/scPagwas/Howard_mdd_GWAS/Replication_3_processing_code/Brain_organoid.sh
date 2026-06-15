#PBS -N Brain_organoid
#PBS -q fat
#PBS -l nodes=fat03
#PBS -l mem=450gb
#PBS -l ncpus=5
cd /share2/pub/chenchg/chenchg/SingleCell/Brain/organoid/Sample_merge/
source activate R4.2.2
Rscript Brain_organoid.R
