#PBS -N merge_h5ad
#PBS -q fat
#PBS -l nodes=fat03
#PBS -l ncpus=10
#PBS -l mem=400gb
cd /share2/pub/chenchg/chenchg/SingleCell/Brain/Human_million_Brain_singlecell/300w_data/
source activate python
python merge_h5ad.py
