#PBS -N convert_h5ad_to_R_matrix
#PBS -q fat
#PBS -l nodes=fat03
#PBS -l ncpus=10
#PBS -l mem=400gb
cd /share2/pub/chenchg/chenchg/SingleCell/Brain/Human_million_Brain_singlecell/300w_data/
source activate python
python convert_h5ad_to_R_matrix.py
