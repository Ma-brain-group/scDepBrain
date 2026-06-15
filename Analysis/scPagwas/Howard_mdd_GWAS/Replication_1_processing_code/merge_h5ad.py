# -*- coding: utf-8 -*-
import scanpy as sc
import numpy as np
import pandas as pd
import os
## set the filepath
os.chdir("/share2/pub/chenchg/chenchg/SingleCell/Brain/Human_million_Brain_singlecell/300w_data/")
Nonneurons = sc.read_h5ad("Nonneurons.h5ad")
Neurons = sc.read_h5ad("Neurons.h5ad")

## calculate the genes numbers and cell numbers
print(Nonneurons.shape)
print(Neurons.shape)

## compare the gene
genes1 = Nonneurons.var_names
genes2 = Neurons.var_names
if genes1.equals(genes2):
    print("equal")
else:
    print("inequal")

## merge the data and save the merged the data
adata_merged = Neurons.concatenate(Nonneurons, batch_key="batch", batch_categories=["Neurons", "Nonneurons"])
adata_merged.obsm['X_tsne'] = adata_merged.obsm['X_tSNE']
adata_merged.obsm['X_umap'] = adata_merged.obsm['X_UMAP']
#adata_merged.write('merged_3M_singlecell_data.h5ad')

## tsne plot
#sc.pl.tsne(adata_merged, color=['supercluster_term'], save='_celltypes_plot.png')
#sc.pl.tsne(adata_merged, color=['batch'], save='_batch_plot.png')

## randomly get 50w cells according to the origion celltypes percent
## Set the random seed for reproducibility
np.random.seed(123)

# get the celltypes information 
cell_types = adata_merged.obs['supercluster_term']
# calculate the cell counts for each celltypes
cell_type_counts = cell_types.value_counts()
total_cells_to_sample = 500000
cell_type_sample_counts = (cell_type_counts / cell_type_counts.sum() * total_cells_to_sample).astype(int)
# construct a list for cell index
sampled_indices = []
for cell_type, sample_count in cell_type_sample_counts.items():
    # get the celltype item
    indices = adata_merged.obs.index[cell_types == cell_type]
    # random select cells by origion celltype percent for each celltype
    sampled_indices.extend(np.random.choice(indices, size=sample_count, replace=False))
adata_sampled = adata_merged[sampled_indices, :].copy()

## print the data size
print(f"抽样后的数据集形状: {adata_sampled.shape}")

## save the randomly selected data
adata_sampled.write('merged_3M_singlecell_data_random_1M.h5ad')













