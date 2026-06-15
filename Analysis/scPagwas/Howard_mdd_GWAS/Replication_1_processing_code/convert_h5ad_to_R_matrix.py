# -*- coding: utf-8 -*-
import scanpy as sc
import pandas as pd
import numpy as np
import seaborn as sns
import pyarrow as pa
import os
import pyarrow.feather as feather
from scipy.sparse import coo_matrix
import scipy.sparse
import csv

#### set the filepath
BaseDirectory = '/share2/pub/chenchg/chenchg/SingleCell/Brain/Human_million_Brain_singlecell/300w_data/'
i="merged_3M_singlecell_data_random_1M"
results_file="/share2/pub/chenchg/chenchg/SingleCell/Brain/Human_million_Brain_singlecell/300w_data/merged_3M_singlecell_data_random_1M.h5ad"

#### import the h5ad file
adata=sc.read_h5ad(results_file)

#### deal with matrix
if scipy.sparse.issparse(adata.X):
    dense_matrix = adata.X.toarray()
else:
    dense_matrix = adata.X
#### 转置矩阵    
dense_matrix = np.transpose(dense_matrix)
df = pd.DataFrame(dense_matrix)

#### 设置矩阵的行名和列名
obs=adata.obs
df.index = adata.var_names
df.columns=obs.index

#### 保存矩阵
out_file=os.path.join(BaseDirectory,i+'_matrix.feather')
feather.write_feather(df, out_file)

#### 保存df前10行
#df_first_10 = df.head(10)
#out_file_first_10 = os.path.join(BaseDirectory, i + '_matrix_first_10.feather')
#feather.write_feather(df_first_10, out_file_first_10)

#### 保存metadata
obs=adata.obs
obs_file=os.path.join(BaseDirectory,i+'_obs.csv')
obs.to_csv(obs_file, index=False, header=True)

#### 保存gene注释信息
var=adata.var
var_file=os.path.join(BaseDirectory,i+'_var_names.csv')
var=pd.DataFrame(var)
var.to_csv(var_file, index=True, header=True)
