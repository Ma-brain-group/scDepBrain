# PGC MDD 2025 European-ancestry — scPagwas pipeline

Scripts to run in numbered order:
```
01_prep_gwas.sh                    # parse PGC EUR tsv.gz → scPagwas format
02_run_scPagwas.R + .sbatch        # initial Pagwas_main reference run
sol1_05_prep_gwas.R + .sbatch      # prepare sparse pathway inputs
sol1_06_chunk_regression.R + .sbatch  # SLURM array — one chunk per task
sol1_07_finalize_v2.R + .sbatch    # merge chunks + Get_CorrectBg_p (patched)
sol1_08_celltype.R + .sbatch       # cell-type & subtype Merge_celltype_p
```

Results:
- `cell_type_pvalues_v2.csv` — 8 major cell types (Merge_celltype_p, final)
- `cell_type_pvalues_v2_subtype.csv` — 20 subtypes
- `cell_type_pvalues_celltype_method.csv` — 8 major × Bootstrap (sensitivity)
