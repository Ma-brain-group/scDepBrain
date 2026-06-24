# PGC MDD 2025 trans-ancestry — scPagwas pipeline

Scripts to run in numbered order:
```
01_prep_gwas.sh                              # parse PGC DIV tsv.gz → scPagwas format
02_run_scPagwas.R + .sbatch                  # initial Pagwas_main reference run
sol1_05_prep_gwas_paper.R + .sbatch          # paper-params (MAF 0.1, marg 10 kb, n_top 1000)
sol1_06_chunk_regression.R + .sbatch         # SLURM array — one chunk per task
sol1_07_finalize_paper.R + .sbatch           # merge + Get_CorrectBg_p (patched)
sol1_08_celltype.R + .sbatch                 # cell-type & subtype Merge_celltype_p
```

Results:
- `cell_type_pvalues_paperparams.csv` — 8 major
- `cell_type_pvalues_paperparams_subtype.csv` — 20 subtypes
- `cell_type_pvalues_celltype_method_paperparams.csv` — Bootstrap sensitivity

**Trans-ancestry caveat:** the strong oligodendrocyte signal from `Merge_celltype_p`
in DIV is **not replicated** by MAGMA-CellTyping, scDRS or LDSC-SEG and likely
reflects a bidirectional-aggregation artefact of the `Merge_celltype_p`
statistic in trans-ancestry meta-analyses (see manuscript Discussion).
