# scPagwas analysis — PGC MDD 2025 European (EUR) and trans-ancestry (DIV)

Final scPagwas pipeline scripts and per-cell-type association results for the two
PGC MDD 2025 meta-analysis GWAS, used as the **cross-cohort and cross-ancestry
replication** for the multi-method cell-type enrichment analysis in our MDD
manuscript.

## GWAS resources

| GWAS | n cases | n controls | N_eff | n SNPs | Genome | Reference |
|---|---|---|---|---|---|---|
| **PGC MDD 2025 European** (EUR) | 412,305 | 1,588,397 | 1,152,656 | 7,363,302 | GRCh37 |  PGC MDD. Cell 2025 |
| **PGC MDD 2025 trans-ancestry** (DIV) | 537,363 | 2,061,567 | 1,487,075 | 5,918,521 | GRCh37 | PGC MDD. Cell 2025 |

## Folder layout

```
results_scPagwas_DIV_EUR_20260623/
├── EUR/
│   ├── scripts/          # 7 R + sbatch (numbered 01 → 08)
│   └── results/          # cell-type level P values (8 broad + 20 subtype)
├── DIV/
│   ├── scripts/
│   └── results/
└── README.md             # this file
```

## Pipeline overview (`sol1` solution)

`Solution 1` is our memory-tractable variant of the scPagwas pipeline,
specifically adapted for the 24 GB / 310k-cell integrated discovery atlas. The
core scPagwas computations are split into per-cell chunks and finalized once
all chunks return. Bug-fix patches applied throughout (Seurat 5 compatibility +
sparse-aware row variance + `findInterval`-based empirical *P*).

| Step | Script | Function | Approx. runtime |
|---|---|---|---|
| 01 | `01_prep_gwas.sh` | Extract from PGC 2025 tsv.gz → scPagwas format (rsid + chr + pos + se + beta + P + N_eff) | 1 min |
| 02 | `02_run_scPagwas.R` + `.sbatch` | `Pagwas_main()` initial run (reference for parameters) | ~30 min |
| 05 | `sol1_05_prep_gwas.R` (EUR) / `sol1_05_prep_gwas_paper.R` (DIV) + `.sbatch` | scPagwas prep, set MAF / window / `n_topgenes` per paper convention | 5 min |
| 06 | `sol1_06_chunk_regression.R` + `.sbatch` (array of N chunks) | `Link_pathway_blocks_gwas` per chunk (~200–400 GB RAM, 2–12 h per chunk) | 4–24 h (parallel) |
| 07 | `sol1_07_finalize_v2.R` (EUR) / `sol1_07_finalize_paper.R` (DIV) + `.sbatch` | Merge chunks, run `Get_CorrectBg_p` with sparse-row-variance patch + `findInterval` empirical *P* (turns 10+ h step into 30 sec) | 1 h |
| 08 | `sol1_08_celltype.R` + `.sbatch` | Compute per-cell-type and per-subtype P values (`Merge_celltype_p` + Bootstrap) | 15 min |

## Key parameters (paper-style)

| Parameter | Value | Notes |
|---|---|---|
| MAF filter | ≥ 0.1 | Stringent; matches Ma *et al.* 2023 Cell Genomics paper |
| Gene window | ±10 kb upstream / downstream | `marg = 10000` |
| `n_topgenes` (BG correction) | 1,000 | Per `Get_CorrectBg_p` |
| Features (TRS) | 2 (top + down) | No upTRS computed (would cost extra memory and is not in published paper) |
| BG control random sampling | `Nrandom = 5`, `Nselect = 200` | Default scPagwas vignette |
| Chunk size | ~62k cells per chunk | 5 chunks total |
| Single-cell input | `MDD_singlecell_data_reannotation.rds` (309,983 cells, 109 donors) | Integrated discovery atlas |

## Result files

### `EUR/results/`
| File | Content |
|---|---|
| `cell_type_pvalues_v2.csv` | 8 major cell types × Merge_celltype_p P values (`$anno`) |
| `cell_type_pvalues_v2_subtype.csv` | 20 subtypes × Merge_celltype_p P values (`$final_anno`) |
| `cell_type_pvalues_celltype_method.csv` | 8 major × Bootstrap (Boot_evaluate) — sensitivity test |

### `DIV/results/`
| File | Content |
|---|---|
| `cell_type_pvalues_paperparams.csv` | 8 major × Merge_celltype_p |
| `cell_type_pvalues_paperparams_subtype.csv` | 20 subtypes × Merge_celltype_p |
| `cell_type_pvalues_celltype_method_paperparams.csv` | 8 major × Bootstrap |

## Notes / caveats

1. **EUR vs DIV file-suffix difference.** EUR uses `_v2.csv` (v2 fix for `Get_CorrectBg_p` named-vector
   bug); DIV uses `_paperparams.csv` (same patches applied, suffix reflects the
   v2-derived "paper-params" finalize variant). The two variants implement identical
   compute logic; the suffix difference is historical.
2. **Trans-ancestry caveat.** scPagwas Merge_celltype_p in DIV showed a strong
   oligodendrocyte signal that is not replicated by MAGMA-CellTyping, scDRS or
   partitioned LDSC-SEG on the same data; see the manuscript Discussion for the
   interpretation that this reflects a bidirectional-aggregation artefact of
   `Merge_celltype_p` in trans-ancestry meta-analyses.
3. **EUR + DIV scPagwas TRSs** are also stored in the per-cell metadata columns
   of the saved scPagwas RDS objects (`scPagwas.TRS.Score`,
   `scPagwas.upTRS.Score`, `scPagwas.downTRS.Score`, `scPagwas.gPAS.score`); these
   were used for the per-cell scPagwas–scDRS correlation analyses described in
   the manuscript.

## Citation

If you reproduce or extend this pipeline, please cite:

- Ma, Y. *et al.* (2023). Polygenic regression uncovers trait-relevant cellular contexts through pathway activation transformation of single-cell RNA sequencing data. *Cell Genomics* 3, 100383.
-  PGC MDD (2025). Trans-ancestry genome-wide study of depression identifies 697 associations implicating cell types and pharmacotherapies. *Cell*, in press.

## Software environment

```
R >= 4.2
Packages: Seurat (>= 5.0), scPagwas, bigreadr, Matrix, data.table, dplyr
Patched functions (loaded inside scripts):
  - scPagwas::Get_CorrectBg_p           (sparse-aware row variance)
  - scPagwas::get_p_from_empi_null      (findInterval empirical P)
  - SeuratObject::.Deprecate            (no-op patch for Seurat 5 compat)
```

Single-cell input rds was the integrated MDD healthy brain discovery atlas
(`MDD_singlecell_data_reannotation.rds`, ~24 GB), available from the corresponding
author on reasonable request.
