# Single-cell genetic mapping of MDD-associated brain cell types

Analysis scripts, intermediate outputs and reproducibility notebooks accompanying:

> **Ma *et al.*** Single-cell genetic mapping links motor-cortical PVALB⁺ inhibitory circuitry to major depressive disorder. *under consideration*, 2026.

# scDepBrain
a online-resource for querying MDD-relevant cell types and subpopulations

🌐 **Interactive web resource:** [https://scdepbrain.su-lab.org/](https://scdepbrain.su-lab.org/) (>5 M cells / 34 brain single-cell and single-nucleus studies) — see also [Ma-brain-group/scDepBrain](https://github.com/Ma-brain-group/scDepBrain).

---

## Overview

This repository contains the full analysis pipeline used to integrate multi-method
GWAS cell-type prioritization with single-cell transcriptomics across three
independent MDD GWAS resources, and to validate the genetic-prioritization
results in an external case-control single-cell cohort.

## Workflow
![Study design](https://github.com/Ma-brain-group/scDepBrain/blob/main/Figures/Figure%201.png)




The project implements four orthogonal GWAS cell-type enrichment methods
(**scPagwas**, **MAGMA-CellTyping**, partitioned **LDSC-SEG** and **scDRS**)
applied to:

- **Howard *et al.* 2019** European-ancestry MDD GWAS (*n* cases = 246,363, *n* controls = 561,190)
- **PGC MDD 2025 European-ancestry** GWAS (*n* cases = 412,305, *n* controls = 1,588,397; *N*_eff = 1,152,656)
- **PGC MDD 2025 trans-ancestry** GWAS (*n* cases = 537,363, *n* controls = 2,061,567; *N*_eff = 1,487,075)

across an integrated discovery atlas of **309,983 nuclei from 109 donors** spanning
**14 anatomical brain regions** and **20 neuronal subtypes**, plus a replication
dataset of 160,222 nuclei from **GSE213982** (37 MDD cases vs 35 controls).

### Highlights

- **PVALB⁺ inhibitory neurons** and **L2/4 cortical excitatory neurons** are
  identified as the convergent MDD-associated populations across four orthogonal
  methods and three independent GWAS.
- **Primary motor cortex (M1C)** shows the highest single-cell MDD trait-relevance
  score, an effect that persists within the PVALB⁺ inhibitory neuron
  subpopulation (i.e. not a cell-composition artefact).
- **CNNM2 mRNA upregulation** in MDD PVALB⁺ inhibitory neurons is independently
  replicated by donor-level pseudobulk DESeq2 in GSE213982
  (log₂FC = +0.21, FDR = 0.018, *n* = 37 cases vs 35 controls).

---

---

## Methods

### Cell-type-level GWAS prioritization (4 orthogonal methods)

| Method | Output | Key parameters | Reference |
|---|---|---|---|
| **scPagwas** | per-cell TRS + cell-type *P* (`Merge_celltype_p`) | MAF ≥ 0.1, gene window ±10 kb, n_topgenes = 1,000 | Ma *et al.* 2023 *Cell Genomics* |
| **MAGMA-CellTyping** | per-cell-type *P* (gene-level covariate) | gene window ±50 kb, extended-MHC excluded, SNP-wise mean model, 1000G EUR LD | Bryois *et al.* 2020 *Nat Genet*, de Leeuw 2015 |
| **partitioned LDSC-SEG** | per-cell-type one-sided coefficient *P* | top 10% specific genes, ±100 kb window, baseline-LD v2.2 + all-gene control, 1-cM LD window, HapMap3 SNPs | Finucane *et al.* 2018 *Nat Genet* |
| **scDRS** | per-cell `norm_score` + Monte-Carlo cell-type `assoc_mcp` | top 1,000 MAGMA-z-weighted genes (cap ±10), 1,000 matched control gene sets | Zhang *et al.* 2022 *Nat Genet* |

### Downstream genetic-risk-gene analyses

- **High-confidence gene-set analysis** — 308 PGC-prioritized high-confidence MDD
  genes scored per-cell via Seurat `AddModuleScore` and as a weighted scDRS gene
  set; cell-type EWCE bootstrap (10,000 reps).
- **Brain-region stratification** — 14 anatomical brain regions (M1C, MTG, S1C,
  A1C, PFC, V1C, CTX, DFC, FC, LA, ACC, SN, CN, CER); validated within neuronal
  subpopulations.
- **Case-control replication (GSE213982)** — donor-level pseudobulk DESeq2
  (Squair *et al.* 2021 *Nat Commun* recommended workflow), *n* = 37 MDD vs 35
  control donors, accounting for sex.

---

## Data resources

### Single-cell input

- **Discovery integrated atlas**: 309,983 nuclei from 109 donors across 14 brain
  regions and 20 neuronal subtypes (see manuscript Methods).
- **Replication**: GSE213982 — 160,222 nuclei from 37 MDD cases and 35
  unaffected controls (Maitra *et al.* 2023 *Nat Commun*).

### GWAS resources

| GWAS | Cases | Controls | *N*_eff | SNPs | Genome | Reference |
|---|---|---|---|---|---|---|
| Howard 2019 (HW) | 246,363 | 561,190 | — | 8,481,297 | GRCh37 | Howard *et al.* 2019 *Nat Neurosci* |
| PGC MDD 2025 EUR | 412,305 | 1,588,397 | 1,152,656 | 7,363,302 | GRCh37 | Adams *et al.* 2024 *Cell* |
| PGC MDD 2025 trans-ancestry (DIV) | 537,363 | 2,061,567 | 1,487,075 | 5,918,521 | GRCh37 | Adams *et al.* 2024 *Cell* |

LDSC-SEG was not applied to the trans-ancestry GWAS owing to LD-structure
mismatch with the 1000 Genomes Phase 3 European reference panel.

---

## Reproducibility

### Software environment

```
R         >= 4.2
Python    >= 3.10 (for scDRS + anndata + scanpy)

R packages : Seurat (>= 5.0), scPagwas, DESeq2, data.table, dplyr, tidyr,
             ggplot2, ggpubr, ggrepel, patchwork, Matrix, SeuratDisk
Python     : scdrs (>= 1.0), anndata, scanpy, pandas, numpy, scipy
External   : MAGMA v1.10, LDSC (Python 2.7 env), PLINK v1.9, 1000G Phase 3
             EUR reference, baseline-LD v2.2, HapMap3 SNP list
```

### Reproduce the main figures

Each subfolder under `github_notebooks/` is self-contained (scripts +
intermediate CSV + final figure). Open the notebook directly on GitHub or run
locally with an R Jupyter kernel:

```bash
git clone <this-repo>
cd <repo>
# Open any notebook in JupyterLab with R kernel
jupyter lab github_notebooks/01_GWAS_4method_integration/
```

For the upstream computations (MAGMA, LDSC, scPagwas, scDRS) the original SLURM
scripts are kept under the corresponding `results_*` folder.


---


### scDepBrain provides an interactive resource for MDD-linked brain cell types and states
To facilitate exploration and reuse of the integrated single-cell atlas and MDD genetic-prioritization results, we developed the single-cell Depression Brain Map (scDepBrain,https://scdepbrain.su-lab.org/), a web-accessible resource for querying MDD-associated cellular programs across brain contexts. scDepBrain organizes more than five million cells and nuclei from 34 published brain single-cell and single-nucleus studies into searchable dataset entries with metadata including species, brain region, developmental stage, disease status, data source, accession number, cell count and omics type. The platform provides dataset-level visualization and analysis modules, including UMAP exploration, selectable gene-expression display, marker-gene inspection, cellular-composition comparison, differential-expression analysis and GO/KEGG pathway enrichment. scDepBrain further enables interactive exploration of genetic-prioritization results, including single-cell TRSs and cell-type-level MDD association results derived from scPagwas. Together, scDepBrain converts the integrated single-cell and genetic analyses into a public, reusable platform for cross-dataset comparison, hypothesis generation and future investigation of the cellular and molecular architecture of MDD.

## scDepBrain workflow
![Workflow](https://github.com/Ma-brain-group/scDepBrain/blob/main/Figures/scDepBrain_workflow.png)





# Citation
If you use this code or the accompanying scDepBrain resource, please cite:

> Ma et al. ***Single-cell genetic mapping links motor-cortical PVALB⁺ inhibitory circuitry to major depressive disorder***, 2026



## License

Code is released under the **MIT License**. Single-cell discovery data are
available from the corresponding author on reasonable request; GSE213982 is
publicly available from GEO. PGC MDD 2025 summary statistics are available
through the PGC consortium (https://pgc.unc.edu/).

## Contact

For questions about the code or analyses, please open an [issue](../../issues)
or contact the corresponding author (see paper).

🌐 **Resource:** [scdepbrain.su-lab.org](https://scdepbrain.su-lab.org/)
📂 **Resource code:** [Ma-brain-group/scDepBrain](https://github.com/Ma-brain-group/scDepBrain)






