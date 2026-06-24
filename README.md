# Single-cell genetic mapping of MDD-associated brain cell types

Analysis scripts, intermediate outputs and reproducibility notebooks accompanying:

> **Ma *et al.*** Single-cell genetic mapping links motor-cortical PVALB⁺ inhibitory circuitry to major depressive disorder. *Nature Neuroscience*, 2026.

🌐 **Interactive web resource:** [https://scdepbrain.su-lab.org/](https://scdepbrain.su-lab.org/) (>5 M cells / 34 brain single-cell and single-nucleus studies) — see also [Ma-brain-group/scDepBrain](https://github.com/Ma-brain-group/scDepBrain).

---

## Overview

This repository contains the full analysis pipeline used to integrate multi-method
GWAS cell-type prioritization with single-cell transcriptomics across three
independent MDD GWAS resources, and to validate the genetic-prioritization
results in an external case-control single-cell cohort.

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

## Repository structure

```
.
├── README.md                                   # this file
├── github_notebooks/                           # GitHub-ready reproducibility notebooks
│   ├── README.md
│   ├── 01_GWAS_4method_integration/            # scPagwas + MAGMA + LDSC + scDRS integration
│   ├── 02_PGC295_geneset_analysis/             # 308 high-confidence MDD gene-set analyses
│   ├── 03_CNNM2_replication/                   # pseudobulk DESeq2 replication
│   ├── 04_scDRS_region_validation/             # scDRS per-cell TRS by brain region
│   └── 05_brain_region_validation/             # cell-type-stratified region analysis
│
├── results_scPagwas_DIV_EUR_20260623/          # scPagwas pipeline (EUR + DIV)
│   ├── EUR/{scripts, results}
│   ├── DIV/{scripts, results}
│   └── README.md
│
├── results_MAGMA_PGC2025_20260617/             # MAGMA gene-based + MAGMA-CellTyping
│   ├── scripts/                                # 4 numbered .sh + .sbatch
│   ├── HW/, EUR/, DIV/                         # MAGMA .genes.out + .genes.raw per GWAS
│   ├── celltype_genesets/                      # top-10% specific gene lists
│   └── enrichment/                             # MAGMA gene-covar/set-annot outputs
│
├── results_LDSC_PGC2025_20260617/              # partitioned LDSC-SEG
│   ├── refs/GRCh37/                            # 1000G + baseline-LD v2.2 + weights
│   ├── annot/                                  # per-cell-type annotations
│   ├── sumstats/                               # munged sumstats
│   └── enrichment_v2/                          # S-LDSC .results
│
├── results_scDRS_20260618/                     # scDRS — per-cell trait-relevance scores
│   ├── gs/                                     # MAGMA z-weighted gene sets
│   ├── scores_1k/                              # full per-cell TRS (1,000 ctrl, 3 GWAS)
│   └── downstream_1k/                          # cell-type association P (assoc_mcp)
│
├── results_PGC295_geneset_20260618/            # PGC 308 high-confidence gene-set analyses
│   ├── gs/, scores/, downstream/               # MAGMA z-weighted scDRS
│   └── figures/                                # AddModuleScore + EWCE + scDRS dot-charts
│
├── results_CNNM2_replication_20260619/         # GSE213982 case-control replication
│   ├── scripts/                                # pseudobulk DESeq2 + per-cell Wilcoxon
│   └── figures/, data/                         # forest plot + per-donor box-plot stats
│
└── results_HW_v2_figures_20260617/             # all published figures + integrated CSVs
    ├── figures/                                # PDF + PNG (Fig 4 + Supp 25-32)
    └── data/                                   # per-figure source CSVs
```

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

### Workflow order

```
GWAS sumstats
   ├─→ [1] MAGMA gene-wise (±50 kb, ext-MHC excluded)
   │       └─→ [2a] MAGMA-CellTyping (gene-covar)
   │       └─→ [2b] scDRS gene-set (top 1000, z-weight)
   ├─→ [3] LDSC munge + partitioned LDSC-SEG (EUR only)
   └─→ [4] scPagwas (chunked solution-1, sol1_05 → sol1_08)
            │
            └─→ Integrated 4-method table → main Fig 4 + Supp 25-31
                 │
                 └─→ Brain-region stratification (scDRS per-cell TRS × 14 regions)
                 └─→ CNNM2 + NEGR1 marker correlation
                 └─→ GSE213982 pseudobulk DESeq2 replication
```

---

## Citation

If you use this code or the accompanying scDepBrain resource, please cite:

> Ma *et al.* Single-cell genetic mapping links motor-cortical PVALB⁺ inhibitory
> circuitry to major depressive disorder. *Nature Neuroscience*, 2026.

Please also cite the underlying method papers:

```
@article{scPagwas2023,
  author  = {Ma, Y. and others},
  title   = {Polygenic regression uncovers trait-relevant cellular contexts
             through pathway activation transformation of single-cell RNA sequencing data},
  journal = {Cell Genomics},
  year    = {2023}, volume = {3}, pages = {100383}
}
@article{MAGMA2015,
  author  = {de Leeuw, C. A. and others},
  title   = {MAGMA: Generalized gene-set analysis of GWAS data},
  journal = {PLoS Computational Biology}, year = {2015}
}
@article{LDSC-SEG2018,
  author  = {Finucane, H. K. and others},
  title   = {Heritability enrichment of specifically expressed genes
             identifies disease-relevant tissues and cell types},
  journal = {Nature Genetics}, year = {2018}, volume = {50}, pages = {621--629}
}
@article{scDRS2022,
  author  = {Zhang, M. J. and others},
  title   = {Polygenic enrichment distinguishes disease associations of
             individual cells in single-cell RNA-seq data},
  journal = {Nature Genetics}, year = {2022}, volume = {54}, pages = {1572--1580}
}
@article{Squair2021,
  author  = {Squair, J. W. and others},
  title   = {Confronting false discoveries in single-cell differential expression},
  journal = {Nature Communications}, year = {2021}
}
@article{Howard2019,
  author  = {Howard, D. M. and others},
  title   = {Genome-wide meta-analysis of depression identifies 102 independent
             variants and highlights the importance of the prefrontal brain regions},
  journal = {Nature Neuroscience}, year = {2019}, volume = {22}, pages = {343--352}
}
@article{PGC2024,
  author  = {Adams, M. J. and others},
  title   = {Genome-wide association study meta-analysis of major depressive disorder},
  journal = {Cell}, year = {2024}
}
```

---

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


# scDepBrain
a online-resource for querying MDD-relevant cell types and subpopulations

## Overview of current study
![Study design](https://github.com/Ma-brain-group/scDepBrain/blob/main/Figures/Figure%201.png)



### scDepBrain provides an interactive resource for MDD-linked brain cell types and states
To facilitate exploration and reuse of the integrated single-cell atlas and MDD genetic-prioritization results, we developed the single-cell Depression Brain Map (scDepBrain,https://scdepbrain.su-lab.org/), a web-accessible resource for querying MDD-associated cellular programs across brain contexts. scDepBrain organizes more than five million cells and nuclei from 34 published brain single-cell and single-nucleus studies into searchable dataset entries with metadata including species, brain region, developmental stage, disease status, data source, accession number, cell count and omics type. The platform provides dataset-level visualization and analysis modules, including UMAP exploration, selectable gene-expression display, marker-gene inspection, cellular-composition comparison, differential-expression analysis and GO/KEGG pathway enrichment. scDepBrain further enables interactive exploration of genetic-prioritization results, including single-cell TRSs and cell-type-level MDD association results derived from scPagwas. Together, scDepBrain converts the integrated single-cell and genetic analyses into a public, reusable platform for cross-dataset comparison, hypothesis generation and future investigation of the cellular and molecular architecture of MDD.

## scDepBrain workflow
![Workflow](https://github.com/Ma-brain-group/scDepBrain/blob/main/Figures/scDepBrain_workflow.png)





# Citation
If you use scRBP in your research, please cite:
> Ma et al. ***Single-cell genetic mapping links motor-cortical PVALB⁺ inhibitory circuitry to major depressive disorder***, 2026







