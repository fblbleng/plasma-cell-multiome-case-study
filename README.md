# Plasma Cell Differentiation — Multi-omics Case Study

A reproducible re-analysis and extension of the plasma cell differentiation dataset from:

> **Moreaux et al., *Blood* (2024)**
> *Integrative single-cell chromatin and transcriptome analysis of human plasma cell differentiation*
> [DOI: 10.1182/blood.2023023253](https://ashpublications.org/blood/article/144/5/496/515800)

GEO accessions:
- **GSE242330** — scRNA-seq
- **GSE242324** — scATAC-seq

## Project goal

This case study demonstrates an end-to-end single-cell multi-omics analysis pipeline,
applied to a publicly available dataset of human plasma cell differentiation. The
workflow reproduces key findings from the original paper and **extends** them with
an additional layer of gene regulatory network inference using SCENIC+.

## Workflow

1. **Data acquisition** — GEO download and inspection
2. **scRNA-seq** — QC, normalization, integration, clustering, cell-type annotation
3. **scATAC-seq** — QC, peak calling, dimensionality reduction, motif enrichment
4. **Multi-modal integration** — label transfer and WNN
5. **Differential analysis** — stage-specific DEGs and DARs
6. **Gene regulatory networks (SCENIC+)** — *the original extension this case study adds*
7. **Visualization & reporting** — publication-quality figures

## Reproducibility

```r
# 1. Install dependencies and create project folders
source("00_setup.R")

# 2. Download data from GEO (~ a few GB)
source("01_download_data.R")

# 3. Inspect what was downloaded and verify file formats
source("02_inspect_data.R")
```

## About this case study

This work is part of my computational biology portfolio. The same pipeline is
applied in my research on long-lived plasma cell survival in bone marrow vs.
spleen niches (proprietary data; available under collaboration).

— Fabiola
