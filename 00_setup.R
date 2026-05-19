# ============================================================
# 00_setup.R
# Plasma Cell Multi-omics Case Study
# Author: Fabiola
# Purpose: Install required packages and set up project structure
# ============================================================

# ---- Project structure ----
dirs <- c(
  "data/raw",
  "data/processed",
  "results/figures",
  "results/tables",
  "scripts"
)

for (d in dirs) {
  if (!dir.exists(d)) {
    dir.create(d, recursive = TRUE)
    message("Created: ", d)
  }
}

# ---- CRAN packages ----
cran_pkgs <- c(
  "Seurat",          # scRNA-seq
  "Signac",          # scATAC-seq
  "harmony",         # batch integration
  "dplyr",
  "ggplot2",
  "patchwork",
  "pheatmap",
  "RColorBrewer",
  "GEOquery",        # GEO data download
  "R.utils",         # gunzip
  "Matrix"
)

missing_cran <- cran_pkgs[!cran_pkgs %in% installed.packages()[, "Package"]]
if (length(missing_cran) > 0) {
  install.packages(missing_cran)
}

# ---- Bioconductor packages ----
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

bioc_pkgs <- c(
  "GEOquery",
  "GenomicRanges",
  "EnsDb.Hsapiens.v86",
  "BSgenome.Hsapiens.UCSC.hg38",
  "JASPAR2020",
  "TFBSTools",
  "motifmatchr",
  "chromVAR"
)

missing_bioc <- bioc_pkgs[!bioc_pkgs %in% installed.packages()[, "Package"]]
if (length(missing_bioc) > 0) {
  BiocManager::install(missing_bioc, update = FALSE, ask = FALSE)
}

# ---- Quick sanity check ----
library(Seurat)
library(Signac)
library(GEOquery)

message("\n=== Setup complete ===")
message("Seurat:  ", as.character(packageVersion("Seurat")))
message("Signac:  ", as.character(packageVersion("Signac")))
message("R:       ", R.version.string)
