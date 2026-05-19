# ============================================================
# 01_download_data.R
# Download GSE242330 (scRNA-seq) and GSE242324 (scATAC-seq)
# from Moreaux et al., Blood 2024
# https://ashpublications.org/blood/article/144/5/496/515800
# ============================================================

library(GEOquery)

# ---- Accessions ----
scrna_acc <- "GSE242330"  # scRNA-seq
scatac_acc <- "GSE242324" # scATAC-seq

raw_dir <- "data/raw"
dir.create(file.path(raw_dir, scrna_acc),  showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(raw_dir, scatac_acc), showWarnings = FALSE, recursive = TRUE)

# ============================================================
# STEP 1 — Get metadata (small download, fast)
# ============================================================

message("\n--- Fetching metadata for ", scrna_acc, " ---")
gse_rna <- getGEO(scrna_acc, GSEMatrix = TRUE, getGPL = FALSE,
                  destdir = file.path(raw_dir, scrna_acc))

message("\n--- Fetching metadata for ", scatac_acc, " ---")
gse_atac <- getGEO(scatac_acc, GSEMatrix = TRUE, getGPL = FALSE,
                   destdir = file.path(raw_dir, scatac_acc))

# Inspect sample metadata
message("\n=== scRNA-seq sample summary ===")
print(pData(gse_rna[[1]])[, c("title", "geo_accession", "source_name_ch1")])

message("\n=== scATAC-seq sample summary ===")
print(pData(gse_atac[[1]])[, c("title", "geo_accession", "source_name_ch1")])

# Save metadata tables for the case study
write.csv(pData(gse_rna[[1]]),  "data/processed/metadata_scRNA.csv",  row.names = FALSE)
write.csv(pData(gse_atac[[1]]), "data/processed/metadata_scATAC.csv", row.names = FALSE)

# ============================================================
# STEP 2 — Download supplementary processed files
# ============================================================
# GEO usually provides processed count matrices in the
# "supplementary files" section. These are what we want
# (raw FASTQ would require 100s of GB and Cell Ranger).
# ============================================================

message("\n--- Downloading scRNA-seq supplementary files ---")
getGEOSuppFiles(scrna_acc,
                baseDir = raw_dir,
                fetch_files = TRUE,
                makeDirectory = FALSE)

message("\n--- Downloading scATAC-seq supplementary files ---")
getGEOSuppFiles(scatac_acc,
                baseDir = raw_dir,
                fetch_files = TRUE,
                makeDirectory = FALSE)

# ============================================================
# STEP 3 — Inventory what we downloaded
# ============================================================

message("\n=== Files in scRNA folder ===")
print(list.files(file.path(raw_dir, scrna_acc), recursive = TRUE))

message("\n=== Files in scATAC folder ===")
print(list.files(file.path(raw_dir, scatac_acc), recursive = TRUE))

message("\nDownload complete.")
message("Next: run 02_inspect_data.R to understand the file formats.")
