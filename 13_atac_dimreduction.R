# ============================================================
# 13_atac_dimreduction.R
#
# Phase 3.2: Dimensionality reduction and clustering of scATAC-seq.
#
# Steps:
#   1. TF-IDF normalization (ATAC-specific, accounts for variable
#      counts per cell)
#   2. Find variable peaks (top features)
#   3. LSI dimensionality reduction (SVD-based, ATAC equivalent of PCA)
#   4. UMAP embedding
#   5. Louvain clustering
#   6. Stage validation: does chromatin alone recover the 4 stages?
#
# Inputs:
#   - data/processed/atac_signac_merged.rds
#
# Outputs:
#   - data/processed/atac_signac_clustered.rds
#   - results/figures/31_atac_lsi_depthcor.{pdf,png}
#   - results/figures/32_atac_umap_by_stage.{pdf,png}
#   - results/figures/33_atac_umap_by_cluster.{pdf,png}
#   - results/figures/34_atac_clusterxstage_heatmap.{pdf,png}
#   - results/tables/12_atac_cluster_stage_composition.csv
# ============================================================

suppressPackageStartupMessages({
  library(Signac)
  library(Seurat)
  library(RSpectra)
  library(Matrix)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

set.seed(42)
out_dir <- "data/processed"
fig_dir <- "results/figures"
tbl_dir <- "results/tables"

stage_colors <- c("MBC"   = "#7d8c6e",
                  "prePB" = "#b88a3e",
                  "PB"    = "#b15835",
                  "PC"    = "#7c5c6b")

# Load the merged object
seu <- readRDS(file.path(out_dir, "atac_signac_merged.rds"))
DefaultAssay(seu) <- "peaks"
message("Loaded: ", ncol(seu), " cells, ", nrow(seu), " peaks")

# TF-IDF
message("\n[1/5] TF-IDF normalization...")
seu <- RunTFIDF(seu)

# Top features
message("\n[2/5] Selecting top features...")
seu <- FindTopFeatures(seu, min.cutoff = "q5")
message("Variable features: ", length(VariableFeatures(seu)))

# Manual SVD via RSpectra (bypasses broken irlba)
message("\n[3/5] Running SVD via RSpectra (irlba bypass)...")
top_features <- VariableFeatures(seu)
mat <- GetAssayData(seu, layer = "data")[top_features, ]
mat <- as(mat, "CsparseMatrix")

svd_res <- RSpectra::svds(t(mat), k = 50)

embeddings <- svd_res$u %*% diag(svd_res$d)
rownames(embeddings) <- colnames(seu)
colnames(embeddings) <- paste0("LSI_", seq_len(50))

loadings <- svd_res$v
rownames(loadings) <- top_features
colnames(loadings) <- paste0("LSI_", seq_len(50))

seu[["lsi"]] <- CreateDimReducObject(
  embeddings = embeddings,
  loadings   = loadings,
  key        = "LSI_",
  assay      = "peaks"
)

message("LSI complete via RSpectra")

# Depth correlation diagnostic
depth_cors <- abs(cor(seu@reductions$lsi@cell.embeddings,
                      seu$nCount_peaks))
depth_problematic <- which(depth_cors > 0.75)
message("Depth-correlated components: ",
        paste(depth_problematic, collapse = ", "))
use_lsi <- setdiff(2:30, depth_problematic)
message("Using LSI components: ", min(use_lsi), "-", max(use_lsi))

# UMAP
message("\n[4/5] UMAP embedding...")
seu <- RunUMAP(seu,
               reduction = "lsi",
               dims = use_lsi,
               verbose = FALSE)

# Clustering
message("\n[5/5] Louvain clustering...")
seu <- FindNeighbors(seu,
                     reduction = "lsi",
                     dims = use_lsi,
                     verbose = FALSE)

for (res in c(0.3, 0.5, 0.8)) {
  seu <- FindClusters(seu,
                      resolution = res,
                      algorithm = 3,
                      verbose = FALSE)
  message("Resolution ", res, ": ",
          length(unique(Idents(seu))), " clusters")
}

seu$atac_clusters <- seu$peaks_snn_res.0.5
Idents(seu) <- "atac_clusters"

# Quick figures
message("\nGenerating UMAP plots...")

p_umap_stage <- DimPlot(seu, group.by = "stage",
                        cols = stage_colors,
                        pt.size = 0.5,
                        raster = FALSE) +
  ggtitle(paste0("scATAC-seq UMAP, by stage (", ncol(seu), " cells)"))

p_umap_cluster <- DimPlot(seu, group.by = "atac_clusters",
                          label = TRUE, label.size = 5,
                          pt.size = 0.5,
                          raster = FALSE) +
  ggtitle("scATAC-seq UMAP, unbiased clusters")

ggsave(file.path(fig_dir, "32_atac_umap_by_stage.png"),
       p_umap_stage, width = 9, height = 7, dpi = 200)
ggsave(file.path(fig_dir, "33_atac_umap_by_cluster.png"),
       p_umap_cluster, width = 9, height = 7, dpi = 200)

# Cluster x stage composition
cluster_stage <- table(seu$atac_clusters, seu$stage)
print(cluster_stage)
cluster_stage_pct <- prop.table(cluster_stage, margin = 1) * 100
print(round(cluster_stage_pct, 1))

write.csv(as.data.frame.matrix(cluster_stage),
          file.path(tbl_dir, "12_atac_cluster_stage_composition.csv"),
          row.names = TRUE)

# Save
saveRDS(seu, file.path(out_dir, "atac_signac_clustered.rds"))

message("\n========================================")
message("Phase 3.2 complete (via RSpectra bypass)!")
message("Clusters: ", length(unique(seu$atac_clusters)))
message("Object saved: atac_signac_clustered.rds")
message("========================================")