# ============================================================
# 05b_retune_umap_17dims.R
# Re-run Harmony, UMAP, and clustering on the THREE all-cells
# objects (A, B, C) using 17 PCs to match the paper.
#
# Original script 05 used 30 PCs - this caused the UMAP to
# spread too much. Moreaux et al. used 17 dimensions per their
# methods. We keep all other parameters identical so the
# regression comparison remains valid.
#
# This is FAST because we skip the expensive parts:
#   - LogNormalize: cached
#   - Variable features: cached
#   - ScaleData: cached
#   - PCA: cached (we just use fewer dims)
# We only re-run: Harmony, UMAP, FindNeighbors, FindClusters,
# and re-make the figures + cross-tabs.
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
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
rep_colors   <- c("1" = "#b15835", "2" = "#7d8c6e")
phase_colors <- c("G1" = "#7d8c6e", "S" = "#b88a3e", "G2M" = "#b15835")

theme_portfolio <- function() {
  theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 12, color = "#2a2622"),
      plot.subtitle = element_text(size = 10, color = "#5c544b"),
      axis.text = element_text(color = "#2a2622", size = 9),
      axis.title = element_text(color = "#5c544b", size = 10),
      panel.grid.minor = element_blank()
    )
}

# ------------------------------------------------------------
# Tunable parameters
# ------------------------------------------------------------
n_pcs   <- 10        # matches paper
res     <- 0.2
n_neigh <- 30        # UMAP default; lower if still too spread
m_dist  <- 0.1       # UMAP default; lower for tighter clusters

message("Re-running with parameters:")
message("  n_pcs       = ", n_pcs, " (matches paper)")
message("  resolution  = ", res)
message("  n.neighbors = ", n_neigh)
message("  min.dist    = ", m_dist)

# ------------------------------------------------------------
# Helper: retune a single object
# ------------------------------------------------------------
retune <- function(filename, label, descriptor) {
  message("\n========================================")
  message("Re-tuning: ", label, " (", filename, ")")
  message("========================================")

  fp <- file.path(out_dir, filename)
  if (!file.exists(fp)) {
    stop("File not found: ", fp)
  }
  seu <- readRDS(fp)
  message("Loaded ", ncol(seu), " cells")

  # Set factor levels
  seu$stage <- factor(seu$stage, levels = c("MBC", "prePB", "PB", "PC"))
  seu$replicate <- as.factor(seu$replicate)
  seu$Phase <- factor(seu$Phase, levels = c("G1", "S", "G2M"))

  # Re-run Harmony on first 17 PCA dims
  message("Re-running Harmony on ", n_pcs, " PCs...")
  seu <- RunHarmony(seu, group.by.vars = "replicate",
                    reduction.use = "pca", dims.use = 1:n_pcs,
                    reduction.save = "harmony", verbose = FALSE)

  # New UMAP on harmony space
  message("Computing UMAP...")
  seu <- RunUMAP(seu,
                 dims = 1:n_pcs,
                 reduction = "harmony",
                 reduction.name = "umap",
                 n.neighbors = n_neigh,
                 min.dist = m_dist,
                 verbose = FALSE)

  # Re-cluster
  message("Re-clustering at res ", res, "...")
  seu <- FindNeighbors(seu, reduction = "harmony", dims = 1:n_pcs, verbose = FALSE)
  seu <- FindClusters(seu, resolution = res, verbose = FALSE)
  seu$seurat_clusters <- Idents(seu)

  message("Clusters: ", length(unique(seu$seurat_clusters)))
  print(table(seu$seurat_clusters))

  # 4-panel figure
  p_rep <- DimPlot(seu, reduction = "umap", group.by = "replicate",
                   cols = rep_colors, pt.size = 0.3) +
    theme_portfolio() + theme(legend.position = "right") +
    labs(title = "by replicate", x = NULL, y = NULL)

  p_stage <- DimPlot(seu, reduction = "umap", group.by = "stage",
                     cols = stage_colors, pt.size = 0.3) +
    theme_portfolio() + theme(legend.position = "right") +
    labs(title = "by stage", x = NULL, y = NULL)

  p_cluster <- DimPlot(seu, reduction = "umap", group.by = "seurat_clusters",
                       label = TRUE, label.size = 3.5, pt.size = 0.3,
                       repel = TRUE) +
    theme_portfolio() + theme(legend.position = "none") +
    labs(title = paste0("by cluster (res ", res, ")"), x = NULL, y = NULL)

  p_phase <- DimPlot(seu, reduction = "umap", group.by = "Phase",
                     cols = phase_colors, pt.size = 0.3) +
    theme_portfolio() + theme(legend.position = "right") +
    labs(title = "by cell cycle phase", x = NULL, y = NULL)

  fig <- (p_rep | p_stage | p_cluster | p_phase) +
    plot_annotation(
      title = paste0(label, " (", n_pcs, " PCs)"),
      subtitle = descriptor,
      theme = theme(
        plot.title = element_text(face = "bold", size = 14, color = "#2a2622"),
        plot.subtitle = element_text(size = 10, color = "#5c544b")
      )
    )
  return(list(seu = seu, fig = fig))
}

# ------------------------------------------------------------
# Run the retune on all three objects
# ------------------------------------------------------------
A <- retune("seurat_A_full_regression.rds",
            "A. Full cell cycle regression (S.Score + G2M.Score)",
            "Matches paper Analysis 1 - removes cell cycle signal")

B <- retune("seurat_B_CCdiff_regression.rds",
            "B. CC.Difference regression (S.Score - G2M.Score)",
            "Removes phase identity while preserving cycling status")

C <- retune("seurat_C_no_regression.rds",
            "C. No cell cycle regression (only percent.mt)",
            "Keeps full cell cycle biology")

# ------------------------------------------------------------
# Save individual figures
# ------------------------------------------------------------
ggsave(file.path(fig_dir, "07A_umap_full_regression.pdf"),
       A$fig, width = 16, height = 4.5)
ggsave(file.path(fig_dir, "07A_umap_full_regression.png"),
       A$fig, width = 16, height = 4.5, dpi = 200)

ggsave(file.path(fig_dir, "07B_umap_ccdiff_regression.pdf"),
       B$fig, width = 16, height = 4.5)
ggsave(file.path(fig_dir, "07B_umap_ccdiff_regression.png"),
       B$fig, width = 16, height = 4.5, dpi = 200)

ggsave(file.path(fig_dir, "07C_umap_no_regression.pdf"),
       C$fig, width = 16, height = 4.5)
ggsave(file.path(fig_dir, "07C_umap_no_regression.png"),
       C$fig, width = 16, height = 4.5, dpi = 200)

# ------------------------------------------------------------
# Methods comparison panel
# ------------------------------------------------------------
mk_small <- function(seu, group, cols, title_text, legend = "none") {
  DimPlot(seu, reduction = "umap", group.by = group,
          cols = cols, pt.size = 0.25) +
    theme_portfolio() +
    theme(legend.position = legend) +
    labs(title = title_text, x = NULL, y = NULL)
}

r1_A <- mk_small(A$seu, "stage", stage_colors, "A. Full regression",   "none")
r1_B <- mk_small(B$seu, "stage", stage_colors, "B. CC.Difference",     "none")
r1_C <- mk_small(C$seu, "stage", stage_colors, "C. No regression",     "right")

r2_A <- mk_small(A$seu, "Phase", phase_colors, "", "none")
r2_B <- mk_small(B$seu, "Phase", phase_colors, "", "none")
r2_C <- mk_small(C$seu, "Phase", phase_colors, "", "right")

methods_comp <- (r1_A | r1_B | r1_C) / (r2_A | r2_B | r2_C) +
  plot_annotation(
    title = "Methods comparison: three cell cycle regression strategies",
    subtitle = paste0("All ", n_pcs, " PCs (matches paper) | Row 1: by stage | Row 2: by phase"),
    theme = theme(
      plot.title = element_text(face = "bold", size = 15, color = "#2a2622"),
      plot.subtitle = element_text(size = 11, color = "#5c544b")
    )
  )

ggsave(file.path(fig_dir, "08_methods_comparison.pdf"),
       methods_comp, width = 14, height = 8)
ggsave(file.path(fig_dir, "08_methods_comparison.png"),
       methods_comp, width = 14, height = 8, dpi = 200)
print(methods_comp)

# ------------------------------------------------------------
# Cluster x stage summaries
# ------------------------------------------------------------
summarize_clustering <- function(seu, label) {
  xtab <- table(Cluster = seu$seurat_clusters, Stage = seu$stage)
  message("\n--- ", label, " ---")
  print(xtab)
  write.csv(as.data.frame.matrix(xtab),
            file.path(tbl_dir, paste0("04_clusterxstage_", label, ".csv")))
  data.frame(approach = label,
             n_clusters = length(unique(seu$seurat_clusters)),
             stringsAsFactors = FALSE)
}

cluster_counts <- bind_rows(
  summarize_clustering(A$seu, "A_full_regression"),
  summarize_clustering(B$seu, "B_CCdiff_regression"),
  summarize_clustering(C$seu, "C_no_regression")
)
message("\nCluster count summary:")
print(cluster_counts)
write.csv(cluster_counts, file.path(tbl_dir, "04_cluster_counts_comparison.csv"),
          row.names = FALSE)

# ------------------------------------------------------------
# Save all three updated objects
# ------------------------------------------------------------
saveRDS(A$seu, file.path(out_dir, "seurat_A_full_regression.rds"))
saveRDS(B$seu, file.path(out_dir, "seurat_B_CCdiff_regression.rds"))
saveRDS(C$seu, file.path(out_dir, "seurat_C_no_regression.rds"))

# A is the headline (paper Analysis 1)
saveRDS(A$seu, file.path(out_dir, "seurat_integrated_clustered.rds"))

message("\n========================================")
message("Retune complete (", n_pcs, " PCs, matches paper)")
message("All three objects updated.")
message("Approach A saved as seurat_integrated_clustered.rds for downstream use.")
message("\nNext: re-run 07_marker_validation.R to refresh marker analysis")
message("on the retuned headline object.")
message("========================================")
