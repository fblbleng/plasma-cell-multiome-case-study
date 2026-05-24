# ============================================================
# 05_regression_comparison.R
# Three-way comparison of cell cycle regression strategies.
#
# The paper (Moreaux et al., Blood 2024) used:
#   - Full regression (S.Score + G2M.Score) for the all-cells
#     trajectory (Analysis 1)
#   - No regression for clustering within prePB + PB (Analysis 2)
#
# We reproduce both and add CC.Difference regression as a
# methodological exploration. This script generates the same
# UMAP + clustering for all three approaches, side by side.
#
# Methods rationale:
#   - LogNormalize (not SCTransform) due to 10x library size
#     variance across stages
#   - Ig and TR genes excluded from variable features in all
#     three approaches
#   - Harmony integration on replicate in all three
#   - Same number of PCs (30) and clustering resolution (0.5)
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(scales)
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
      panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold", color = "#2a2622")
    )
}

# ------------------------------------------------------------
# 1. Load the cell-cycle-scored object
# ------------------------------------------------------------
seu_base <- readRDS(file.path(out_dir, "seurat_with_cellcycle.rds"))
seu_base$stage <- factor(seu_base$stage, levels = c("MBC", "prePB", "PB", "PC"))
seu_base$replicate <- as.factor(seu_base$replicate)
seu_base$Phase <- factor(seu_base$Phase, levels = c("G1", "S", "G2M"))
seu_base$CC.Difference <- seu_base$S.Score - seu_base$G2M.Score
message("Base object: ", ncol(seu_base), " cells, ", nrow(seu_base), " genes")

# ------------------------------------------------------------
# 2. Exclude Ig and TR genes
# ------------------------------------------------------------
ig_genes <- grep("^IG[HKL][VDJC]", rownames(seu_base), value = TRUE)
ig_genes <- unique(c(ig_genes, grep("^JCHAIN$|^JCH$", rownames(seu_base), value = TRUE)))
tr_genes <- grep("^TR[ABGD][VDJC]", rownames(seu_base), value = TRUE)
genes_exclude <- unique(c(ig_genes, tr_genes))
message("Excluded from variable features: ", length(genes_exclude), " genes")

# ------------------------------------------------------------
# 3. Pipeline function that takes regression vars
# ------------------------------------------------------------
# Wraps LogNormalize -> variable features -> scale (with given
# regressors) -> PCA -> Harmony -> UMAP -> clustering.
# Returns the integrated Seurat object.
# ------------------------------------------------------------
run_pipeline <- function(seu, regressors, label, n_pcs = 30, res = 0.5) {
  message("\n========================================")
  message("PIPELINE: ", label)
  message("Regressors: ", paste(regressors, collapse = ", "))
  message("========================================")

  DefaultAssay(seu) <- "RNA"

  # LogNormalize
  seu <- NormalizeData(seu, normalization.method = "LogNormalize",
                       scale.factor = 10000, verbose = FALSE)

  # Variable features per replicate
  seu_list <- SplitObject(seu, split.by = "replicate")
  seu_list <- lapply(seu_list, function(obj) {
    obj <- FindVariableFeatures(obj, selection.method = "vst",
                                 nfeatures = 3000, verbose = FALSE)
    VariableFeatures(obj) <- setdiff(VariableFeatures(obj), genes_exclude)
    obj
  })
  vf_combined <- SelectIntegrationFeatures(seu_list, nfeatures = 3000)
  vf_combined <- setdiff(vf_combined, genes_exclude)
  VariableFeatures(seu) <- vf_combined
  message("Variable features: ", length(vf_combined))

  # Scale with regressors
  seu <- ScaleData(seu, features = vf_combined,
                   vars.to.regress = regressors, verbose = FALSE)

  # PCA
  seu <- RunPCA(seu, features = vf_combined, npcs = 50, verbose = FALSE)

  # Harmony
  seu <- RunHarmony(seu, group.by.vars = "replicate",
                    reduction.use = "pca", reduction.save = "harmony",
                    verbose = FALSE)

  # UMAP
  seu <- RunUMAP(seu, dims = 1:n_pcs, reduction = "harmony",
                 reduction.name = "umap", verbose = FALSE)

  # Clustering
  seu <- FindNeighbors(seu, reduction = "harmony", dims = 1:n_pcs, verbose = FALSE)
  seu <- FindClusters(seu, resolution = res, verbose = FALSE)
  seu$seurat_clusters <- Idents(seu)

  message(label, " clusters: ", length(unique(seu$seurat_clusters)))
  return(seu)
}

# ------------------------------------------------------------
# 4. Run all three pipelines
# ------------------------------------------------------------
seu_A <- run_pipeline(seu_base,
                     regressors = c("percent.mt", "S.Score", "G2M.Score"),
                     label = "A_full_regression")

seu_B <- run_pipeline(seu_base,
                     regressors = c("percent.mt", "CC.Difference"),
                     label = "B_CCdifference_regression")

seu_C <- run_pipeline(seu_base,
                     regressors = c("percent.mt"),
                     label = "C_no_regression")

# ------------------------------------------------------------
# 5. Save all three integrated objects
# ------------------------------------------------------------
saveRDS(seu_A, file.path(out_dir, "seurat_A_full_regression.rds"))
saveRDS(seu_B, file.path(out_dir, "seurat_B_CCdiff_regression.rds"))
saveRDS(seu_C, file.path(out_dir, "seurat_C_no_regression.rds"))
message("\nSaved all three integrated objects.")

# ------------------------------------------------------------
# 6. Helper to make a 4-panel UMAP figure for one approach
# ------------------------------------------------------------
make_umap_panel <- function(seu, title_text, subtitle_text) {
  p_rep <- DimPlot(seu, reduction = "umap", group.by = "replicate",
                   cols = rep_colors, pt.size = 0.3) +
    theme_portfolio() + theme(legend.position = "right") +
    labs(title = "by replicate", x = NULL, y = NULL)

  p_stage <- DimPlot(seu, reduction = "umap", group.by = "stage",
                     cols = stage_colors, pt.size = 0.3) +
    theme_portfolio() + theme(legend.position = "right") +
    labs(title = "by stage", x = NULL, y = NULL)

  p_cluster <- DimPlot(seu, reduction = "umap", group.by = "seurat_clusters",
                       label = TRUE, label.size = 3.5, pt.size = 0.3, repel = TRUE) +
    theme_portfolio() + theme(legend.position = "none") +
    labs(title = "by cluster (res 0.5)", x = NULL, y = NULL)

  p_phase <- DimPlot(seu, reduction = "umap", group.by = "Phase",
                     cols = phase_colors, pt.size = 0.3) +
    theme_portfolio() + theme(legend.position = "right") +
    labs(title = "by cell cycle phase", x = NULL, y = NULL)

  (p_rep | p_stage | p_cluster | p_phase) +
    plot_annotation(
      title = title_text,
      subtitle = subtitle_text,
      theme = theme(
        plot.title = element_text(face = "bold", size = 14, color = "#2a2622"),
        plot.subtitle = element_text(size = 10, color = "#5c544b")
      )
    )
}

# Individual full-detail figures per approach
fig_A <- make_umap_panel(seu_A,
                         "A. Full cell cycle regression (S.Score + G2M.Score)",
                         "Matches paper Analysis 1 - reveals trajectory by removing all cycle signal")
ggsave(file.path(fig_dir, "07A_umap_full_regression.pdf"), fig_A, width = 16, height = 4.5)
ggsave(file.path(fig_dir, "07A_umap_full_regression.png"), fig_A, width = 16, height = 4.5, dpi = 200)

fig_B <- make_umap_panel(seu_B,
                         "B. CC.Difference regression (S.Score - G2M.Score)",
                         "Removes phase identity (S vs G2M) while preserving cycling-vs-quiescent distinction")
ggsave(file.path(fig_dir, "07B_umap_ccdiff_regression.pdf"), fig_B, width = 16, height = 4.5)
ggsave(file.path(fig_dir, "07B_umap_ccdiff_regression.png"), fig_B, width = 16, height = 4.5, dpi = 200)

fig_C <- make_umap_panel(seu_C,
                         "C. No cell cycle regression (only percent.mt)",
                         "Keeps all cell cycle biology - shows how proliferation shapes the trajectory")
ggsave(file.path(fig_dir, "07C_umap_no_regression.pdf"), fig_C, width = 16, height = 4.5)
ggsave(file.path(fig_dir, "07C_umap_no_regression.png"), fig_C, width = 16, height = 4.5, dpi = 200)

# ------------------------------------------------------------
# 7. The METHODS COMPARISON figure: 3x2 grid (3 approaches x stage/phase)
# ------------------------------------------------------------
mk_small <- function(seu, group, cols, title_text, legend = "none") {
  DimPlot(seu, reduction = "umap", group.by = group,
          cols = cols, pt.size = 0.25) +
    theme_portfolio() +
    theme(legend.position = legend) +
    labs(title = title_text, x = NULL, y = NULL)
}

# Row 1: stage coloring
r1_A <- mk_small(seu_A, "stage", stage_colors,
                 "A. Full regression",   legend = "none")
r1_B <- mk_small(seu_B, "stage", stage_colors,
                 "B. CC.Difference",     legend = "none")
r1_C <- mk_small(seu_C, "stage", stage_colors,
                 "C. No regression",     legend = "right")

# Row 2: phase coloring
r2_A <- mk_small(seu_A, "Phase", phase_colors, "", legend = "none")
r2_B <- mk_small(seu_B, "Phase", phase_colors, "", legend = "none")
r2_C <- mk_small(seu_C, "Phase", phase_colors, "", legend = "right")

methods_comp <- (r1_A | r1_B | r1_C) / (r2_A | r2_B | r2_C) +
  plot_annotation(
    title = "Methods comparison: three cell cycle regression strategies",
    subtitle = "Row 1: colored by stage | Row 2: colored by cell cycle phase | Same cells, same Harmony, same UMAP params",
    theme = theme(
      plot.title = element_text(face = "bold", size = 15, color = "#2a2622"),
      plot.subtitle = element_text(size = 11, color = "#5c544b")
    )
  )

ggsave(file.path(fig_dir, "08_methods_comparison.pdf"), methods_comp,
       width = 14, height = 8)
ggsave(file.path(fig_dir, "08_methods_comparison.png"), methods_comp,
       width = 14, height = 8, dpi = 200)
print(methods_comp)

# ------------------------------------------------------------
# 8. Cross-tabulations + cluster counts summary
# ------------------------------------------------------------
summarize_clustering <- function(seu, label) {
  xtab <- table(Cluster = seu$seurat_clusters, Stage = seu$stage)
  message("\n--- ", label, " cluster x stage ---")
  print(xtab)
  write.csv(as.data.frame.matrix(xtab),
            file.path(tbl_dir, paste0("04_clusterxstage_", label, ".csv")))
  return(data.frame(approach = label,
                    n_clusters = length(unique(seu$seurat_clusters)),
                    stringsAsFactors = FALSE))
}

cluster_counts <- bind_rows(
  summarize_clustering(seu_A, "A_full_regression"),
  summarize_clustering(seu_B, "B_CCdiff_regression"),
  summarize_clustering(seu_C, "C_no_regression")
)
message("\nCluster count summary:")
print(cluster_counts)
write.csv(cluster_counts, file.path(tbl_dir, "04_cluster_counts_comparison.csv"),
          row.names = FALSE)

# ------------------------------------------------------------
# 9. The "headline" object for downstream
# ------------------------------------------------------------
# Per paper Analysis 1, we use Approach A (full regression)
# as the canonical all-cells analysis for the headline figure
# and downstream differential expression.
# ------------------------------------------------------------
saveRDS(seu_A, file.path(out_dir, "seurat_integrated_clustered.rds"))
message("\nApproach A saved as seurat_integrated_clustered.rds for downstream use")

message("\n========================================")
message("Stage 2.3 (Three-way regression comparison) complete!")
message("\nFigures:")
message("  - 07A_umap_full_regression       (paper Analysis 1)")
message("  - 07B_umap_ccdiff_regression     (our exploration)")
message("  - 07C_umap_no_regression         (paper Analysis 2 baseline)")
message("  - 08_methods_comparison          <-- METHODS SHOWCASE")
message("\nNext:")
message("  - 06_subset_prePB_PB.R           (paper Analysis 2)")
message("  - 07_marker_validation.R         (annotate headline UMAP)")
message("========================================")
