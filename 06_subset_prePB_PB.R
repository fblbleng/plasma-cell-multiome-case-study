# ============================================================
# 06_subset_prePB_PB.R   (v2 - with Seurat v5 JoinLayers fix)
# Sub-analysis of proliferating cells (prePB + PB only).
#
# Matches Moreaux et al. Analysis 2:
#   - Subset to prePB + PB cells
#   - Re-normalize and re-cluster WITHOUT cell cycle regression
#   - Clustering resolution 0.2
#   - Goal: discover regulatory programs WITHIN proliferating
#     cells (the paper found AP-1 / BATF heterogeneity here)
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

stage_colors <- c("prePB" = "#b88a3e", "PB" = "#b15835")
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
# 1. Load and subset
# ------------------------------------------------------------
seu_base <- readRDS(file.path(out_dir, "seurat_with_cellcycle.rds"))
seu_base$stage <- factor(seu_base$stage, levels = c("MBC", "prePB", "PB", "PC"))
seu_base$replicate <- as.factor(seu_base$replicate)
seu_base$Phase <- factor(seu_base$Phase, levels = c("G1", "S", "G2M"))

seu <- subset(seu_base, subset = stage %in% c("prePB", "PB"))
seu$stage <- droplevels(seu$stage)
message("Subset: ", ncol(seu), " cells")
print(table(seu$stage, seu$replicate))

# ------------------------------------------------------------
# 2. Seurat v5: join layers BEFORE downstream analysis
# ------------------------------------------------------------
# The base object has split counts layers (counts.NPCD_rep1,
# counts.NPCD_rep2). FindAllMarkers etc. require these joined
# into a single 'counts' layer. We do this before normalization
# so the entire pipeline is consistent.
# ------------------------------------------------------------
DefaultAssay(seu) <- "RNA"
seu <- JoinLayers(seu, assay = "RNA")
message("\nLayers after JoinLayers():")
print(Layers(seu, assay = "RNA"))

# ------------------------------------------------------------
# 3. Exclude Ig and TR
# ------------------------------------------------------------
ig_genes <- grep("^IG[HKL][VDJC]", rownames(seu), value = TRUE)
ig_genes <- unique(c(ig_genes, grep("^JCHAIN$|^JCH$", rownames(seu), value = TRUE)))
tr_genes <- grep("^TR[ABGD][VDJC]", rownames(seu), value = TRUE)
genes_exclude <- unique(c(ig_genes, tr_genes))

# ------------------------------------------------------------
# 4. LogNormalize
# ------------------------------------------------------------
message("\nLogNormalize on subset...")
seu <- NormalizeData(seu, normalization.method = "LogNormalize",
                     scale.factor = 10000, verbose = FALSE)

# ------------------------------------------------------------
# 5. Variable features (per replicate, then combine)
# ------------------------------------------------------------
# After JoinLayers, SplitObject won't separate the counts again
# (the split is metadata-based, not layer-based). FindVariableFeatures
# on the merged object is fine; we just compute it once.
# ------------------------------------------------------------
message("\nFinding variable features...")
seu <- FindVariableFeatures(seu, selection.method = "vst",
                             nfeatures = 2500, verbose = FALSE)
vf <- VariableFeatures(seu)
vf_clean <- setdiff(vf, genes_exclude)
VariableFeatures(seu) <- vf_clean
message("Variable features: ", length(vf), " -> ", length(vf_clean),
        " after Ig/TR exclusion")

# ------------------------------------------------------------
# 6. Scale WITHOUT cell cycle regression
# ------------------------------------------------------------
message("\nScaling (regressing only percent.mt - NO cell cycle)...")
seu <- ScaleData(seu, features = vf_clean,
                 vars.to.regress = "percent.mt", verbose = FALSE)

# ------------------------------------------------------------
# 7. PCA, Harmony, UMAP, clustering
# ------------------------------------------------------------
message("\nPCA + Harmony + UMAP + clustering...")
seu <- RunPCA(seu, features = vf_clean, npcs = 50, verbose = FALSE)

elbow <- ElbowPlot(seu, ndims = 50) +
  theme_portfolio() +
  labs(title = "PCA elbow (prePB + PB subset)")
ggsave(file.path(fig_dir, "09_subset_pca_elbow.pdf"), elbow, width = 6, height = 4)

n_pcs <- 25

seu <- RunHarmony(seu, group.by.vars = "replicate",
                  reduction.use = "pca", reduction.save = "harmony",
                  verbose = FALSE)
seu <- RunUMAP(seu, dims = 1:n_pcs, reduction = "harmony",
               reduction.name = "umap", verbose = FALSE)
seu <- FindNeighbors(seu, reduction = "harmony", dims = 1:n_pcs, verbose = FALSE)
seu <- FindClusters(seu, resolution = 0.2, verbose = FALSE)
seu$subset_clusters <- Idents(seu)

message("\nSubset clusters at res 0.2:")
print(table(seu$subset_clusters))

# ------------------------------------------------------------
# 8. Visualization - 4 panels
# ------------------------------------------------------------
p_rep <- DimPlot(seu, reduction = "umap", group.by = "replicate",
                 cols = rep_colors, pt.size = 0.4) +
  theme_portfolio() + theme(legend.position = "right") +
  labs(title = "Subset: by replicate")

p_stage <- DimPlot(seu, reduction = "umap", group.by = "stage",
                   cols = stage_colors, pt.size = 0.4) +
  theme_portfolio() + theme(legend.position = "right") +
  labs(title = "Subset: by stage")

p_cluster <- DimPlot(seu, reduction = "umap", group.by = "subset_clusters",
                     label = TRUE, label.size = 4, pt.size = 0.4, repel = TRUE) +
  theme_portfolio() + theme(legend.position = "none") +
  labs(title = "Subset: by cluster (res 0.2)")

p_phase <- DimPlot(seu, reduction = "umap", group.by = "Phase",
                   cols = phase_colors, pt.size = 0.4) +
  theme_portfolio() + theme(legend.position = "right") +
  labs(title = "Subset: by cell cycle phase")

subset_fig <- (p_rep | p_stage) / (p_cluster | p_phase) +
  plot_annotation(
    title = "Sub-analysis of proliferating cells (prePB + PB)",
    subtitle = "No cell cycle regression | res 0.2 | matches paper Analysis 2",
    theme = theme(
      plot.title = element_text(face = "bold", size = 15, color = "#2a2622"),
      plot.subtitle = element_text(size = 11, color = "#5c544b")
    )
  )

ggsave(file.path(fig_dir, "10_subset_umap.pdf"), subset_fig, width = 12, height = 10)
ggsave(file.path(fig_dir, "10_subset_umap.png"), subset_fig, width = 12, height = 10, dpi = 200)
print(subset_fig)

# ------------------------------------------------------------
# 9. Find cluster markers in the subset
# ------------------------------------------------------------
# Layers are already joined (step 2), so FindAllMarkers works.
# ------------------------------------------------------------
message("\nFinding markers per subset cluster (this takes ~2-5 min)...")
Idents(seu) <- "subset_clusters"
markers_subset <- FindAllMarkers(seu,
                                  only.pos = TRUE,
                                  min.pct = 0.25,
                                  logfc.threshold = 0.25,
                                  verbose = FALSE)

if (nrow(markers_subset) == 0) {
  warning("No markers identified. This is rare - check that JoinLayers ran.")
} else {
  top_markers <- markers_subset %>%
    group_by(cluster) %>%
    slice_max(order_by = avg_log2FC, n = 20) %>%
    ungroup()

  write.csv(markers_subset, file.path(tbl_dir, "05_subset_all_markers.csv"),
            row.names = FALSE)
  write.csv(top_markers, file.path(tbl_dir, "05_subset_top20_per_cluster.csv"),
            row.names = FALSE)

  message("\nTop 5 markers per subset cluster:")
  print(top_markers %>%
          group_by(cluster) %>%
          slice_head(n = 5) %>%
          select(cluster, gene, avg_log2FC, p_val_adj))
}

# ------------------------------------------------------------
# 10. AP-1 / BATF / FOS family visualization
# ------------------------------------------------------------
ap1_genes <- c("FOS", "FOSB", "FOSL1", "FOSL2",
               "JUN", "JUNB", "JUND",
               "BATF", "BATF2", "BATF3",
               "ATF3", "ATF4")
ap1_genes <- ap1_genes[ap1_genes %in% rownames(seu)]
message("\nAP-1 family genes available: ", length(ap1_genes))
message("  ", paste(ap1_genes, collapse = ", "))

if (length(ap1_genes) > 0) {
  p_ap1_dot <- DotPlot(seu, features = ap1_genes, group.by = "subset_clusters",
                       cols = c("#f3ede3", "#b15835"), dot.scale = 8) +
    theme_portfolio() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "AP-1 family expression across subset clusters",
         subtitle = "FOS / JUN / BATF / ATF families",
         x = NULL, y = NULL)

  ggsave(file.path(fig_dir, "11_subset_AP1_dotplot.pdf"),
         p_ap1_dot, width = 9, height = 5)
  ggsave(file.path(fig_dir, "11_subset_AP1_dotplot.png"),
         p_ap1_dot, width = 9, height = 5, dpi = 200)

  top_ap1 <- intersect(c("FOS", "JUN", "JUND", "BATF", "BATF3", "ATF3"),
                       rownames(seu))
  if (length(top_ap1) > 0) {
    p_feat <- FeaturePlot(seu, features = top_ap1, reduction = "umap",
                          cols = c("#f3ede3", "#b15835"),
                          ncol = 3, pt.size = 0.3) &
      theme_portfolio() & NoAxes()
    ggsave(file.path(fig_dir, "12_subset_AP1_featureplot.pdf"),
           p_feat, width = 12, height = ceiling(length(top_ap1) / 3) * 4)
    ggsave(file.path(fig_dir, "12_subset_AP1_featureplot.png"),
           p_feat, width = 12, height = ceiling(length(top_ap1) / 3) * 4, dpi = 200)
  }
}

# ------------------------------------------------------------
# 11. Cluster x stage cross-tab
# ------------------------------------------------------------
xtab_sub <- table(Cluster = seu$subset_clusters, Stage = seu$stage)
message("\nSubset cluster x stage:")
print(xtab_sub)
write.csv(as.data.frame.matrix(xtab_sub),
          file.path(tbl_dir, "05_subset_clusterxstage.csv"))

# ------------------------------------------------------------
# 12. Save
# ------------------------------------------------------------
saveRDS(seu, file.path(out_dir, "seurat_prePB_PB_subset.rds"))

message("\n========================================")
message("Stage 2.4 (prePB + PB subset) complete!")
message("Subset cells: ", ncol(seu))
message("Subset clusters at res 0.2: ", length(unique(seu$subset_clusters)))
message("========================================")
