# ============================================================
# 07_marker_validation.R   (v2)
# Marker gene validation + cluster annotation for the headline
# UMAP. Uses Approach B (CC.Difference regression) as the
# headline because it preserves the proliferation axis while
# removing within-cycle phase identity.
#
# This script:
#   - Loads the headline integrated object
#   - Annotates each cluster with its dominant stage
#   - Validates canonical B-cell -> PC markers across stages
#     and clusters (dot plot + feature plot)
#   - Finds top markers per cluster (full table, not truncated)
#   - Generates a top-markers-per-cluster heatmap
#   - Specifically compares the two PC clusters to find what
#     distinguishes them biologically
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(pheatmap)
  library(RColorBrewer)
})

set.seed(42)
out_dir <- "data/processed"
fig_dir <- "results/figures"
tbl_dir <- "results/tables"

stage_colors <- c("MBC"   = "#7d8c6e",
                  "prePB" = "#b88a3e",
                  "PB"    = "#b15835",
                  "PC"    = "#7c5c6b")

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
# 1. Load the headline (Approach B) integrated object
# ------------------------------------------------------------
seu <- readRDS(file.path(out_dir, "seurat_integrated_clustered.rds"))
seu$stage <- factor(seu$stage, levels = c("MBC", "prePB", "PB", "PC"))
Idents(seu) <- "seurat_clusters"
message("Loaded headline object: ", ncol(seu), " cells, ",
        length(unique(seu$seurat_clusters)), " clusters")
message("(Headline = Approach B: CC.Difference regression)")

# ------------------------------------------------------------
# 2. Cluster annotation by dominant stage
# ------------------------------------------------------------
xtab     <- table(Cluster = seu$seurat_clusters, Stage = seu$stage)
xtab_pct <- prop.table(xtab, margin = 1) * 100

cluster_annot <- data.frame(
  cluster        = rownames(xtab),
  dominant_stage = colnames(xtab)[apply(xtab, 1, which.max)],
  pct            = round(apply(xtab_pct, 1, max), 0),
  n_cells        = rowSums(xtab),
  stringsAsFactors = FALSE
)
cluster_annot$annotation <- paste0("C", cluster_annot$cluster, "_",
                                    cluster_annot$dominant_stage, "_",
                                    cluster_annot$pct, "%")
message("\nCluster annotations:")
print(cluster_annot)
write.csv(cluster_annot,
          file.path(tbl_dir, "06_cluster_annotations.csv"),
          row.names = FALSE)

# Apply annotation to cells (the Seurat v5-safe way)
new_annot <- cluster_annot$annotation[
  match(as.character(seu$seurat_clusters), cluster_annot$cluster)
]
names(new_annot) <- colnames(seu)
seu <- AddMetaData(seu, metadata = new_annot, col.name = "annotation")
seu$annotation <- factor(seu$annotation, levels = cluster_annot$annotation)
message("Annotated ", sum(!is.na(seu$annotation)), " / ", ncol(seu), " cells")

# Cluster x stage table
write.csv(as.data.frame.matrix(xtab),
          file.path(tbl_dir, "06_headline_clusterxstage.csv"))

# ------------------------------------------------------------
# 3. Annotated headline UMAP
# ------------------------------------------------------------
cluster_palette <- c(
  "#7d8c6e", "#b88a3e", "#b15835", "#7c5c6b",
  "#c08d4f", "#94634a", "#5c7060", "#9b3d2a",
  "#a8763e", "#6b8a72", "#8a5570"
)[seq_along(unique(seu$seurat_clusters))]
names(cluster_palette) <- as.character(sort(unique(seu$seurat_clusters)))

p_umap_stage <- DimPlot(seu, reduction = "umap", group.by = "stage",
                        cols = stage_colors, pt.size = 0.3) +
  theme_portfolio() + theme(legend.position = "right") +
  labs(title = "Headline UMAP: by stage")

p_umap_annot <- DimPlot(seu, reduction = "umap", group.by = "annotation",
                        label = TRUE, label.size = 3, pt.size = 0.3,
                        repel = TRUE) +
  theme_portfolio() + theme(legend.position = "right") +
  labs(title = "Headline UMAP: annotated clusters")

annotated_fig <- p_umap_stage | p_umap_annot
ggsave(file.path(fig_dir, "13_annotated_umap.pdf"),
       annotated_fig, width = 14, height = 5.5)
ggsave(file.path(fig_dir, "13_annotated_umap.png"),
       annotated_fig, width = 14, height = 5.5, dpi = 200)

# ------------------------------------------------------------
# 4. Canonical marker dot plots
# ------------------------------------------------------------
markers <- list(
  Memory_B    = c("CD19", "MS4A1", "CD27", "IGHM", "IGHD", "CXCR5"),
  Activation  = c("MKI67", "TOP2A", "AICDA"),
  PrePB_PB    = c("XBP1", "IRF4", "PRDM1", "MZB1"),
  PC_terminal = c("SDC1", "TNFRSF17", "CD38", "JCHAIN"),
  UPR         = c("ATF6", "HSPA5", "DDIT3")
)

all_markers <- unique(unlist(markers))
all_markers <- all_markers[all_markers %in% rownames(seu)]
missing <- setdiff(unique(unlist(markers)), all_markers)
message("Canonical markers present: ", length(all_markers), " / ",
        length(unique(unlist(markers))))
if (length(missing) > 0) message("Missing: ", paste(missing, collapse = ", "))

DefaultAssay(seu) <- "RNA"
p_markers_stage <- DotPlot(seu, features = all_markers, group.by = "stage",
                            cols = c("#f3ede3", "#b15835"), dot.scale = 7) +
  theme_portfolio() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Canonical B-cell -> PC markers by stage", x = NULL, y = NULL)

p_markers_cluster <- DotPlot(seu, features = all_markers,
                              group.by = "annotation",
                              cols = c("#f3ede3", "#b15835"), dot.scale = 7) +
  theme_portfolio() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Same markers by annotated cluster", x = NULL, y = NULL)

markers_combined <- p_markers_stage / p_markers_cluster +
  plot_annotation(
    title = "Marker gene validation",
    subtitle = "Top: by experimental stage label | Bottom: by annotated cluster",
    theme = theme(
      plot.title = element_text(face = "bold", size = 15, color = "#2a2622"),
      plot.subtitle = element_text(size = 11, color = "#5c544b")
    )
  )

ggsave(file.path(fig_dir, "14_marker_dotplot.pdf"),
       markers_combined, width = 12, height = 8)
ggsave(file.path(fig_dir, "14_marker_dotplot.png"),
       markers_combined, width = 12, height = 8, dpi = 200)

# ------------------------------------------------------------
# 5. FeaturePlot of showcase markers
# ------------------------------------------------------------
showcase_markers <- intersect(
  c("CD19", "CD27", "MKI67", "XBP1", "PRDM1", "IRF4",
    "MZB1", "SDC1", "TNFRSF17"),
  rownames(seu)
)

p_feat <- FeaturePlot(seu, features = showcase_markers, reduction = "umap",
                      cols = c("#f3ede3", "#b15835"),
                      ncol = 3, pt.size = 0.3) &
  theme_portfolio() & NoAxes()

ggsave(file.path(fig_dir, "15_featureplot_markers.pdf"),
       p_feat, width = 12, height = ceiling(length(showcase_markers) / 3) * 3.5)
ggsave(file.path(fig_dir, "15_featureplot_markers.png"),
       p_feat, width = 12, height = ceiling(length(showcase_markers) / 3) * 3.5,
       dpi = 200)

# ------------------------------------------------------------
# 6. Find ALL markers per cluster (full, not truncated)
# ------------------------------------------------------------
message("\nFinding markers per cluster (~2-5 min)...")
Idents(seu) <- "seurat_clusters"
top_markers_all <- FindAllMarkers(seu, only.pos = TRUE,
                                   min.pct = 0.25,
                                   logfc.threshold = 0.25,
                                   verbose = FALSE)

# Top N per cluster - we keep more this time (top 30)
top30_per_cluster <- top_markers_all %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 30) %>%
  ungroup()

write.csv(top_markers_all,
          file.path(tbl_dir, "06_all_markers_per_cluster.csv"),
          row.names = FALSE)
write.csv(top30_per_cluster,
          file.path(tbl_dir, "06_top30_markers_per_cluster.csv"),
          row.names = FALSE)

# Print top 10 per cluster, BUT not truncated by tibble
message("\nTop 10 markers per cluster (full table):")
for (cl in levels(top30_per_cluster$cluster)) {
  cat("\n--- Cluster ", cl, " ---\n", sep = "")
  top10 <- top30_per_cluster %>%
    filter(cluster == cl) %>%
    slice_head(n = 10) %>%
    select(gene, avg_log2FC, p_val_adj)
  print(as.data.frame(top10))
}

# ------------------------------------------------------------
# 7. Top marker heatmap
# ------------------------------------------------------------
# Take top 5 per cluster, plot a scaled-expression heatmap
# averaged per cluster. This is a standard, publication-quality
# figure to summarize cluster identity.
# ------------------------------------------------------------
top5_per_cluster <- top_markers_all %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 10) %>%
  ungroup()
top5_genes <- unique(top5_per_cluster$gene)
top5_genes <- top5_genes[top5_genes %in% rownames(seu)]

# Average expression per cluster
avg_expr <- AverageExpression(seu, features = top5_genes,
                               group.by = "annotation",
                               assays = "RNA")$RNA
# Scale across clusters
avg_scaled <- t(scale(t(as.matrix(avg_expr))))
avg_scaled[avg_scaled > 2.5]  <- 2.5
avg_scaled[avg_scaled < -2.5] <- -2.5

# Build pheatmap with portfolio palette
heat_colors <- colorRampPalette(c("#7d8c6e", "#f3ede3", "#b15835"))(50)

pdf(file.path(fig_dir, "16_top_marker_heatmap.pdf"),
    width = 7, height = 11)
pheatmap(avg_scaled,
         color = heat_colors,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize_row = 8,
         fontsize_col = 9,
         border_color = "white",
         main = "Top 5 marker genes per cluster (scaled, average expression)",
         treeheight_row = 25,
         treeheight_col = 20)
dev.off()

png(file.path(fig_dir, "16_top_marker_heatmap.png"),
    width = 7, height = 11, units = "in", res = 200)
pheatmap(avg_scaled,
         color = heat_colors,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize_row = 8,
         fontsize_col = 9,
         border_color = "white",
         main = "Top 10 marker genes per cluster (scaled, average expression)",
         treeheight_row = 25,
         treeheight_col = 20)
dev.off()

# ------------------------------------------------------------
# 8. Focused comparison: the two PC clusters
# ------------------------------------------------------------
# We have two PC-dominated clusters. What distinguishes them?
# ------------------------------------------------------------
pc_clusters <- cluster_annot$cluster[cluster_annot$dominant_stage == "PC"]
message("\nPC clusters: ", paste(pc_clusters, collapse = ", "))

if (length(pc_clusters) >= 2) {
  message("\nFinding DEGs between PC clusters ",
          pc_clusters[1], " vs ", pc_clusters[2], "...")

  pc_de <- FindMarkers(seu,
                       ident.1 = pc_clusters[1],
                       ident.2 = pc_clusters[2],
                       min.pct = 0.25,
                       logfc.threshold = 0.25,
                       verbose = FALSE)
  pc_de$gene <- rownames(pc_de)
  pc_de <- pc_de %>% arrange(p_val_adj, desc(abs(avg_log2FC)))

  write.csv(pc_de,
            file.path(tbl_dir, "06_PCcluster_comparison.csv"),
            row.names = FALSE)

  message("\nTop 15 genes enriched in C", pc_clusters[1],
          " vs C", pc_clusters[2], ":")
  print(as.data.frame(
    pc_de %>% filter(avg_log2FC > 0) %>% slice_head(n = 15) %>%
      select(gene, avg_log2FC, p_val_adj)
  ))

  message("\nTop 15 genes enriched in C", pc_clusters[2],
          " vs C", pc_clusters[1], ":")
  print(as.data.frame(
    pc_de %>% filter(avg_log2FC < 0) %>% slice_head(n = 15) %>%
      select(gene, avg_log2FC, p_val_adj)
  ))

  # Quick volcano plot
  pc_de$direction <- "n.s."
  pc_de$direction[pc_de$avg_log2FC >  0.5 & pc_de$p_val_adj < 0.01] <- paste0("Up in C", pc_clusters[1])
  pc_de$direction[pc_de$avg_log2FC < -0.5 & pc_de$p_val_adj < 0.01] <- paste0("Up in C", pc_clusters[2])

  p_volcano <- ggplot(pc_de,
                      aes(x = avg_log2FC, y = -log10(p_val_adj + 1e-300),
                          color = direction)) +
    geom_point(alpha = 0.6, size = 1.2) +
    scale_color_manual(values = setNames(
      c("#b15835", "#7d8c6e", "#cccccc"),
      c(paste0("Up in C", pc_clusters[1]),
        paste0("Up in C", pc_clusters[2]),
        "n.s."))) +
    geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed",
               color = "#8a7f73", linewidth = 0.3) +
    geom_hline(yintercept = -log10(0.01), linetype = "dashed",
               color = "#8a7f73", linewidth = 0.3) +
    theme_portfolio() +
    labs(title = paste0("Differential expression: C", pc_clusters[1],
                        " vs C", pc_clusters[2]),
         subtitle = "Two PC-dominated clusters from the headline UMAP",
         x = "log2 fold change",
         y = "-log10(adjusted p-value)",
         color = NULL) +
    theme(legend.position = "top")

  ggsave(file.path(fig_dir, "17_PCcluster_volcano.pdf"),
         p_volcano, width = 8, height = 6)
  ggsave(file.path(fig_dir, "17_PCcluster_volcano.png"),
         p_volcano, width = 8, height = 6, dpi = 200)
}

# ------------------------------------------------------------
# 9. Save annotated object
# ------------------------------------------------------------
saveRDS(seu, file.path(out_dir, "seurat_integrated_annotated.rds"))

message("\n========================================")
message("Stage 2.5 (Marker validation) complete!")
message("\nFigures:")
message("  - 13_annotated_umap")
message("  - 14_marker_dotplot")
message("  - 15_featureplot_markers")
message("  - 16_top_marker_heatmap")
message("  - 17_PCcluster_volcano   (the two PC clusters compared)")
message("\nTables:")
message("  - 06_cluster_annotations.csv")
message("  - 06_headline_clusterxstage.csv")
message("  - 06_all_markers_per_cluster.csv")
message("  - 06_top30_markers_per_cluster.csv")
message("  - 06_PCcluster_comparison.csv")
message("========================================")
