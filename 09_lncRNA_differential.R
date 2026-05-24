# ============================================================
# 09_lncRNA_differential.R
#
# Differential expression of lncRNAs across:
#   1. Differentiation stages (MBC, prePB, PB, PC)
#   2. Annotated clusters (the 6 clusters from Phase 1)
#
# Builds on:
#   - data/processed/seurat_integrated_annotated.rds  (the headline object)
#   - data/processed/lncRNA_detected_symbols.rds      (from script 08)
#
# Outputs:
#   - DE tables for stage and cluster comparisons
#   - Volcano plots for the most informative contrasts
#   - Heatmap of top stage-marker lncRNAs
#   - Annotated UMAP highlighting top lncRNA expression patterns
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(pheatmap)
  library(ggrepel)
  library(Matrix)
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
      plot.title    = element_text(face = "bold", size = 12, color = "#2a2622"),
      plot.subtitle = element_text(size = 10, color = "#5c544b"),
      axis.text     = element_text(color = "#2a2622", size = 9),
      axis.title    = element_text(color = "#5c544b", size = 10),
      panel.grid.minor = element_blank()
    )
}

# ------------------------------------------------------------
# 1. Load object and lncRNA panel
# ------------------------------------------------------------
seu <- readRDS(file.path(out_dir, "seurat_integrated_annotated.rds"))
seu$stage <- factor(seu$stage, levels = c("MBC", "prePB", "PB", "PC"))
DefaultAssay(seu) <- "RNA"

detected_lncRNAs <- readRDS(file.path(out_dir, "lncRNA_detected_symbols.rds"))
message("Loaded: ", ncol(seu), " cells, ", length(detected_lncRNAs), " lncRNAs detected")

# Restrict downstream tests to lncRNAs only
# This is faster than testing all genes and then filtering.

# ------------------------------------------------------------
# 2. Stage-level differential expression (one vs rest)
# ------------------------------------------------------------
# For each of the 4 stages, find lncRNAs enriched vs all
# other stages. This identifies stage-specific lncRNA markers.
# ------------------------------------------------------------
message("\n=== Stage-level DE (one stage vs all others) ===")

Idents(seu) <- "stage"

stage_de_results <- list()
for (s in levels(seu$stage)) {
  message("  Testing stage: ", s)
  de_s <- FindMarkers(
    seu,
    ident.1 = s,
    features = detected_lncRNAs,
    only.pos = FALSE,
    min.pct = 0.10,
    logfc.threshold = 0.25,
    verbose = FALSE
  )
  if (nrow(de_s) > 0) {
    de_s$gene  <- rownames(de_s)
    de_s$stage <- s
    stage_de_results[[s]] <- de_s
  }
}

stage_de <- bind_rows(stage_de_results) %>%
  arrange(stage, p_val_adj)

write.csv(stage_de,
          file.path(tbl_dir, "08_lncRNA_DE_per_stage.csv"),
          row.names = FALSE)

# Top hits per stage
message("\nTop 5 stage-enriched lncRNAs per stage:")
top5_per_stage <- stage_de %>%
  filter(avg_log2FC > 0) %>%
  group_by(stage) %>%
  arrange(p_val_adj) %>%
  slice_head(n = 5) %>%
  select(stage, gene, avg_log2FC, p_val_adj, pct.1, pct.2)
print(as.data.frame(top5_per_stage))

# ------------------------------------------------------------
# 3. Stage-marker lncRNA heatmap
# ------------------------------------------------------------
# Top 10 enriched lncRNAs per stage -> heatmap of average
# expression per stage.
# ------------------------------------------------------------
top10_per_stage <- stage_de %>%
  filter(avg_log2FC > 0, p_val_adj < 0.01) %>%
  group_by(stage) %>%
  arrange(p_val_adj) %>%
  slice_head(n = 10) %>%
  ungroup()

heatmap_genes <- unique(top10_per_stage$gene)
heatmap_genes <- heatmap_genes[heatmap_genes %in% rownames(seu)]
message("\nHeatmap genes (unique top10 per stage): ", length(heatmap_genes))

if (length(heatmap_genes) >= 3) {
  avg_expr <- AverageExpression(
    seu, features = heatmap_genes,
    group.by = "stage", assays = "RNA"
  )$RNA
  avg_scaled <- t(scale(t(as.matrix(avg_expr))))
  avg_scaled[avg_scaled >  2.5] <-  2.5
  avg_scaled[avg_scaled < -2.5] <- -2.5

  heat_colors <- colorRampPalette(c("#7d8c6e", "#f3ede3", "#b15835"))(50)

  pdf(file.path(fig_dir, "20_lncRNA_stage_heatmap.pdf"),
      width = 6, height = max(7, length(heatmap_genes) * 0.22))
  pheatmap(avg_scaled,
           color = heat_colors,
           cluster_rows = TRUE, cluster_cols = FALSE,
           show_rownames = TRUE, show_colnames = TRUE,
           fontsize_row = 7.5, fontsize_col = 10,
           border_color = "white",
           main = "Stage-enriched lncRNAs (top 10 per stage)",
           treeheight_row = 25)
  dev.off()

  png(file.path(fig_dir, "20_lncRNA_stage_heatmap.png"),
      width = 6, height = max(7, length(heatmap_genes) * 0.22),
      units = "in", res = 200)
  pheatmap(avg_scaled,
           color = heat_colors,
           cluster_rows = TRUE, cluster_cols = FALSE,
           show_rownames = TRUE, show_colnames = TRUE,
           fontsize_row = 7.5, fontsize_col = 10,
           border_color = "white",
           main = "Stage-enriched lncRNAs (top 10 per stage)",
           treeheight_row = 25)
  dev.off()
}

# ------------------------------------------------------------
# 4. Cluster-level DE (lncRNAs marking the 6 clusters)
# ------------------------------------------------------------
message("\n=== Cluster-level DE (annotated clusters from Phase 1) ===")

Idents(seu) <- "seurat_clusters"

cluster_de_results <- list()
for (cl in levels(seu$seurat_clusters)) {
  message("  Testing cluster: ", cl)
  de_cl <- FindMarkers(
    seu,
    ident.1 = cl,
    features = detected_lncRNAs,
    only.pos = TRUE,         # we want CLUSTER-MARKER lncRNAs
    min.pct = 0.10,
    logfc.threshold = 0.25,
    verbose = FALSE
  )
  if (nrow(de_cl) > 0) {
    de_cl$gene    <- rownames(de_cl)
    de_cl$cluster <- cl
    cluster_de_results[[as.character(cl)]] <- de_cl
  }
}

cluster_de <- bind_rows(cluster_de_results) %>%
  arrange(cluster, p_val_adj)

write.csv(cluster_de,
          file.path(tbl_dir, "08_lncRNA_DE_per_cluster.csv"),
          row.names = FALSE)

message("\nTop 5 cluster-marker lncRNAs per cluster:")
top5_per_cluster <- cluster_de %>%
  group_by(cluster) %>%
  arrange(p_val_adj) %>%
  slice_head(n = 5) %>%
  select(cluster, gene, avg_log2FC, p_val_adj, pct.1, pct.2)
print(as.data.frame(top5_per_cluster))

# ------------------------------------------------------------
# 5. Volcano plots, stage-level (one per stage)
# ------------------------------------------------------------
make_volcano <- function(de_data, stage_label, color_up, color_down = "#7d8c6e") {
  de_data <- de_data %>%
    mutate(
      p_plot = pmax(p_val_adj, 1e-300),
      neg_log10_p = -log10(p_plot),
      direction = case_when(
        avg_log2FC >  0.5 & p_val_adj < 0.01 ~ paste0("Up in ", stage_label),
        avg_log2FC < -0.5 & p_val_adj < 0.01 ~ paste0("Down in ", stage_label),
        TRUE ~ "n.s."
      )
    )

  to_label <- de_data %>%
    filter(direction != "n.s.") %>%
    arrange(p_val_adj) %>%
    slice_head(n = 10)

  cat_colors <- setNames(
    c(color_up, color_down, "#cfc7b8"),
    c(paste0("Up in ", stage_label),
      paste0("Down in ", stage_label),
      "n.s.")
  )

  ggplot(de_data,
         aes(x = avg_log2FC, y = neg_log10_p, color = direction)) +
    geom_vline(xintercept = c(-0.5, 0.5),
               linetype = "dashed", color = "#8a7f73", linewidth = 0.3) +
    geom_hline(yintercept = -log10(0.01),
               linetype = "dashed", color = "#8a7f73", linewidth = 0.3) +
    geom_point(data = de_data %>% filter(direction == "n.s."),
               alpha = 0.3, size = 1) +
    geom_point(data = de_data %>% filter(direction != "n.s."),
               alpha = 0.85, size = 1.6) +
    ggrepel::geom_text_repel(
      data = to_label, aes(label = gene),
      size = 2.8, max.overlaps = 20, box.padding = 0.4,
      color = "#2a2622", fontface = "italic", show.legend = FALSE
    ) +
    scale_color_manual(values = cat_colors, name = NULL) +
    theme_portfolio() +
    labs(
      title = paste0("lncRNAs differentially expressed in ", stage_label),
      x = "log2 fold change",
      y = "-log10(adjusted p-value)"
    ) +
    theme(legend.position = "top")
}

volcano_plots <- list()
for (s in levels(seu$stage)) {
  de_s <- stage_de %>% filter(stage == s)
  if (nrow(de_s) > 5) {
    volcano_plots[[s]] <- make_volcano(de_s, s, color_up = stage_colors[s])
  }
}

# Combined 2x2 volcano panel
if (length(volcano_plots) == 4) {
  combined_volcano <- (volcano_plots$MBC | volcano_plots$prePB) /
                      (volcano_plots$PB  | volcano_plots$PC) +
    plot_annotation(
      title = "lncRNA differential expression across plasma cell differentiation",
      subtitle = "One stage vs. all other stages combined",
      theme = theme(
        plot.title = element_text(face = "bold", size = 14, color = "#2a2622"),
        plot.subtitle = element_text(size = 11, color = "#5c544b")
      )
    )

  ggsave(file.path(fig_dir, "21_lncRNA_volcanos_per_stage.pdf"),
         combined_volcano, width = 14, height = 12)
  ggsave(file.path(fig_dir, "21_lncRNA_volcanos_per_stage.png"),
         combined_volcano, width = 14, height = 12, dpi = 200)
}

# ------------------------------------------------------------
# 6. Highlight a few "candidate hit" lncRNAs on the UMAP
# ------------------------------------------------------------
# Take the strongest stage-marker lncRNA per stage and show
# it as a FeaturePlot, so the reader can see the spatial
# pattern in the headline UMAP.
# ------------------------------------------------------------
candidate_lncRNAs <- top10_per_stage %>%
  filter(avg_log2FC > 0) %>%
  group_by(stage) %>%
  arrange(p_val_adj) %>%
  slice_head(n = 1) %>%
  pull(gene)

candidate_lncRNAs <- intersect(candidate_lncRNAs, rownames(seu))
message("\nCandidate lncRNAs for UMAP showcase: ",
        paste(candidate_lncRNAs, collapse = ", "))

if (length(candidate_lncRNAs) >= 3) {
  p_feat <- FeaturePlot(seu, features = candidate_lncRNAs,
                        reduction = "umap",
                        cols = c("#f3ede3", "#b15835"),
                        ncol = 2, pt.size = 0.3) &
    theme_portfolio() & NoAxes()

  ggsave(file.path(fig_dir, "22_lncRNA_candidates_umap.pdf"),
         p_feat, width = 10,
         height = ceiling(length(candidate_lncRNAs) / 2) * 3.5)
  ggsave(file.path(fig_dir, "22_lncRNA_candidates_umap.png"),
         p_feat, width = 10,
         height = ceiling(length(candidate_lncRNAs) / 2) * 3.5,
         dpi = 200)
}

# ------------------------------------------------------------
# 7. Done
# ------------------------------------------------------------
message("\n========================================")
message("Script 09 (lncRNA differential expression) complete!")
message("\nOutputs:")
message("  Tables:")
message("    08_lncRNA_DE_per_stage.csv")
message("    08_lncRNA_DE_per_cluster.csv")
message("  Figures:")
message("    20_lncRNA_stage_heatmap      (top lncRNAs per stage)")
message("    21_lncRNA_volcanos_per_stage (4-panel volcano)")
message("    22_lncRNA_candidates_umap    (FeaturePlot of top hits)")
message("\nNext: 10_hdWGCNA_setup.R")
message("========================================")
