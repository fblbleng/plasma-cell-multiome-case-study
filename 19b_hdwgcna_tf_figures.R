# ============================================================
# 19b_hdwgcna_tf_figures.R
#
# Figure generation for the hdWGCNA TF regulatory network.
# Run AFTER ConstructTFNetwork2 + AssignTFRegulons, with `seu`
# and `regulons` in memory (or load the final object below).
#
# Correction vs the original script: the regulon importance
# column is 'Gain' (xgboost gain), not 'score'.
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(hdWGCNA)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

if (requireNamespace("conflicted", quietly = TRUE)) {
  conflicted::conflicts_prefer(dplyr::filter, .quiet = TRUE)
  conflicted::conflicts_prefer(dplyr::select, .quiet = TRUE)
  conflicted::conflicts_prefer(dplyr::mutate, .quiet = TRUE)
  conflicted::conflicts_prefer(dplyr::desc, .quiet = TRUE)
  conflicted::conflicts_prefer(dplyr::slice, .quiet = TRUE)
  conflicted::conflicts_prefer(stats::cor, .quiet = TRUE)
  conflicted::conflicts_prefer(base::intersect, .quiet = TRUE)
}

fig_dir <- "results/figures"
tbl_dir <- "results/tables"
out_dir <- "data/processed"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tbl_dir, showWarnings = FALSE, recursive = TRUE)

# If not already in memory, load:
# seu <- readRDS(file.path(out_dir, "seurat_hdwgcna_tf_network.rds"))
# regulons <- GetTFRegulons(seu)

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

module_colors <- c(
  "yellow"    = "#e8b800",
  "brown"     = "#964b00",
  "turquoise" = "#4cb5b0",
  "blue"      = "#4472c4"
)

modules <- GetModules(seu)

# ------------------------------------------------------------
# Regulon summary (corrected: use Gain, not score)
# ------------------------------------------------------------
# Detect the importance column name robustly
imp_col <- intersect(c("Gain", "score", "importance"), colnames(regulons))[1]
message("Using importance column: ", imp_col)

regulon_summary <- regulons |>
  dplyr::group_by(tf) |>
  dplyr::summarize(
    n_targets = dplyr::n(),
    mean_gain = mean(.data[[imp_col]], na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::arrange(dplyr::desc(n_targets))

write.csv(regulon_summary,
          file.path(tbl_dir, "17_regulon_summary.csv"),
          row.names = FALSE)

# ------------------------------------------------------------
# Figure 46: top 30 TFs by regulon size, colored by mean gain
# ------------------------------------------------------------
message("Figure 46: top TFs by regulon size...")

top30 <- head(regulon_summary, 30)

p46 <- ggplot(top30,
              aes(x = reorder(tf, n_targets), y = n_targets, fill = mean_gain)) +
  geom_col(width = 0.75) +
  coord_flip() +
  scale_fill_gradient(low = "#e9d8c4", high = "#b15835",
                      name = "Mean\nxgboost gain") +
  theme_portfolio() +
  labs(
    title    = "Top 30 transcription factors by regulon size",
    subtitle = "Bar length: number of predicted targets. Color: mean predictive strength (xgboost gain).",
    x = NULL, y = "Number of target genes"
  )

ggsave(file.path(fig_dir, "46_tf_count_summary.pdf"), p46,
       width = 8, height = 9)
ggsave(file.path(fig_dir, "46_tf_count_summary.png"), p46,
       width = 8, height = 9, dpi = 200)

# ------------------------------------------------------------
# Figure 47: regulon bar plots for key TFs (cross-modal)
# ------------------------------------------------------------
message("Figure 47: key TF regulon bar plots...")

candidate_tfs <- c("PRDM1", "BATF3", "IRF1", "IRF7", "SPIB",
                   "STAT1", "STAT3", "KLF6", "KLF2", "EGR1",
                   "POU2F2", "XBP1", "NFKB1", "ETV4", "TFDP1")
available_tfs <- intersect(candidate_tfs, unique(regulons$tf))
message("  Key TFs available: ", paste(available_tfs, collapse = ", "))

if (length(available_tfs) >= 2) {
  display_tfs <- head(available_tfs, 6)
  bar_plots <- list()
  for (tf_name in display_tfs) {
    tryCatch({
      p <- RegulonBarPlot(seu, selected_tf = tf_name) +
        theme_portfolio() +
        theme(plot.title = element_text(size = 10))
      bar_plots[[tf_name]] <- p
    }, error = function(e) message("  Skipping ", tf_name, ": ", e$message))
  }
  if (length(bar_plots) > 0) {
    n_rows <- ceiling(length(bar_plots) / 2)
    combined <- wrap_plots(bar_plots, ncol = 2) +
      plot_annotation(
        title    = "Top regulons of key transcription factors",
        subtitle = "Target genes ranked by xgboost importance, signed by correlation"
      )
    ggsave(file.path(fig_dir, "47_top_regulons_barplots.pdf"),
           combined, width = 12, height = 4 * n_rows)
    ggsave(file.path(fig_dir, "47_top_regulons_barplots.png"),
           combined, width = 12, height = 4 * n_rows, dpi = 200)
  }
}

# ------------------------------------------------------------
# Figure 48: TF regulons x hdWGCNA module heatmap (headline)
# ------------------------------------------------------------
message("Figure 48: TF x module overlap heatmap...")

regulons_with_module <- regulons |>
  dplyr::left_join(
    modules |> dplyr::select(gene_name, module),
    by = c("gene" = "gene_name")
  ) |>
  dplyr::filter(!is.na(module), module != "grey")

overlap_summary <- regulons_with_module |>
  dplyr::group_by(tf, module) |>
  dplyr::summarize(n_targets = dplyr::n(), .groups = "drop")

tfs_with_signal <- overlap_summary |>
  dplyr::group_by(tf) |>
  dplyr::summarize(max_t = max(n_targets), .groups = "drop") |>
  dplyr::filter(max_t >= 5) |>
  dplyr::pull(tf)

overlap_filtered <- overlap_summary |>
  dplyr::filter(tf %in% tfs_with_signal)

write.csv(overlap_filtered,
          file.path(tbl_dir, "17_regulon_module_overlap.csv"),
          row.names = FALSE)
message("  TFs with >= 5 targets in any module: ", length(tfs_with_signal))

if (nrow(overlap_filtered) > 0) {
  top_tfs <- overlap_filtered |>
    dplyr::group_by(tf) |>
    dplyr::summarize(total = sum(n_targets), .groups = "drop") |>
    dplyr::arrange(dplyr::desc(total)) |>
    dplyr::slice_head(n = 25) |>
    dplyr::pull(tf)
  
  heatmap_df <- overlap_filtered |>
    dplyr::filter(tf %in% top_tfs) |>
    dplyr::mutate(module = factor(module,
                                  levels = c("yellow", "brown", "turquoise", "blue")))
  
  tf_order <- heatmap_df |>
    dplyr::group_by(tf) |>
    dplyr::slice_max(n_targets, n = 1, with_ties = FALSE) |>
    dplyr::arrange(module, dplyr::desc(n_targets)) |>
    dplyr::pull(tf)
  heatmap_df$tf <- factor(heatmap_df$tf, levels = rev(tf_order))
  
  p48 <- ggplot(heatmap_df, aes(x = module, y = tf, fill = n_targets)) +
    geom_tile(color = "white", linewidth = 0.4) +
    geom_text(aes(label = n_targets), color = "#2a2622", size = 3) +
    scale_fill_gradient(low = "#faf7f2", high = "#b15835",
                        name = "Targets in\nmodule") +
    theme_minimal(base_size = 9) +
    theme(
      plot.title    = element_text(face = "bold", size = 12, color = "#2a2622"),
      plot.subtitle = element_text(size = 10, color = "#5c544b"),
      axis.text.y   = element_text(size = 9),
      panel.grid    = element_blank()
    ) +
    labs(
      title    = "TF regulons mapped to hdWGCNA modules",
      subtitle = "Each cell: number of a TF's predicted targets in each co-expression module",
      x = NULL, y = NULL
    )
  
  ggsave(file.path(fig_dir, "48_regulon_module_overlap.pdf"),
         p48, width = 7, height = max(7, length(top_tfs) * 0.3))
  ggsave(file.path(fig_dir, "48_regulon_module_overlap.png"),
         p48, width = 7, height = max(7, length(top_tfs) * 0.3), dpi = 200)
}

message("\n========================================")
message("Figures complete:")
message("  46_tf_count_summary.png")
message("  47_top_regulons_barplots.png")
message("  48_regulon_module_overlap.png")
message("\nTransfer the 3 PNGs to fblbleng.github.io/images/case-study/")
message("========================================")