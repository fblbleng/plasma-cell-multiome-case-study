# ============================================================
# 15b_motif_figures_refined.R
#
# Refined visualizations of the motif enrichment results.
# Reuses enrichment_full from script 15 (or loads from CSV)
# without re-running AddMotifs or FindMotifs.
#
# Changes from script 15:
#   - Heatmap: row z-scored to show RELATIVE enrichment across
#     transitions (instead of absolute, which is dominated by
#     housekeeping TFs like KLF/SP/NRF1)
#   - Heatmap: selects top motifs by SPECIFICITY (which motifs
#     differ most between transitions), not by absolute enrichment
#   - Volcanos: top 10 labeled per transition with smart
#     repulsion to avoid label collision
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
  library(ggrepel)
  library(scales)
})

if (requireNamespace("conflicted", quietly = TRUE)) {
  conflicted::conflicts_prefer(dplyr::filter, .quiet = TRUE)
  conflicted::conflicts_prefer(dplyr::select, .quiet = TRUE)
  conflicted::conflicts_prefer(dplyr::mutate, .quiet = TRUE)
  conflicted::conflicts_prefer(dplyr::desc, .quiet = TRUE)
}

set.seed(42)

fig_dir <- "results/figures"
tbl_dir <- "results/tables"

transitions <- list(
  T1_MBC_to_prePB = list(from = "MBC",   to = "prePB"),
  T2_prePB_to_PB  = list(from = "prePB", to = "PB"),
  T3_PB_to_PC     = list(from = "PB",    to = "PC")
)

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
# 1. Load enrichment results
# ------------------------------------------------------------
# Try to use the in-memory object first; fall back to CSV
if (!exists("enrichment_full")) {
  enrichment_full <- read.csv(
    file.path(tbl_dir, "14_motif_enrichment_per_transition.csv"),
    stringsAsFactors = FALSE
  )
  message("Loaded enrichment_full from CSV")
} else {
  message("Using enrichment_full from current R session")
}

# ------------------------------------------------------------
# 2. Refined heatmap: row z-scored, transition-specific motifs
# ------------------------------------------------------------
message("\nBuilding refined heatmap...")

# Pivot to motif x transition (opening only)
enrich_open <- enrichment_full |>
  dplyr::filter(direction == "opening") |>
  dplyr::select(motif.name, transition, log2_fold_enrichment, pvalue)

# Wide form
fe_wide <- enrich_open |>
  dplyr::select(motif.name, transition, log2_fold_enrichment) |>
  tidyr::pivot_wider(names_from = transition,
                     values_from = log2_fold_enrichment,
                     values_fill = 0)

fe_mat <- as.matrix(fe_wide[, -1])
rownames(fe_mat) <- fe_wide$motif.name

# Drop motifs that aren't enriched anywhere
fe_mat <- fe_mat[apply(fe_mat, 1, max) > log2(1.5), ]

# Row z-score
fe_zscored <- t(scale(t(fe_mat)))
fe_zscored[is.na(fe_zscored)] <- 0   # constant rows

# Compute specificity: variance across transitions
# Higher variance = more transition-specific
motif_specificity <- apply(fe_zscored, 1, function(x) {
  sd(x, na.rm = TRUE)
})

# Get top motifs PER transition by z-score (positive deviation = enriched
# in that transition relative to others)
top_per_t <- list()
for (tname in colnames(fe_zscored)) {
  # Order by z-score for this transition, take top 10
  motifs_ranked <- names(sort(fe_zscored[, tname], decreasing = TRUE))
  top_per_t[[tname]] <- head(motifs_ranked, 12)
}

# Union of top specific motifs across all transitions
selected_motifs <- unique(unlist(top_per_t))
message("Heatmap motifs selected (transition-specific): ",
        length(selected_motifs))

# Subset the z-scored matrix
fe_zscored_sub <- fe_zscored[selected_motifs, , drop = FALSE]

# Order rows: assign each motif to its top transition, then sort within
top_transition <- apply(fe_zscored_sub, 1, which.max)
fe_zscored_sub <- fe_zscored_sub[
  order(top_transition, -apply(fe_zscored_sub, 1, max)), ,
  drop = FALSE
]

# Long format for ggplot
heatmap_df <- fe_zscored_sub |>
  as.data.frame() |>
  tibble::rownames_to_column("motif") |>
  tidyr::pivot_longer(-motif, names_to = "transition", values_to = "z")

heatmap_df$motif <- factor(heatmap_df$motif,
                           levels = rev(rownames(fe_zscored_sub)))
heatmap_df$transition <- factor(heatmap_df$transition,
                                levels = names(transitions),
                                labels = c("MBC -> prePB",
                                           "prePB -> PB",
                                           "PB -> PC"))

# Add an annotation for the most-enriched transition per motif
top_annot <- data.frame(
  motif = rownames(fe_zscored_sub),
  top_in = factor(names(transitions)[top_transition[rownames(fe_zscored_sub)]],
                  levels = names(transitions))
)
heatmap_df <- heatmap_df |>
  dplyr::left_join(top_annot, by = "motif")

# Better color: diverging scale for z-scores
p_heatmap <- ggplot(heatmap_df,
                    aes(x = transition, y = motif, fill = z)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.2f", z)),
            color = "#2a2622", size = 2.5) +
  scale_fill_gradient2(low = "#3b6286", mid = "#faf7f2", high = "#b15835",
                       midpoint = 0, name = "z-score",
                       limits = c(-1.3, 1.3),
                       oob = scales::squish) +
  theme_minimal(base_size = 9) +
  theme(
    plot.title    = element_text(face = "bold", size = 12, color = "#2a2622"),
    plot.subtitle = element_text(size = 10, color = "#5c544b"),
    axis.text.y   = element_text(size = 7.5),
    axis.text.x   = element_text(angle = 0, hjust = 0.5),
    panel.grid    = element_blank()
  ) +
  labs(title    = "Transition-specific motif enrichment",
       subtitle = "Row z-scored log2 fold enrichment. High z = relatively more enriched in that transition.",
       x = NULL, y = NULL)

ggsave(file.path(fig_dir, "38b_motif_heatmap_zscored.pdf"),
       p_heatmap, width = 7,
       height = max(8, length(selected_motifs) * 0.20))
ggsave(file.path(fig_dir, "38b_motif_heatmap_zscored.png"),
       p_heatmap, width = 7,
       height = max(8, length(selected_motifs) * 0.20),
       dpi = 200)
message("Saved: 38b_motif_heatmap_zscored")

# ------------------------------------------------------------
# 3. Refined volcanos: top 10 labels per transition with ggrepel
# ------------------------------------------------------------
message("\nBuilding refined volcanos with smart labels...")

# Add log10p to the data (cap at 50 for visualization)
enrich_with_p <- enrichment_full |>
  dplyr::filter(direction == "opening") |>
  dplyr::mutate(
    neg_log10p = pmin(-log10(pvalue + 1e-300), 50)
  )

volcano_list <- list()
for (tname in names(transitions)) {
  df <- enrich_with_p |>
    dplyr::filter(transition == tname)
  
  # Pick top 10 by combined score: significance AND magnitude
  # This avoids labelling only ones at the p-value ceiling
  top_labels <- df |>
    dplyr::filter(significant) |>
    dplyr::mutate(combined_score = -log10(pvalue + 1e-300) *
                    log2_fold_enrichment) |>
    dplyr::arrange(dplyr::desc(combined_score)) |>
    dplyr::slice_head(n = 10)
  
  p <- ggplot(df, aes(x = log2_fold_enrichment, y = neg_log10p)) +
    geom_point(aes(color = significant, alpha = significant),
               size = 1.3) +
    scale_color_manual(values = c("FALSE" = "gray75", "TRUE" = "#b15835")) +
    scale_alpha_manual(values = c("FALSE" = 0.35, "TRUE" = 0.85)) +
    ggrepel::geom_text_repel(
      data = top_labels,
      aes(label = motif.name),
      size = 3.0,
      color = "#2a2622",
      max.overlaps = 20,
      box.padding = 0.4,
      min.segment.length = 0.2,
      segment.color = "gray50",
      segment.size = 0.3,
      seed = 42
    ) +
    geom_vline(xintercept = log2(1.5), linetype = "dashed",
               color = "gray50", linewidth = 0.3) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed",
               color = "gray50", linewidth = 0.3) +
    theme_portfolio() +
    labs(title    = paste0(transitions[[tname]]$from, " -> ",
                           transitions[[tname]]$to),
         subtitle = paste0(sum(df$significant), " significant motifs"),
         x = "log2 fold enrichment", y = "-log10(p-value)") +
    theme(legend.position = "none")
  
  volcano_list[[tname]] <- p
}

volcano_combined <- wrap_plots(volcano_list, ncol = 3) +
  plot_annotation(
    title = "TF motif enrichment in opening peaks per transition",
    subtitle = "Top 10 labels per panel ranked by (-log10 p) x log2FE"
  )

ggsave(file.path(fig_dir, "39b_motif_volcano_per_transition.pdf"),
       volcano_combined, width = 15, height = 6)
ggsave(file.path(fig_dir, "39b_motif_volcano_per_transition.png"),
       volcano_combined, width = 15, height = 6, dpi = 200)
message("Saved: 39b_motif_volcano_per_transition")

# ------------------------------------------------------------
# 4. Summary of transition-specific TFs (table)
# ------------------------------------------------------------
specificity_table <- data.frame(
  motif = rownames(fe_zscored),
  T1_MBC_to_prePB = fe_zscored[, "T1_MBC_to_prePB"],
  T2_prePB_to_PB  = fe_zscored[, "T2_prePB_to_PB"],
  T3_PB_to_PC     = fe_zscored[, "T3_PB_to_PC"],
  specificity     = motif_specificity,
  top_transition  = names(transitions)[
    apply(fe_zscored, 1, which.max)
  ]
) |>
  dplyr::arrange(dplyr::desc(specificity))

write.csv(specificity_table,
          file.path(tbl_dir, "14b_motif_specificity_ranked.csv"),
          row.names = FALSE)

message("\nTop 15 most transition-specific motifs:")
print(head(specificity_table, 15))

message("\n========================================")
message("Refined figures saved:")
message("  38b_motif_heatmap_zscored")
message("  39b_motif_volcano_per_transition")
message("\nNew table:")
message("  14b_motif_specificity_ranked.csv")
message("========================================")