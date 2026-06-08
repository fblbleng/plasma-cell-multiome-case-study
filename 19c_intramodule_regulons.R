# ============================================================
# 19c_intramodule_regulons.R
#
# For each hdWGCNA module, extract the regulatory subnetwork
# where BOTH the TF and its target gene are members of that
# module. This gives the intra-module regulatory circuit for
# each program (yellow / brown / turquoise / blue), removing
# the module-size bias that dominated the global TF x module
# heatmap.
#
# Run after AssignTFRegulons, with `seu` and `regulons` in
# memory.
#
# Note on column names (this dataset's regulons table):
#   - TF column      : 'tf'
#   - target column  : 'gene'
#   - importance     : 'Gain'
#
# Outputs:
#   - results/tables/17_intramodule_regulon_edges.csv
#   - results/tables/17_intramodule_top_regulators.csv
#   - results/figures/49_intramodule_top_regulators.{pdf,png}
#   - results/figures/50_intramodule_networks.{pdf,png}
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(hdWGCNA)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(igraph)
  library(ggraph)
})

if (requireNamespace("conflicted", quietly = TRUE)) {
  conflicted::conflicts_prefer(dplyr::filter, .quiet = TRUE)
  conflicted::conflicts_prefer(dplyr::select, .quiet = TRUE)
  conflicted::conflicts_prefer(dplyr::mutate, .quiet = TRUE)
  conflicted::conflicts_prefer(dplyr::desc, .quiet = TRUE)
  conflicted::conflicts_prefer(dplyr::slice, .quiet = TRUE)
  conflicted::conflicts_prefer(base::intersect, .quiet = TRUE)
}

fig_dir <- "results/figures"
tbl_dir <- "results/tables"
out_dir <- "data/processed"

# If not in memory:
# seu <- readRDS(file.path(out_dir, "seurat_hdwgcna_tf_network.rds"))
# regulons <- GetTFRegulons(seu)

module_colors <- c(
  "yellow"    = "#e8b800",
  "brown"     = "#964b00",
  "turquoise" = "#4cb5b0",
  "blue"      = "#4472c4"
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
# 1. Module assignment lookup
# ------------------------------------------------------------
modules <- GetModules(seu)
gene_to_module <- setNames(as.character(modules$module), modules$gene_name)

# Detect column names in the regulons table
tf_col   <- intersect(c("tf", "TF"), colnames(regulons))[1]
gene_col <- intersect(c("gene", "gene_name", "target"), colnames(regulons))[1]
imp_col  <- intersect(c("Gain", "score", "importance"), colnames(regulons))[1]
message("Regulon columns -> TF: ", tf_col,
        " | target: ", gene_col, " | importance: ", imp_col)

# ------------------------------------------------------------
# 2. Annotate each regulon edge with TF module and target module
# ------------------------------------------------------------
reg <- regulons |>
  dplyr::rename(tf = !!tf_col, gene = !!gene_col, Gain = !!imp_col) |>
  dplyr::mutate(
    tf_module     = gene_to_module[tf],
    target_module = gene_to_module[gene]
  )

# ------------------------------------------------------------
# 3. Keep only intra-module edges: TF and target in SAME module
# ------------------------------------------------------------
intramodule <- reg |>
  dplyr::filter(
    !is.na(tf_module),
    !is.na(target_module),
    tf_module == target_module,
    tf_module != "grey",
    tf != gene
  ) |>
  dplyr::rename(module = tf_module) |>
  dplyr::select(tf, gene, module, Gain)

message("\nIntra-module regulon edges (TF + target in same module):")
print(table(intramodule$module))

write.csv(intramodule,
          file.path(tbl_dir, "17_intramodule_regulon_edges.csv"),
          row.names = FALSE)

# ------------------------------------------------------------
# 4. Top regulators within each module
# ------------------------------------------------------------
# Rank module-member TFs by how many module-member targets they
# regulate (intra-module out-degree).
# ------------------------------------------------------------
top_regulators <- intramodule |>
  dplyr::group_by(module, tf) |>
  dplyr::summarize(
    n_intramodule_targets = dplyr::n(),
    mean_gain = mean(Gain, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::arrange(module, dplyr::desc(n_intramodule_targets))

write.csv(top_regulators,
          file.path(tbl_dir, "17_intramodule_top_regulators.csv"),
          row.names = FALSE)

message("\nTop intra-module regulators per module:")
for (mod in c("yellow", "brown", "turquoise", "blue")) {
  cat("\n---", mod, "module ---\n")
  mod_top <- top_regulators |>
    dplyr::filter(module == mod) |>
    head(10)
  print(as.data.frame(mod_top[, c("tf", "n_intramodule_targets", "mean_gain")]))
}

# ------------------------------------------------------------
# 5. Figure 49: top regulators per module (faceted bar chart)
# ------------------------------------------------------------
message("\nFigure 49: top intra-module regulators...")

top_n_per_module <- top_regulators |>
  dplyr::group_by(module) |>
  dplyr::slice_head(n = 10) |>
  dplyr::ungroup() |>
  dplyr::mutate(
    module = factor(module,
                    levels = c("yellow", "brown", "turquoise", "blue"))
  )

p49 <- ggplot(top_n_per_module,
              aes(x = reorder(tf, n_intramodule_targets),
                  y = n_intramodule_targets, fill = module)) +
  geom_col(width = 0.72) +
  coord_flip() +
  facet_wrap(~ module, scales = "free", ncol = 2) +
  scale_fill_manual(values = module_colors, guide = "none") +
  theme_portfolio() +
  theme(strip.text = element_text(face = "bold", size = 11, color = "#2a2622")) +
  labs(
    title    = "Top transcription factors within each co-expression module",
    subtitle = "Intra-module regulons: TF and target both module members. Bar = number of module-member targets.",
    x = NULL, y = "Module-member targets"
  )

ggsave(file.path(fig_dir, "49_intramodule_top_regulators.pdf"),
       p49, width = 11, height = 9)
ggsave(file.path(fig_dir, "49_intramodule_top_regulators.png"),
       p49, width = 11, height = 9, dpi = 200)

# ------------------------------------------------------------
# 6. Figure 50: intra-module regulatory networks (one per module)
# ------------------------------------------------------------
# For each module, draw the regulatory subnetwork among the top
# regulators and their module-member targets.
# ------------------------------------------------------------
message("Figure 50: intra-module regulatory networks...")

make_module_network <- function(mod, max_tfs = 6, max_targets_per_tf = 8) {
  mod_edges <- intramodule |>
    dplyr::filter(module == mod)
  if (nrow(mod_edges) == 0) return(NULL)
  
  # Top TFs in this module
  top_mod_tfs <- top_regulators |>
    dplyr::filter(module == mod) |>
    dplyr::slice_head(n = max_tfs) |>
    dplyr::pull(tf)
  
  # For each top TF, keep its strongest targets
  sub_edges <- mod_edges |>
    dplyr::filter(tf %in% top_mod_tfs) |>
    dplyr::group_by(tf) |>
    dplyr::slice_max(Gain, n = max_targets_per_tf, with_ties = FALSE) |>
    dplyr::ungroup()
  
  if (nrow(sub_edges) == 0) return(NULL)
  
  g <- igraph::graph_from_data_frame(
    sub_edges[, c("tf", "gene", "Gain")],
    directed = TRUE
  )
  
  # Node type: is it one of the top TFs?
  V(g)$is_tf <- V(g)$name %in% top_mod_tfs
  
  ggraph(g, layout = "fr") +
    geom_edge_link(aes(width = Gain),
                   color = module_colors[mod], alpha = 0.4,
                   arrow = arrow(length = unit(2, "mm"), type = "closed"),
                   end_cap = circle(3, "mm")) +
    geom_node_point(aes(size = ifelse(is_tf, 6, 3),
                        shape = ifelse(is_tf, "TF", "target")),
                    color = module_colors[mod]) +
    geom_node_text(aes(label = name,
                       fontface = ifelse(is_tf, "bold", "plain")),
                   size = 3, repel = TRUE, color = "#2a2622") +
    scale_edge_width(range = c(0.3, 2), guide = "none") +
    scale_size_identity() +
    scale_shape_manual(values = c("TF" = 15, "target" = 16), guide = "none") +
    theme_void() +
    labs(title = paste0(mod, " module")) +
    theme(plot.title = element_text(face = "bold", size = 12,
                                    color = module_colors[mod], hjust = 0.5))
}

net_plots <- list()
for (mod in c("yellow", "brown", "turquoise", "blue")) {
  np <- tryCatch(make_module_network(mod),
                 error = function(e) {
                   message("  ", mod, " network skipped: ", e$message)
                   NULL
                 })
  if (!is.null(np)) net_plots[[mod]] <- np
}

if (length(net_plots) > 0) {
  combined_nets <- wrap_plots(net_plots, ncol = 2) +
    plot_annotation(
      title = "Intra-module regulatory circuits",
      subtitle = "Top module-member TFs (squares) and their module-member targets (circles)",
      theme = theme(
        plot.title = element_text(face = "bold", size = 14, color = "#2a2622"),
        plot.subtitle = element_text(size = 11, color = "#5c544b")
      )
    )
  n_rows <- ceiling(length(net_plots) / 2)
  ggsave(file.path(fig_dir, "50_intramodule_networks.pdf"),
         combined_nets, width = 12, height = 6 * n_rows)
  ggsave(file.path(fig_dir, "50_intramodule_networks.png"),
         combined_nets, width = 12, height = 6 * n_rows, dpi = 200)
}

message("\n========================================")
message("Intra-module regulon analysis complete!")
message("\nTables:")
message("  17_intramodule_regulon_edges.csv")
message("  17_intramodule_top_regulators.csv")
message("\nFigures:")
message("  49_intramodule_top_regulators (faceted bars per module)")
message("  50_intramodule_networks       (regulatory circuit per module)")
message("========================================")