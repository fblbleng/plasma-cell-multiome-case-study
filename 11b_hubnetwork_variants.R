# ============================================================
# 11b_hubnetwork_variants.R
#
# Regenerates HubGeneNetworkPlot with multiple parameter sets
# so we can choose the cleanest variant for the case study.
#
# Variant a: more hub genes per module, no edge sampling
# Variant b: hubs only, no "other" genes
# Variant c: balanced, sampled edges
#
# Loads the saved hdWGCNA modules object; no re-running of
# the slow network construction needed.
#
# Inputs:
#   - data/processed/seurat_hdWGCNA_modules.rds
#
# Outputs:
#   - results/figures/28a_hubnet_dense.{pdf,png}
#   - results/figures/28b_hubnet_hubs_only.{pdf,png}
#   - results/figures/28c_hubnet_clean.{pdf,png}
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(hdWGCNA)
  library(WGCNA)
  library(igraph)
  library(ggraph)
  library(tidygraph)
})

# Memory limit for parallel TOM transfer (required by HubGeneNetworkPlot)
options(future.globals.maxSize = 8 * 1024^3)

set.seed(42)
out_dir <- "data/processed"
fig_dir <- "results/figures"

# ------------------------------------------------------------
# 1. Load the modules object from script 11
# ------------------------------------------------------------
seu <- readRDS(file.path(out_dir, "seurat_hdWGCNA_modules.rds"))
message("Loaded: ", ncol(seu), " cells, ",
        length(unique(GetModules(seu)$module)),
        " modules (including grey)")

# ------------------------------------------------------------
# Variant a: more hubs, no edge sampling (denser)
# ------------------------------------------------------------
message("\nVariant a: 10 hubs + 5 other, all edges...")
pdf(file.path(fig_dir, "28a_hubnet_dense.pdf"), width = 12, height = 12)
HubGeneNetworkPlot(
  seu,
  mods             = "all",
  n_hubs           = 10,
  n_other          = 5,
  sample_edges     = FALSE,
  edge.alpha       = 0.15,
  vertex.label.cex = 0.6,
  hub.vertex.size  = 5,
  other.vertex.size = 1.5,
  wgcna_name       = "plasma_cell"
)
dev.off()

png(file.path(fig_dir, "28a_hubnet_dense.png"),
    width = 12, height = 12, units = "in", res = 200)
HubGeneNetworkPlot(
  seu, mods = "all", n_hubs = 10, n_other = 5,
  sample_edges = FALSE, edge.alpha = 0.15,
  vertex.label.cex = 0.6, hub.vertex.size = 5,
  other.vertex.size = 1.5, wgcna_name = "plasma_cell"
)
dev.off()

# ------------------------------------------------------------
# Variant b: hubs only (cleanest)
# ------------------------------------------------------------
message("\nVariant b: hubs only (no 'other' nodes)...")
pdf(file.path(fig_dir, "28b_hubnet_hubs_only.pdf"), width = 10, height = 10)
HubGeneNetworkPlot(
  seu,
  mods             = "all",
  n_hubs           = 12,
  n_other          = 0,
  sample_edges     = FALSE,
  edge.alpha       = 0.2,
  vertex.label.cex = 0.7,
  hub.vertex.size  = 6,
  other.vertex.size = 1,
  wgcna_name       = "plasma_cell"
)
dev.off()

png(file.path(fig_dir, "28b_hubnet_hubs_only.png"),
    width = 10, height = 10, units = "in", res = 200)
HubGeneNetworkPlot(
  seu, mods = "all", n_hubs = 12, n_other = 0,
  sample_edges = FALSE, edge.alpha = 0.2,
  vertex.label.cex = 0.7, hub.vertex.size = 6,
  other.vertex.size = 1, wgcna_name = "plasma_cell"
)
dev.off()

# ------------------------------------------------------------
# Variant c: balanced (chosen for the case study)
# ------------------------------------------------------------
message("\nVariant c: 8 hubs + 2 other, sampled edges...")
pdf(file.path(fig_dir, "28c_hubnet_clean.pdf"), width = 9, height = 9)
HubGeneNetworkPlot(
  seu,
  mods             = "all",
  n_hubs           = 8,
  n_other          = 2,
  sample_edges     = TRUE,
  edge_prop        = 0.75,
  edge.alpha       = 0.25,
  vertex.label.cex = 0.65,
  hub.vertex.size  = 5,
  other.vertex.size = 1.5,
  wgcna_name       = "plasma_cell"
)
dev.off()

png(file.path(fig_dir, "28c_hubnet_clean.png"),
    width = 9, height = 9, units = "in", res = 200)
HubGeneNetworkPlot(
  seu, mods = "all", n_hubs = 8, n_other = 2,
  sample_edges = TRUE, edge_prop = 0.75, edge.alpha = 0.25,
  vertex.label.cex = 0.65, hub.vertex.size = 5,
  other.vertex.size = 1.5, wgcna_name = "plasma_cell"
)
dev.off()

message("\n========================================")
message("Three variants generated:")
message("  28a_hubnet_dense       (10 hubs + 5 other, all edges)")
message("  28b_hubnet_hubs_only   (12 hubs each, no 'other', all edges)")
message("  28c_hubnet_clean       (8 hubs + 2 other, 75% edges, CHOSEN)")
message("========================================")
