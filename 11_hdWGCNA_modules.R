# ============================================================
# 11_hdWGCNA_modules.R   (v3 - with all network visualizations)
#
# Constructs the hdWGCNA co-expression network with beta = 4,
# identifies modules, computes module eigengenes, finds hub
# genes, and generates network visualizations:
#   - kME distribution per module (PlotKMEs)
#   - Per-module circular network plots (ModuleNetworkPlot)
#   - Unified force-directed network (HubGeneNetworkPlot)
#
# Inputs:
#   - data/processed/seurat_hdWGCNA_setup.rds
#
# Outputs:
#   - data/processed/seurat_hdWGCNA_modules.rds
#   - results/figures/24_dendrogram
#   - results/figures/25_module_sizes
#   - results/figures/26_module_eigengenes_umap
#   - results/figures/27_kME_per_module
#   - results/figures/28_hub_network_unified
#   - results/figures/ModuleNetworks/ (one PDF per module)
#   - results/tables/10_module_assignments.csv
#   - results/tables/10_hub_genes_per_module.csv
#   - results/tables/10_modules_lncRNA_content.csv
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(hdWGCNA)
  library(WGCNA)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(cowplot)
  library(igraph)
  library(ggraph)
  library(tidygraph)
})

allowWGCNAThreads(nThreads = 4)
set.seed(42)

# hdWGCNA's network plotting passes the TOM matrix to parallel
# workers via the future package. With ~5,500 genes the TOM
# exceeds the default 500 MiB limit. Bump it to 8 GB.
options(future.globals.maxSize = 8 * 1024^3)

out_dir <- "data/processed"
fig_dir <- "results/figures"
tbl_dir <- "results/tables"
modnet_dir <- file.path(fig_dir, "ModuleNetworks")
dir.create(modnet_dir, showWarnings = FALSE, recursive = TRUE)

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
# 1. Load
# ------------------------------------------------------------
seu <- readRDS(file.path(out_dir, "seurat_hdWGCNA_setup.rds"))
seu$stage <- factor(seu$stage, levels = c("MBC", "prePB", "PB", "PC"))
message("Loaded hdWGCNA setup object: ", ncol(seu), " cells")

# ------------------------------------------------------------
# 2. Construct the network with beta = 4
# ------------------------------------------------------------
message("\nConstructing network with beta = 4 (this takes ~3-5 minutes)...")
seu <- ConstructNetwork(
  seu,
  soft_power        = 4,
  setDatExpr        = FALSE,
  networkType       = "signed",
  TOMType           = "signed",
  detectCutHeight   = 0.995,
  minModuleSize     = 30,
  mergeCutHeight    = 0.25,
  numericLabels     = FALSE,
  pamRespectsDendro = FALSE,
  overwrite_tom     = TRUE,
  wgcna_name        = "plasma_cell"
)

# ------------------------------------------------------------
# 3. Module diagnostics
# ------------------------------------------------------------
modules <- GetModules(seu)
message("\nNetwork construction complete!")
message("Total genes in modules: ", nrow(modules))
message("Modules identified:")
print(table(modules$module))

module_sizes <- modules %>%
  filter(module != "grey") %>%
  count(module, sort = TRUE)
message("\nModule sizes (excluding grey/unassigned):")
print(module_sizes)

non_grey_modules <- setdiff(unique(modules$module), "grey")
n_modules <- length(non_grey_modules)

# ------------------------------------------------------------
# 4. Dendrogram
# ------------------------------------------------------------
pdf(file.path(fig_dir, "24_dendrogram.pdf"), width = 10, height = 5)
PlotDendrogram(seu, main = "Co-expression module dendrogram (beta = 4)")
dev.off()

png(file.path(fig_dir, "24_dendrogram.png"),
    width = 10, height = 5, units = "in", res = 200)
PlotDendrogram(seu, main = "Co-expression module dendrogram (beta = 4)")
dev.off()

# ------------------------------------------------------------
# 5. Module eigengenes
# ------------------------------------------------------------
message("\nComputing module eigengenes...")
seu <- ModuleEigengenes(
  seu,
  group.by.vars = "replicate",
  verbose       = FALSE
)

hMEs <- GetMEs(seu, harmonized = TRUE)
message("Eigengenes computed for ", ncol(hMEs), " modules")

# ------------------------------------------------------------
# 6. Module connectivity (kME) and hub genes
# ------------------------------------------------------------
message("\nComputing module connectivity (kME)...")
seu <- ModuleConnectivity(
  seu,
  group.by   = "stage",
  group_name = c("MBC", "prePB", "PB", "PC")
)

hub_genes_top10 <- GetHubGenes(seu, n_hubs = 10)
hub_genes_top25 <- GetHubGenes(seu, n_hubs = 25)

write.csv(hub_genes_top25,
          file.path(tbl_dir, "10_hub_genes_per_module.csv"),
          row.names = FALSE)

message("\nTop 5 hub genes per module:")
print(
  hub_genes_top10 %>%
    group_by(module) %>%
    slice_head(n = 5) %>%
    select(module, gene_name, kME) %>%
    arrange(module)
)

# ------------------------------------------------------------
# 7. PlotKMEs: kME distribution per module
# ------------------------------------------------------------
message("\nGenerating PlotKMEs figure...")

p_kme <- PlotKMEs(seu, ncol = 4, n_hubs = 10)

ncol_kme <- 4
nrow_kme <- ceiling(n_modules / ncol_kme)

ggsave(file.path(fig_dir, "27_kME_per_module.pdf"),
       p_kme, width = ncol_kme * 3.5, height = nrow_kme * 2.8,
       limitsize = FALSE)
ggsave(file.path(fig_dir, "27_kME_per_module.png"),
       p_kme, width = ncol_kme * 3.5, height = nrow_kme * 2.8,
       dpi = 150, limitsize = FALSE)

# ------------------------------------------------------------
# 8. ModuleNetworkPlot: one circular network per module
# ------------------------------------------------------------
# Saves one PDF per module to results/figures/ModuleNetworks/.
# Each PDF shows the top 10 hub genes in the inner circle
# and the next 15 in the outer circle, connected by TOM edges.
# ------------------------------------------------------------
message("\nGenerating per-module network plots (one PDF per module)...")
ModuleNetworkPlot(
  seu,
  outdir       = modnet_dir,
  n_inner      = 10,
  n_outer      = 15,
  n_conns      = Inf,
  mods         = "all",
  plot_size    = c(7, 7),
  edge.alpha   = 0.25,
  edge.width   = 1,
  vertex.label.cex = 1,
  vertex.size  = 6,
  wgcna_name   = "plasma_cell"
)
message("  Saved to: ", modnet_dir)

# ------------------------------------------------------------
# 9. HubGeneNetworkPlot: unified network of all modules
# ------------------------------------------------------------
# One figure showing the top hub genes from all modules in
# a single force-directed layout. Good summary figure for the
# case study writeup.
# ------------------------------------------------------------
message("\nGenerating unified hub-gene network plot...")

p_hubnet <- HubGeneNetworkPlot(
  seu,
  mods            = "all",
  n_hubs          = 6,
  n_other         = 3,
  sample_edges    = TRUE,
  edge_prop       = 0.5,
  edge.alpha      = 0.25,
  vertex.label.cex = 0.5,
  hub.vertex.size  = 4,
  other.vertex.size = 1,
  wgcna_name      = "plasma_cell",
  return_graph    = FALSE
)

# HubGeneNetworkPlot prints the plot directly; we capture & save
pdf(file.path(fig_dir, "28c_hubnet_clean.pdf"),
    width = 10, height = 10)
HubGeneNetworkPlot(
  seu,
  mods             = "all",
  n_hubs           = 8,
  n_other          = 2,
  sample_edges     = TRUE,
  edge_prop        = 0.75,
  edge.alpha       = 0.25,
  vertex.label.cex = 1,
  hub.vertex.size  = 2,
  other.vertex.size = 1,
  wgcna_name       = "plasma_cell"
)
dev.off()

png(file.path(fig_dir, "28c_hubnet_clean.png"),
    width = 10, height = 10, units = "in", res = 200)
HubGeneNetworkPlot(
  seu,
  mods             = "all",
  n_hubs           = 8,
  n_other          = 2,
  sample_edges     = TRUE,
  edge_prop        = 0.75,
  edge.alpha       = 0.25,
  vertex.label.cex = 1,
  hub.vertex.size  = 2,
  other.vertex.size = 1,
  wgcna_name       = "plasma_cell"
)
dev.off()

# ------------------------------------------------------------
# 10. Save module assignments
# ------------------------------------------------------------
modules_export <- modules %>%
  select(gene_name, module, color) %>%
  arrange(module, gene_name)
write.csv(modules_export,
          file.path(tbl_dir, "10_module_assignments.csv"),
          row.names = FALSE)

# lncRNA content per module
detected_lncRNAs <- readRDS(file.path(out_dir, "lncRNA_detected_symbols.rds"))
modules_with_lnc <- modules %>%
  mutate(is_lncRNA = gene_name %in% detected_lncRNAs) %>%
  group_by(module) %>%
  summarize(
    n_genes      = n(),
    n_lncRNAs    = sum(is_lncRNA),
    pct_lncRNA   = round(100 * n_lncRNAs / n_genes, 1),
    lncRNA_examples = paste(head(gene_name[is_lncRNA], 5), collapse = ", ")
  ) %>%
  arrange(desc(n_lncRNAs))

write.csv(modules_with_lnc,
          file.path(tbl_dir, "10_modules_lncRNA_content.csv"),
          row.names = FALSE)

message("\nlncRNA content per module:")
print(modules_with_lnc)

# ------------------------------------------------------------
# 11. Module sizes figure
# ------------------------------------------------------------
p_sizes <- module_sizes %>%
  ggplot(aes(x = reorder(module, n), y = n, fill = module)) +
  geom_col(width = 0.7) +
  scale_fill_identity() +
  coord_flip() +
  geom_text(aes(label = n), hjust = -0.15, size = 3.2, color = "#2a2622") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  theme_portfolio() +
  labs(
    title    = "Co-expression module sizes",
    subtitle = paste0(nrow(module_sizes), " modules identified, excluding 'grey' (unassigned)"),
    x = NULL,
    y = "Number of genes"
  ) +
  theme(legend.position = "none")

ggsave(file.path(fig_dir, "25_module_sizes.pdf"), p_sizes,
       width = 8, height = max(4, nrow(module_sizes) * 0.3))
ggsave(file.path(fig_dir, "25_module_sizes.png"), p_sizes,
       width = 8, height = max(4, nrow(module_sizes) * 0.3), dpi = 200)

# ------------------------------------------------------------
# 12. Module eigengene activity on UMAP
# ------------------------------------------------------------
message("\nGenerating module eigengene UMAP overlays...")

p_meumap <- ModuleFeaturePlot(
  seu,
  features = "hMEs",
  order    = "shuffle",
  raster   = FALSE
)

ncol_use2 <- ifelse(n_modules <= 4, 2,
             ifelse(n_modules <= 12, 3, 4))

p_meumap_combined <- wrap_plots(p_meumap, ncol = ncol_use2) +
  plot_annotation(
    title    = "Module eigengene activity across the differentiation trajectory",
    subtitle = paste0("Each panel shows one module's expression pattern (",
                      n_modules, " modules total, harmonized for batch)"),
    theme = theme(
      plot.title    = element_text(face = "bold", size = 14, color = "#2a2622"),
      plot.subtitle = element_text(size = 10, color = "#5c544b")
    )
  )

height_use <- ceiling(n_modules / ncol_use2) * 3.5
ggsave(file.path(fig_dir, "26_module_eigengenes_umap.pdf"),
       p_meumap_combined, width = ncol_use2 * 4, height = height_use,
       limitsize = FALSE)
ggsave(file.path(fig_dir, "26_module_eigengenes_umap.png"),
       p_meumap_combined, width = ncol_use2 * 4, height = height_use,
       dpi = 150, limitsize = FALSE)

# ------------------------------------------------------------
# 13. Save
# ------------------------------------------------------------
saveRDS(seu, file.path(out_dir, "seurat_hdWGCNA_modules.rds"))

message("\n========================================")
message("Script 11 (hdWGCNA modules) complete!")
message("\nNumbers:")
message("  Modules identified:    ", n_modules)
message("  Total genes assigned:  ", sum(modules$module != "grey"))
message("  Genes in 'grey':       ", sum(modules$module == "grey"))
message("\nlncRNA distribution:")
print(modules_with_lnc %>% select(module, n_genes, n_lncRNAs, pct_lncRNA))
message("\nFigures:")
message("  24_dendrogram               (clustering tree)")
message("  25_module_sizes")
message("  26_module_eigengenes_umap   (each module on UMAP)")
message("  27_kME_per_module           (kME distribution + hub gene labels)")
message("  28_hub_network_unified      (force-directed network across modules)")
message("  ModuleNetworks/<color>.pdf  (one circular network per module)")
message("\nTables:")
message("  10_module_assignments.csv")
message("  10_hub_genes_per_module.csv")
message("  10_modules_lncRNA_content.csv")
message("\nNext: 12_hdWGCNA_analysis.R (stage-level module activity, lncRNA integration)")
message("========================================")
