# ============================================================
# 09b_lncRNA_featureplot_curated.R
#
# Replaces the auto-selected candidates in figure 22 with a
# curated panel of biologically meaningful lncRNAs, each
# marking a distinct stage of plasma cell differentiation.
#
# Curated selection (one informative lncRNA per stage):
#   CTA-250D10.23 -> MBC marker (avg_log2FC > 5 in DE)
#   MIR155HG      -> prePB-sub-cluster (cluster 5) marker
#   NEAT1         -> PB UPR activation marker
#   WT1-AS        -> PC marker, with specific cluster-2 (IgG-PC) bias
#
# These four together tell the full differentiation story
# at a glance.
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
})

set.seed(42)

out_dir <- "data/processed"
fig_dir <- "results/figures"

seu <- readRDS(file.path(out_dir, "seurat_integrated_annotated.rds"))
seu$stage <- factor(seu$stage, levels = c("MBC", "prePB", "PB", "PC"))
DefaultAssay(seu) <- "RNA"

# Curated panel
curated_lncRNAs <- c("CTA-250D10.23", "MIR155HG", "NEAT1", "WT1-AS")
curated_present <- intersect(curated_lncRNAs, rownames(seu))
message("Curated lncRNAs present in data: ",
        paste(curated_present, collapse = ", "))

theme_portfolio <- function() {
  theme_minimal(base_size = 11) +
    theme(
      plot.title    = element_text(face = "bold", size = 13, color = "#2a2622"),
      plot.subtitle = element_text(size = 9, color = "#5c544b", face = "italic"),
      panel.grid.minor = element_blank()
    )
}

# Pretty labels with biological context
labels <- list(
  "CTA-250D10.23" = list(title = "CTA-250D10.23",
                        subtitle = "MBC marker"),
  "MIR155HG"      = list(title = "MIR155HG",
                        subtitle = "B-cell activation lncRNA (host of miR-155)"),
  "NEAT1"         = list(title = "NEAT1",
                        subtitle = "PB / UPR activation marker"),
  "WT1-AS"        = list(title = "WT1-AS",
                        subtitle = "PC marker (IgG-PC bias)")
)

# Build individual plots
plots <- list()
for (gene in curated_present) {
  lab <- labels[[gene]]
  p <- FeaturePlot(seu, features = gene, reduction = "umap",
                   cols = c("#f3ede3", "#b15835"), pt.size = 0.3,
                   order = TRUE) +  # high-expressing cells on top
    theme_portfolio() +
    NoAxes() +
    labs(title = lab$title, subtitle = lab$subtitle)
  plots[[gene]] <- p
}

# 2x2 grid
if (length(plots) == 4) {
  combined <- (plots[[1]] | plots[[2]]) / (plots[[3]] | plots[[4]]) +
    plot_annotation(
      title = "Curated lncRNA markers across plasma cell differentiation",
      subtitle = "One stage-specific lncRNA per differentiation stage",
      theme = theme(
        plot.title = element_text(face = "bold", size = 15, color = "#2a2622"),
        plot.subtitle = element_text(size = 11, color = "#5c544b")
      )
    )
} else {
  combined <- wrap_plots(plots, ncol = 2) +
    plot_annotation(
      title = "Curated lncRNA markers across plasma cell differentiation",
      theme = theme(
        plot.title = element_text(face = "bold", size = 15, color = "#2a2622")
      )
    )
}

ggsave(file.path(fig_dir, "22_lncRNA_candidates_umap.pdf"),
       combined, width = 12, height = 9)
ggsave(file.path(fig_dir, "22_lncRNA_candidates_umap.png"),
       combined, width = 12, height = 9, dpi = 200)

print(combined)

message("\nFigure 22 regenerated with curated lncRNAs:")
message("  ", paste(curated_present, collapse = ", "))
message("File: results/figures/22_lncRNA_candidates_umap.png")
