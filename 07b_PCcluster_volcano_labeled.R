# ============================================================
# 07b_PCcluster_volcano_labeled.R
# Regenerates the C2 vs C4 PC cluster volcano with gene
# labels on top hits and biological annotation.
#
# Standalone: assumes the seurat_integrated_annotated.rds and
# the PC comparison table already exist from script 07.
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
})

out_dir <- "data/processed"
fig_dir <- "results/figures"
tbl_dir <- "results/tables"

theme_portfolio <- function() {
  theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 13, color = "#2a2622"),
      plot.subtitle = element_text(size = 10, color = "#5c544b"),
      axis.text = element_text(color = "#2a2622", size = 10),
      axis.title = element_text(color = "#5c544b", size = 11),
      panel.grid.minor = element_blank()
    )
}

# ------------------------------------------------------------
# 1. Load the saved DE result
# ------------------------------------------------------------
pc_de_path <- file.path(tbl_dir, "06_PCcluster_comparison.csv")
if (!file.exists(pc_de_path)) {
  stop("Run script 07 first - PC cluster comparison table not found")
}
pc_de <- read.csv(pc_de_path, stringsAsFactors = FALSE)
message("Loaded ", nrow(pc_de), " DE genes between PC clusters")

# ------------------------------------------------------------
# 2. Identify gene categories for coloring & labeling
# ------------------------------------------------------------
# Isotype constant-region genes (the dominant signal)
igg_genes <- grep("^IGHG", pc_de$gene, value = TRUE)
iga_genes <- grep("^IGHA", pc_de$gene, value = TRUE)
ighm      <- pc_de$gene[pc_de$gene == "IGHM"]
ig_other  <- c("IGJ", "IGHJ4", "JCHAIN", "JCH")
ig_other  <- intersect(ig_other, pc_de$gene)

# Interferon-stimulated genes (the C2 signature we found)
ifn_genes <- pc_de$gene[grepl(
  "^IFI|^IFIT|^ISG|^MX[12]$|^RSAD2$|^OAS[123L]$|^STAT1$|^IRF7$",
  pc_de$gene)]

# Cap p-values to avoid Inf in -log10
pc_de$p_plot <- pmax(pc_de$p_val_adj, 1e-300)
pc_de$neg_log10_p <- -log10(pc_de$p_plot)

# Classify
pc_de$category <- "Other"
pc_de$category[pc_de$gene %in% c(igg_genes)]          <- "IgG genes"
pc_de$category[pc_de$gene %in% c(iga_genes, ighm,
                                  ig_other)]            <- "IgA / IgM genes"
pc_de$category[pc_de$gene %in% ifn_genes]              <- "IFN-stimulated"

# Significance
pc_de$sig <- pc_de$p_val_adj < 0.01 & abs(pc_de$avg_log2FC) > 0.5

# ------------------------------------------------------------
# 3. Pick which genes to label
# ------------------------------------------------------------
# All isotype + IFN genes that pass threshold, plus top 5 in
# each direction not already in those categories.
to_label <- pc_de %>%
  filter(sig) %>%
  filter(category %in% c("IgG genes", "IgA / IgM genes", "IFN-stimulated") |
         (avg_log2FC > 0 & rank(-avg_log2FC) <= 5) |
         (avg_log2FC < 0 & rank( avg_log2FC) <= 5))

message("\nGenes to label: ", nrow(to_label))
print(to_label %>% select(gene, category, avg_log2FC, p_val_adj))

# ------------------------------------------------------------
# 4. Plot
# ------------------------------------------------------------
cat_colors <- c(
  "IgG genes"        = "#b15835",   # terracotta
  "IgA / IgM genes"  = "#7d8c6e",   # sage
  "IFN-stimulated"   = "#7c5c6b",   # plum
  "Other"            = "#cfc7b8"    # muted background
)

p <- ggplot(pc_de,
            aes(x = avg_log2FC, y = neg_log10_p, color = category)) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed",
             color = "#8a7f73", linewidth = 0.3) +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed",
             color = "#8a7f73", linewidth = 0.3) +
  geom_point(data = pc_de %>% filter(category == "Other"),
             alpha = 0.4, size = 1) +
  geom_point(data = pc_de %>% filter(category != "Other"),
             alpha = 0.9, size = 2) +
  ggrepel::geom_text_repel(
    data = to_label,
    aes(label = gene),
    size = 3.2,
    max.overlaps = 30,
    box.padding = 0.5,
    point.padding = 0.3,
    segment.color = "#8a7f73",
    segment.size = 0.3,
    color = "#2a2622",
    fontface = "italic",
    show.legend = FALSE
  ) +
  scale_color_manual(values = cat_colors,
                     name = NULL,
                     breaks = c("IgG genes", "IgA / IgM genes",
                                "IFN-stimulated", "Other")) +
  theme_portfolio() +
  labs(
    title = "PC cluster C2 vs C4: isotype switch + IFN signature",
    subtitle = "Right: enriched in C2 (IgG-PC, IFN-primed)  |  Left: enriched in C4 (IgA/IgM-PC)",
    x = "log2 fold change  (C2 vs C4)",
    y = "-log10(adjusted p-value)"
  ) +
  theme(legend.position = "top")

ggsave(file.path(fig_dir, "17_PCcluster_volcano.pdf"),
       p, width = 9, height = 7)
ggsave(file.path(fig_dir, "17_PCcluster_volcano.png"),
       p, width = 9, height = 7, dpi = 200)
print(p)

message("\nVolcano plot regenerated with labels and category colors.")
message("File: results/figures/17_PCcluster_volcano.pdf")
