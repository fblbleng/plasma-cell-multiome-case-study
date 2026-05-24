# ============================================================
# 05a_cellcycle_diagnostic.R
# Cell cycle scoring diagnostic (run BEFORE script 05).
# Decide whether to regress out cell cycle by visualizing
# how strongly it confounds the biology.
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

set.seed(42)
out_dir <- "data/processed"
fig_dir <- "results/figures"
tbl_dir <- "results/tables"

stage_colors <- c("MBC"   = "#7d8c6e",
                  "prePB" = "#b88a3e",
                  "PB"    = "#b15835",
                  "PC"    = "#7c5c6b")

phase_colors <- c("G1" = "#7d8c6e", "S" = "#b88a3e", "G2M" = "#b15835")

theme_portfolio <- function() {
  theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 12, color = "#2a2622"),
      plot.subtitle = element_text(size = 10, color = "#5c544b"),
      axis.text = element_text(color = "#2a2622", size = 9),
      axis.title = element_text(color = "#5c544b", size = 10),
      panel.grid.minor = element_blank()
    )
}

# ------------------------------------------------------------
# 1. Load QC-filtered object
# ------------------------------------------------------------
seu <- readRDS(file.path(out_dir, "seurat_merged_qc_filtered.rds"))
seu$stage <- factor(seu$stage, levels = c("MBC", "prePB", "PB", "PC"))
seu$replicate <- as.factor(seu$replicate)
message("Loaded: ", ncol(seu), " cells, ", nrow(seu), " genes")

# ------------------------------------------------------------
# 2. LogNormalize and cell cycle scoring
# ------------------------------------------------------------
message("\nLogNormalize + variable features + scaling...")
seu <- NormalizeData(seu, verbose = FALSE)
seu <- FindVariableFeatures(seu, nfeatures = 3000, verbose = FALSE)

# Exclude Ig and TR genes from variable features
ig_genes <- grep("^IG[HKL][VDJC]", rownames(seu), value = TRUE)
ig_genes <- c(ig_genes, grep("^JCHAIN$|^JCH$", rownames(seu), value = TRUE))
tr_genes <- grep("^TR[ABGD][VDJC]", rownames(seu), value = TRUE)
genes_exclude <- unique(c(ig_genes, tr_genes))

vf <- VariableFeatures(seu)
VariableFeatures(seu) <- setdiff(vf, genes_exclude)
message("Variable features after Ig/TR exclusion: ",
        length(VariableFeatures(seu)))

# Score cell cycle (using Seurat's built-in Tirosh et al. lists)
s_genes   <- intersect(cc.genes.updated.2019$s.genes,   rownames(seu))
g2m_genes <- intersect(cc.genes.updated.2019$g2m.genes, rownames(seu))
message("S phase markers available:  ", length(s_genes),  "/43")
message("G2M phase markers available: ", length(g2m_genes), "/54")

seu <- CellCycleScoring(seu,
                        s.features   = s_genes,
                        g2m.features = g2m_genes,
                        set.ident = FALSE)
seu$Phase <- factor(seu$Phase, levels = c("G1", "S", "G2M"))

# ------------------------------------------------------------
# 3. Diagnostic plots
# ------------------------------------------------------------

# Plot 1: Phase composition per stage (the "do I need regression?" plot)
cc_df <- as.data.frame(prop.table(table(seu$stage, seu$Phase), margin = 1) * 100)
colnames(cc_df) <- c("Stage", "Phase", "Percent")
cc_df$Phase <- factor(cc_df$Phase, levels = c("G1", "S", "G2M"))

p1 <- ggplot(cc_df, aes(x = Stage, y = Percent, fill = Phase)) +
  geom_bar(stat = "identity", color = "white", linewidth = 0.3) +
  geom_text(aes(label = paste0(round(Percent), "%")),
            position = position_stack(vjust = 0.5),
            color = "white", size = 3, fontface = "bold") +
  scale_fill_manual(values = phase_colors) +
  theme_portfolio() +
  labs(title = "Cell cycle phase composition per stage",
       subtitle = "If prePB / PB are dominantly cycling, regression is recommended",
       x = NULL, y = "% of cells") +
  theme(legend.position = "right")

# Plot 2: S.Score violin per stage
p2 <- ggplot(seu@meta.data, aes(x = stage, y = S.Score, fill = stage)) +
  geom_violin(scale = "width", alpha = 0.85,
              color = "#2a2622", linewidth = 0.3) +
  geom_boxplot(width = 0.15, outlier.size = 0.3, outlier.alpha = 0.3,
               fill = "white", color = "#2a2622", linewidth = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "#8a7f73", linewidth = 0.4) +
  scale_fill_manual(values = stage_colors) +
  theme_portfolio() +
  labs(title = "S phase score per stage", x = NULL, y = "S.Score") +
  theme(legend.position = "none")

# Plot 3: G2M.Score violin per stage
p3 <- ggplot(seu@meta.data, aes(x = stage, y = G2M.Score, fill = stage)) +
  geom_violin(scale = "width", alpha = 0.85,
              color = "#2a2622", linewidth = 0.3) +
  geom_boxplot(width = 0.15, outlier.size = 0.3, outlier.alpha = 0.3,
               fill = "white", color = "#2a2622", linewidth = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "#8a7f73", linewidth = 0.4) +
  scale_fill_manual(values = stage_colors) +
  theme_portfolio() +
  labs(title = "G2/M phase score per stage", x = NULL, y = "G2M.Score") +
  theme(legend.position = "none")

# Plot 4: S.Score vs G2M.Score scatter (the classic Tirosh plot)
p4 <- ggplot(seu@meta.data, aes(x = S.Score, y = G2M.Score, color = stage)) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "#8a7f73", linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed",
             color = "#8a7f73", linewidth = 0.3) +
  geom_point(size = 0.4, alpha = 0.5) +
  scale_color_manual(values = stage_colors) +
  theme_portfolio() +
  labs(title = "Cell cycle position by stage",
       subtitle = "Top-right = actively cycling | Bottom-left = quiescent",
       x = "S.Score", y = "G2M.Score") +
  theme(legend.position = "right") +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))

# ------------------------------------------------------------
# 4. Pre-regression diagnostic UMAP
# ------------------------------------------------------------
# Compute a quick UMAP WITHOUT regressing cell cycle, just to see
# whether cells cluster primarily by Phase rather than by Stage.
# This is the most informative diagnostic.
# ------------------------------------------------------------
message("\nComputing diagnostic UMAP (no cell cycle regression)...")

# Scale data regressing only mt (NOT cell cycle), for comparison
seu <- ScaleData(seu, vars.to.regress = "percent.mt", verbose = FALSE)
seu <- RunPCA(seu, npcs = 30, verbose = FALSE)
seu <- RunUMAP(seu, dims = 1:30, reduction.name = "umap_diag", verbose = FALSE)

p5_stage <- DimPlot(seu, reduction = "umap_diag", group.by = "stage",
                    cols = stage_colors, pt.size = 0.4) +
  theme_portfolio() +
  labs(title = "Diagnostic UMAP - by stage",
       subtitle = "No cell cycle regression") +
  theme(legend.position = "right")

p5_phase <- DimPlot(seu, reduction = "umap_diag", group.by = "Phase",
                    cols = phase_colors, pt.size = 0.4) +
  theme_portfolio() +
  labs(title = "Diagnostic UMAP - by cell cycle phase",
       subtitle = "If Phase makes clear separate clusters here, regression is needed") +
  theme(legend.position = "right")

# ------------------------------------------------------------
# 5. Assemble the diagnostic dashboard
# ------------------------------------------------------------
top_row    <- p1 | p4
middle_row <- p2 | p3
bottom_row <- p5_stage | p5_phase

dashboard <- top_row / middle_row / bottom_row +
  plot_annotation(
    title = "Cell cycle diagnostic - do we need to regress?",
    subtitle = "Examine all 6 panels before deciding. Decision logic in comments below.",
    theme = theme(
      plot.title = element_text(face = "bold", size = 16, color = "#2a2622"),
      plot.subtitle = element_text(size = 11, color = "#5c544b")
    )
  )

ggsave(file.path(fig_dir, "04a_cellcycle_diagnostic.pdf"),
       dashboard, width = 13, height = 13)
ggsave(file.path(fig_dir, "04a_cellcycle_diagnostic.png"),
       dashboard, width = 13, height = 13, dpi = 200)
print(dashboard)

# ------------------------------------------------------------
# 6. Quantitative summary
# ------------------------------------------------------------
phase_table <- as.data.frame.matrix(table(seu$stage, seu$Phase))
phase_pct   <- as.data.frame.matrix(
  round(prop.table(table(seu$stage, seu$Phase), margin = 1) * 100, 1)
)

message("\n=== Cell cycle phase counts per stage ===")
print(phase_table)
message("\n=== Cell cycle phase percentages per stage ===")
print(phase_pct)

# Save phase metadata so we don't recompute
saveRDS(seu, file.path(out_dir, "seurat_with_cellcycle.rds"))

write.csv(phase_pct, file.path(tbl_dir, "03a_cellcycle_percentages.csv"))

# ------------------------------------------------------------
# DECISION LOGIC (for you to evaluate):
# ------------------------------------------------------------
#
# REGRESS cell cycle if:
#   - prePB and/or PB are >40% in S+G2M phase  (strong cycling)
#   - The diagnostic UMAP (p5_phase) shows visible clusters
#     forming by Phase rather than by Stage
#   - S.Score / G2M.Score violins separate cleanly by stage
#
# DON'T regress if:
#   - All stages have similar cycling distributions
#   - Phase is randomly distributed in the diagnostic UMAP
#   - You want to PRESERVE cell-cycle biology
#     (e.g., if "highly proliferating prePB" is a key finding)
#
# THIRD OPTION: regress G2M.Score - S.Score
#   - Removes phase identity (G1 vs S vs G2M) but preserves
#     overall proliferation signal
#   - Useful when proliferation IS part of the biology
#
# ------------------------------------------------------------

message("\n========================================")
message("Diagnostic complete.")
message("Open results/figures/04a_cellcycle_diagnostic.pdf and examine.")
message("Tell me what you see and we'll choose the right approach for script 05.")
message("========================================")
