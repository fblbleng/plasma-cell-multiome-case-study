# ============================================================
# 10_hdWGCNA_setup.R   (v2)
#
# Sets up hdWGCNA on the integrated plasma cell Seurat object,
# builds metacells, and runs the soft-thresholding sweep to
# select the network construction power (beta).
#
# Changes from v1:
#   - Ig and TR genes excluded from the gene panel (consistent
#     with the variable-feature exclusion in Phase 1)
#   - MetacellsByGroups call uses `layer = "counts"` for
#     Seurat v5 compatibility (passes both `slot` and `layer`)
#
# Requires SeuratObject 5.0.2 for hdWGCNA compatibility.
# See: https://github.com/smorabit/hdWGCNA/issues/408
#
# Inputs:
#   - data/processed/seurat_integrated_annotated.rds
#   - data/processed/lncRNA_detected_symbols.rds
#
# Outputs:
#   - data/processed/seurat_hdWGCNA_setup.rds
#   - results/figures/23_softthreshold_diagnostic
#   - results/tables/09_softthreshold_diagnostic.csv
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(hdWGCNA)
  library(WGCNA)
  library(harmony)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(cowplot)
})

allowWGCNAThreads(nThreads = 4)
set.seed(42)

out_dir <- "data/processed"
fig_dir <- "results/figures"
tbl_dir <- "results/tables"

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
# 1. Load object + lncRNA panel
# ------------------------------------------------------------
seu <- readRDS(file.path(out_dir, "seurat_integrated_annotated.rds"))
seu$stage <- factor(seu$stage, levels = c("MBC", "prePB", "PB", "PC"))
DefaultAssay(seu) <- "RNA"
message("Seurat object: ", ncol(seu), " cells, ", nrow(seu), " genes")

detected_lncRNAs <- readRDS(file.path(out_dir, "lncRNA_detected_symbols.rds"))
message("Detected lncRNAs available: ", length(detected_lncRNAs))

# ------------------------------------------------------------
# 2. Build gene panel: variable features + lncRNAs, MINUS Ig and TR
# ------------------------------------------------------------
# We exclude Ig and TR genes for the same reason as Phase 1:
# they are so dominantly expressed in PB/PC that they would
# drive co-expression modules around antibody class identity
# rather than around regulatory programs.
# ------------------------------------------------------------
ig_genes <- grep("^IG[HKL][VDJC]", rownames(seu), value = TRUE)
ig_genes <- unique(c(ig_genes,
                     grep("^JCHAIN$|^JCH$", rownames(seu), value = TRUE)))
tr_genes <- grep("^TR[ABGD][VDJC]", rownames(seu), value = TRUE)
genes_exclude <- unique(c(ig_genes, tr_genes))
message("\nIg/TR genes to exclude: ", length(genes_exclude))

vf <- VariableFeatures(seu)
vf_clean <- setdiff(vf, genes_exclude)
lncRNA_clean <- setdiff(detected_lncRNAs, genes_exclude)
gene_panel <- intersect(union(vf_clean, lncRNA_clean), rownames(seu))

message("Gene panel composition:")
message("  Variable features (clean):  ", length(vf_clean),
        "  (", length(vf), " - ", length(intersect(vf, genes_exclude)),
        " Ig/TR)")
message("  Detected lncRNAs (clean):   ", length(lncRNA_clean))
message("  Final gene panel (union):   ", length(gene_panel))
message("    lncRNAs in panel:         ",
        length(intersect(gene_panel, detected_lncRNAs)))

# ------------------------------------------------------------
# 3. Initialize hdWGCNA
# ------------------------------------------------------------
message("\nSetting up hdWGCNA...")
seu <- SetupForWGCNA(
  seu,
  gene_select = "custom",
  features    = gene_panel,
  wgcna_name  = "plasma_cell"
)

# ------------------------------------------------------------
# 4. Build metacells (Seurat v5 compatible call)
# ------------------------------------------------------------
# hdWGCNA's MetacellsByGroups expects both `slot` and `layer`
# in some versions. We pass both to be safe and explicit.
#
# group.by = c("stage", "replicate") groups by both, so each
# metacell is from a single stage AND single replicate. This
# prevents metacells from blurring the differentiation
# trajectory or the batch structure.
# ------------------------------------------------------------
message("\nConstructing metacells (k=25, grouped by stage + replicate)...")
seu <- MetacellsByGroups(
  seu,
  group.by        = c("stage", "replicate"),
  ident.group     = "stage",
  reduction       = "harmony",
  k               = 25,
  max_shared      = 10,
  min_cells       = 50,
  target_metacells = 1000,
  max_iter        = 5000,
  mode            = "average",
  assay           = "RNA",
  slot            = "counts",
  layer           = "counts",
  verbose         = FALSE,
  wgcna_name      = "plasma_cell"
)

seu <- NormalizeMetacells(seu)

metacell_obj <- GetMetacellObject(seu)
message("\nMetacell counts per stage:")
print(table(metacell_obj$stage))
message("Total metacells: ", ncol(metacell_obj))

# ------------------------------------------------------------
# 5. Set the expression matrix for network building
# ------------------------------------------------------------
seu <- SetDatExpr(
  seu,
  group_name      = c("MBC", "prePB", "PB", "PC"),
  group.by        = "stage",
  assay           = "RNA",
  slot            = "data",
  layer           = "data"
)

# ------------------------------------------------------------
# 6. Soft-thresholding sweep
# ------------------------------------------------------------
message("\nTesting soft-thresholding powers (this takes ~3-5 minutes)...")
seu <- TestSoftPowers(
  seu,
  networkType = "signed"
)

power_table <- GetPowerTable(seu)
message("\nSoft-thresholding diagnostic table:")
print(power_table)

write.csv(power_table,
          file.path(tbl_dir, "09_softthreshold_diagnostic.csv"),
          row.names = FALSE)

# ------------------------------------------------------------
# 7. Diagnostic plot
# ------------------------------------------------------------
plot_list <- PlotSoftPowers(seu)

diag_plot <- wrap_plots(plot_list, ncol = 2) +
  plot_annotation(
    title    = "Soft-thresholding diagnostic for hdWGCNA",
    subtitle = "Pick the smallest beta where scale-free R^2 > 0.8 (red horizontal line)",
    theme = theme(
      plot.title    = element_text(face = "bold", size = 13, color = "#2a2622"),
      plot.subtitle = element_text(size = 10, color = "#5c544b")
    )
  )

ggsave(file.path(fig_dir, "23_softthreshold_diagnostic.pdf"),
       diag_plot, width = 12, height = 9)
ggsave(file.path(fig_dir, "23_softthreshold_diagnostic.png"),
       diag_plot, width = 12, height = 9, dpi = 200)
print(diag_plot)

# ------------------------------------------------------------
# 8. Recommend beta
# ------------------------------------------------------------
recommended_beta <- GetPowerTable(seu) %>%
  filter(SFT.R.sq > 0.8) %>%
  arrange(Power) %>%
  slice_head(n = 1) %>%
  pull(Power)

if (length(recommended_beta) == 0) {
  recommended_beta <- GetPowerTable(seu) %>%
    arrange(desc(SFT.R.sq)) %>%
    slice_head(n = 1) %>%
    pull(Power)
  message("\nNote: no beta reached R^2 > 0.8. Falling back to best available.")
}

message("\n========================================")
message("RECOMMENDED BETA: ", recommended_beta)
message("========================================")

# ------------------------------------------------------------
# 9. Save
# ------------------------------------------------------------
saveRDS(seu, file.path(out_dir, "seurat_hdWGCNA_setup.rds"))

message("\n========================================")
message("Script 10 (hdWGCNA setup) complete!")
message("\nGene panel composition:")
message("  Variable features (clean):  ", length(vf_clean))
message("  Detected lncRNAs (clean):   ", length(lncRNA_clean))
message("  Final gene panel:           ", length(gene_panel))
message("  lncRNAs in panel:           ",
        length(intersect(gene_panel, detected_lncRNAs)))
message("\nMetacells:")
print(table(metacell_obj$stage))
message("  Total: ", ncol(metacell_obj))
message("\nRecommended beta: ", recommended_beta)
message("\nFigure: 23_softthreshold_diagnostic")
message("Object saved: seurat_hdWGCNA_setup.rds")
message("Next: 11_hdWGCNA_modules.R (uses beta = ", recommended_beta, ")")
message("========================================")
