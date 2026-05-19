# ============================================================
# 04_qc_and_filter.R   (v3 - protein-coding mt + stage-aware filter)
# Quality control for BD Rhapsody scRNA-seq data.
#
# Design decisions (documented for the case study narrative):
#
# 1. %mt computed on the 13 mitochondrial protein-coding genes
#    only (MT-ATP*, MT-CO*, MT-CYB, MT-ND*). The 2 mt-rRNAs
#    (MT-RNR1, MT-RNR2) inflate %mt by ~20-30% in this dataset
#    and aren't informative about cell health.
#
# 2. Stage-aware percentile thresholds (not fixed cutoffs).
#    Plasmablasts naturally have elevated %mt due to high
#    metabolic activity and FACS-sorting stress. A fixed 15%
#    cutoff removes 60%+ of PB cells — that's bias, not QC.
#    We use the 95th percentile *within each stage* instead.
#
# 3. Biological sanity check via %immunoglobulin: should rise
#    monotonically MBC -> prePB -> PB -> PC. This validates
#    the upstream sample tag demultiplexing.
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

stage_colors <- c("MBC"   = "#7d8c6e",   # sage
                  "prePB" = "#b88a3e",   # ochre
                  "PB"    = "#b15835",   # terracotta
                  "PC"    = "#7c5c6b")   # muted plum

# ------------------------------------------------------------
# 1. Load raw Seurat objects
# ------------------------------------------------------------
seu_rep1 <- readRDS(file.path(out_dir, "seurat_rep1_raw.rds"))
seu_rep2 <- readRDS(file.path(out_dir, "seurat_rep2_raw.rds"))

seu <- merge(seu_rep1, seu_rep2, project = "PC_diff")
message("Merged object:")
print(seu)

# ------------------------------------------------------------
# 2. Identify mt protein-coding genes explicitly
# ------------------------------------------------------------
# Pattern matches MT-ATPn, MT-COn, MT-CYB, MT-NDn (the 13
# protein-coding genes). Excludes MT-RNR* (rRNAs) and MT-T*
# (tRNAs), which inflate %mt without being health-informative.
mt_pc_genes <- grep("^MT-(ATP|CO|CYB|ND)", rownames(seu), value = TRUE)
message("\nMitochondrial protein-coding genes used for %mt: ",
        length(mt_pc_genes))
print(mt_pc_genes)

if (length(mt_pc_genes) < 13) {
  warning("Found ", length(mt_pc_genes), " mt protein-coding genes, expected 13.")
}

# ------------------------------------------------------------
# 3. Compute QC metrics
# ------------------------------------------------------------
seu[["percent.mt"]]   <- PercentageFeatureSet(seu, features = mt_pc_genes)
seu[["percent.ribo"]] <- PercentageFeatureSet(seu, pattern  = "^RP[SL]")
seu[["percent.ig"]]   <- PercentageFeatureSet(seu, pattern  = "^IG[HKL][VDJC]")

seu$stage <- factor(seu$stage, levels = c("MBC", "prePB", "PB", "PC"))

qc_summary <- seu@meta.data %>%
  group_by(stage) %>%
  summarize(
    n_cells = n(),
    median_nFeature = median(nFeature_RNA),
    median_nCount   = median(nCount_RNA),
    median_pct_mt   = round(median(percent.mt), 2),
    median_pct_ribo = round(median(percent.ribo), 2),
    median_pct_ig   = round(median(percent.ig), 2),
    .groups = "drop"
  )
message("\nQC summary per stage (mt = protein-coding only):")
print(qc_summary)
write.csv(qc_summary, file.path(tbl_dir, "02_qc_summary_per_stage.csv"),
          row.names = FALSE)

message("\nBiological sanity check (%Ig should rise MBC -> PC):")
message("  ", paste(qc_summary$median_pct_ig, collapse = " -> "))

# ------------------------------------------------------------
# 4. Theme + plotting helpers
# ------------------------------------------------------------
theme_portfolio <- function() {
  theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 12, color = "#2a2622"),
      plot.subtitle = element_text(size = 10, color = "#5c544b"),
      axis.text = element_text(color = "#2a2622"),
      axis.title = element_text(color = "#5c544b", size = 10),
      legend.position = "none",
      panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold", color = "#2a2622")
    )
}

make_violin <- function(data, metric, ylab, log_y = FALSE) {
  p <- ggplot(data, aes(x = stage, y = .data[[metric]], fill = stage)) +
    geom_violin(scale = "width", trim = TRUE, alpha = 0.85,
                color = "#2a2622", linewidth = 0.3) +
    geom_boxplot(width = 0.15, outlier.size = 0.3, outlier.alpha = 0.3,
                 fill = "white", color = "#2a2622", linewidth = 0.3) +
    scale_fill_manual(values = stage_colors) +
    labs(x = NULL, y = ylab) +
    theme_portfolio()
  if (log_y) p <- p + scale_y_log10(labels = comma)
  p
}

# ------------------------------------------------------------
# 5. Pre-filter QC dashboard
# ------------------------------------------------------------
md_pre <- seu@meta.data

p1 <- make_violin(md_pre, "nFeature_RNA", "Genes per cell")        + ggtitle("Gene complexity")
p2 <- make_violin(md_pre, "nCount_RNA",   "UMIs per cell", TRUE)   + ggtitle("UMI counts (log10)")
p3 <- make_violin(md_pre, "percent.mt",   "% mt (protein-coding)") + ggtitle("Mitochondrial content")
p4 <- make_violin(md_pre, "percent.ig",   "% immunoglobulin")      + ggtitle("Immunoglobulin content")

p5 <- ggplot(md_pre, aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt)) +
  geom_point(size = 0.4, alpha = 0.5) +
  facet_wrap(~ stage, nrow = 1) +
  scale_x_log10(labels = comma) +
  scale_color_gradientn(colors = c("#7d8c6e", "#b88a3e", "#b15835"),
                        name = "% mt") +
  labs(x = "UMIs per cell (log10)", y = "Genes per cell",
       title = "Gene complexity vs UMI depth, by stage") +
  theme_portfolio() +
  theme(legend.position = "right",
        strip.background = element_rect(fill = "#f3ede3", color = NA))

qc_dashboard_pre <- (p1 | p2 | p3 | p4) / p5 +
  plot_layout(heights = c(1, 1.2)) +
  plot_annotation(
    title = "scRNA-seq quality control dashboard",
    subtitle = "GSE242330 - Moreaux et al., Blood 2024 | Pre-filtering",
    theme = theme(
      plot.title = element_text(face = "bold", size = 15, color = "#2a2622"),
      plot.subtitle = element_text(size = 11, color = "#5c544b")
    )
  )

ggsave(file.path(fig_dir, "02_qc_dashboard_prefilter.pdf"),
       qc_dashboard_pre, width = 13, height = 8)
ggsave(file.path(fig_dir, "02_qc_dashboard_prefilter.png"),
       qc_dashboard_pre, width = 13, height = 8, dpi = 200)
print(qc_dashboard_pre)

# ------------------------------------------------------------
# 6. Stage-aware thresholds
# ------------------------------------------------------------
thresholds <- seu@meta.data %>%
  group_by(stage) %>%
  summarize(
    nFeature_lo = quantile(nFeature_RNA, 0.01),
    nFeature_hi = quantile(nFeature_RNA, 0.99),
    nCount_lo   = quantile(nCount_RNA,   0.01),
    nCount_hi   = quantile(nCount_RNA,   0.99),
    mt_hi       = quantile(percent.mt,   0.95),
    .groups = "drop"
  )
message("\nStage-aware thresholds (1st/99th for nFeat/nCount, 95th for mt):")
print(thresholds)
write.csv(thresholds, file.path(tbl_dir, "02_qc_thresholds.csv"),
          row.names = FALSE)

# ------------------------------------------------------------
# 7. Apply filters
# ------------------------------------------------------------
md <- seu@meta.data
md$Cell <- rownames(md)
md <- md %>%
  left_join(thresholds, by = "stage") %>%
  mutate(
    pass_nFeature = nFeature_RNA >= nFeature_lo & nFeature_RNA <= nFeature_hi,
    pass_nCount   = nCount_RNA   >= nCount_lo   & nCount_RNA   <= nCount_hi,
    pass_mt       = percent.mt   <= mt_hi,
    pass_all      = pass_nFeature & pass_nCount & pass_mt
  )

filter_report <- md %>%
  group_by(stage, replicate) %>%
  summarize(
    n_in        = n(),
    n_out_nFeat = sum(!pass_nFeature),
    n_out_nCount = sum(!pass_nCount),
    n_out_mt    = sum(!pass_mt),
    n_pass      = sum(pass_all),
    pct_pass    = round(100 * mean(pass_all), 1),
    .groups = "drop"
  )
message("\nFiltering report:")
print(filter_report)
write.csv(filter_report, file.path(tbl_dir, "02_filter_report.csv"),
          row.names = FALSE)

overall <- md %>%
  group_by(stage) %>%
  summarize(n_in = n(), n_pass = sum(pass_all),
            pct_pass = round(100 * mean(pass_all), 1), .groups = "drop")
message("\nOverall retention per stage:")
print(overall)

cells_keep <- md$Cell[md$pass_all]
seu <- subset(seu, cells = cells_keep)
message("\n>>> Cells retained: ", ncol(seu), " / ", nrow(md),
        " (", round(100 * ncol(seu) / nrow(md), 1), "%)")

# ------------------------------------------------------------
# 8. Post-filter QC dashboard
# ------------------------------------------------------------
md_post <- seu@meta.data
md_post$stage <- factor(md_post$stage, levels = c("MBC", "prePB", "PB", "PC"))

p1b <- make_violin(md_post, "nFeature_RNA", "Genes per cell")        + ggtitle("Gene complexity")
p2b <- make_violin(md_post, "nCount_RNA",   "UMIs per cell", TRUE)   + ggtitle("UMI counts (log10)")
p3b <- make_violin(md_post, "percent.mt",   "% mt (protein-coding)") + ggtitle("Mitochondrial content")
p4b <- make_violin(md_post, "percent.ig",   "% immunoglobulin")      + ggtitle("Immunoglobulin content")

qc_dashboard_post <- (p1b | p2b | p3b | p4b) +
  plot_annotation(
    title = "scRNA-seq QC dashboard - post-filtering",
    subtitle = paste0("Retained ", ncol(seu), " cells (",
                      round(100 * ncol(seu) / nrow(md), 1),
                      "%) across 4 differentiation stages"),
    theme = theme(
      plot.title = element_text(face = "bold", size = 15, color = "#2a2622"),
      plot.subtitle = element_text(size = 11, color = "#5c544b")
    )
  )

ggsave(file.path(fig_dir, "03_qc_dashboard_postfilter.pdf"),
       qc_dashboard_post, width = 13, height = 4.5)
ggsave(file.path(fig_dir, "03_qc_dashboard_postfilter.png"),
       qc_dashboard_post, width = 13, height = 4.5, dpi = 200)
print(qc_dashboard_post)

# ------------------------------------------------------------
# 9. Save filtered object
# ------------------------------------------------------------
saveRDS(seu, file.path(out_dir, "seurat_merged_qc_filtered.rds"))

message("\n========================================")
message("Stage 2.2 (QC) complete!")
message("Cells before QC: ", nrow(md))
message("Cells after QC:  ", ncol(seu))
message("Retention:       ", round(100 * ncol(seu) / nrow(md), 1), "%")
message("========================================")
