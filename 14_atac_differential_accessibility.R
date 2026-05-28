# ============================================================
# 14_atac_differential_accessibility.R   (final, all fixes)
#
# Phase 3.3: Differential accessibility via sequential pair
# transitions along the differentiation trajectory.
#
# Comparisons:
#   T1: MBC vs prePB   (activation entry)
#   T2: prePB vs PB    (plasmablast specification)
#   T3: PB vs PC       (terminal differentiation)
#
# Each transition produces:
#   - "opening" peaks (gained accessibility moving forward)
#   - "closing" peaks (lost accessibility moving forward)
#
# Method: logistic regression on per-cell accessibility,
# with nCount_peaks as latent covariate.
#
# Final version fixes:
#   - conflicted package namespace resolution
#   - explicit dplyr:: prefixes throughout
#   - robust peak annotation via genes() + forced UCSC seqlevels
#   - saves da_results as checkpoint to avoid re-running FindMarkers
#
# Inputs:
#   - data/processed/atac_signac_clustered.rds
#
# Outputs:
#   - results/tables/13_atac_da_transitions_full.csv
#   - results/tables/13_atac_da_top_per_transition.csv
#   - results/tables/13_atac_da_summary.csv
#   - results/tables/13_atac_da_peak_annotation.csv
#   - results/figures/35_atac_da_volcano_per_transition.{pdf,png}
#   - results/figures/36_atac_da_transition_summary.{pdf,png}
#   - results/figures/37_atac_da_peak_annotation.{pdf,png}
#   - data/processed/atac_da_results_raw.rds       (checkpoint)
#   - data/processed/atac_da_peaks_per_transition.rds (GRanges)
# ============================================================

suppressPackageStartupMessages({
  library(Signac)
  library(Seurat)
  library(GenomicRanges)
  library(EnsDb.Hsapiens.v86)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

# Resolve namespace conflicts that block Signac and dplyr operations.
# Without this, the conflicted package (loaded transitively) prevents
# Signac's FindMarkers from using intersect() and breaks downstream
# dplyr pipelines.
if (requireNamespace("conflicted", quietly = TRUE)) {
  conflicted::conflicts_prefer(base::intersect, .quiet = TRUE)
  conflicted::conflicts_prefer(base::union, .quiet = TRUE)
  conflicted::conflicts_prefer(base::setdiff, .quiet = TRUE)
  conflicted::conflicts_prefer(dplyr::filter, .quiet = TRUE)
  conflicted::conflicts_prefer(dplyr::select, .quiet = TRUE)
  conflicted::conflicts_prefer(dplyr::rename, .quiet = TRUE)
  conflicted::conflicts_prefer(dplyr::mutate, .quiet = TRUE)
  conflicted::conflicts_prefer(dplyr::desc, .quiet = TRUE)
  conflicted::conflicts_prefer(dplyr::lag, .quiet = TRUE)
}

set.seed(42)

out_dir <- "data/processed"
fig_dir <- "results/figures"
tbl_dir <- "results/tables"

stage_colors <- c("MBC"   = "#7d8c6e",
                  "prePB" = "#b88a3e",
                  "PB"    = "#b15835",
                  "PC"    = "#7c5c6b")

# The three transitions
transitions <- list(
  T1_MBC_to_prePB  = list(from = "MBC",   to = "prePB"),
  T2_prePB_to_PB   = list(from = "prePB", to = "PB"),
  T3_PB_to_PC      = list(from = "PB",    to = "PC")
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
# 1. Load clustered object
# ------------------------------------------------------------
message("Loading clustered ATAC object...")
seu <- readRDS(file.path(out_dir, "atac_signac_clustered.rds"))
DefaultAssay(seu) <- "peaks"
message("Cells: ", ncol(seu), "  Peaks: ", nrow(seu))
print(table(seu$stage))

Idents(seu) <- "stage"

# ------------------------------------------------------------
# 2. Run differential accessibility per transition
# ------------------------------------------------------------
# Skip if checkpoint exists. FindMarkers is slow (~30-45 min per
# transition); we never want to re-run it accidentally.
# ------------------------------------------------------------
checkpoint_path <- file.path(out_dir, "atac_da_results_raw.rds")

if (file.exists(checkpoint_path)) {
  message("\nFound checkpoint at ", checkpoint_path,
          ". Loading instead of re-running FindMarkers.")
  message("(Delete the file if you want to recompute.)")
  da_results <- readRDS(checkpoint_path)
} else {
  message("\nRunning differential accessibility per transition...")
  message("(approximately 30-45 minutes per transition; total ~90 min)")
  
  da_results <- list()
  
  for (tname in names(transitions)) {
    t_info <- transitions[[tname]]
    from_stage <- t_info$from
    to_stage   <- t_info$to
    
    message("\n--- ", tname, " (", from_stage, " -> ", to_stage, ") ---")
    
    res <- FindMarkers(
      seu,
      ident.1         = to_stage,    # forward state
      ident.2         = from_stage,  # backward state
      test.use        = "LR",
      latent.vars     = "nCount_peaks",
      min.pct         = 0.05,
      only.pos        = FALSE,
      logfc.threshold = 0.1
    )
    
    res$peak        <- rownames(res)
    res$transition  <- tname
    res$from_stage  <- from_stage
    res$to_stage    <- to_stage
    res$direction   <- ifelse(res$avg_log2FC > 0, "opening", "closing")
    res$significant <- res$p_val_adj < 0.05 & abs(res$avg_log2FC) > 0.25
    
    da_results[[tname]] <- res
    
    message("  Peaks tested: ", nrow(res))
    message("  OPENING (gained in ", to_stage, "): ",
            sum(res$significant & res$direction == "opening"))
    message("  CLOSING (lost in ", to_stage, "): ",
            sum(res$significant & res$direction == "closing"))
  }
  
  saveRDS(da_results, checkpoint_path)
  message("\nCheckpoint saved: ", checkpoint_path)
}

# ------------------------------------------------------------
# 3. Combined results table
# ------------------------------------------------------------
da_full <- do.call(rbind, da_results)
write.csv(da_full,
          file.path(tbl_dir, "13_atac_da_transitions_full.csv"),
          row.names = FALSE)

# Top 200 opening per transition
top_per_transition <- da_full |>
  dplyr::filter(significant == TRUE, direction == "opening") |>
  dplyr::group_by(transition) |>
  dplyr::arrange(p_val_adj, dplyr::desc(avg_log2FC)) |>
  dplyr::slice_head(n = 200) |>
  dplyr::ungroup()
write.csv(top_per_transition,
          file.path(tbl_dir, "13_atac_da_top_per_transition.csv"),
          row.names = FALSE)

# Summary
summary_table <- da_full |>
  dplyr::group_by(transition) |>
  dplyr::summarize(
    peaks_tested = dplyr::n(),
    opening      = sum(significant == TRUE & direction == "opening"),
    closing      = sum(significant == TRUE & direction == "closing"),
    median_log2FC_opening = round(
      median(avg_log2FC[significant == TRUE & direction == "opening"]), 2),
    median_log2FC_closing = round(
      median(avg_log2FC[significant == TRUE & direction == "closing"]), 2),
    .groups = "drop"
  )
write.csv(summary_table,
          file.path(tbl_dir, "13_atac_da_summary.csv"),
          row.names = FALSE)
message("\nSummary:")
print(summary_table)

# ------------------------------------------------------------
# 4. Convert significant peaks to GRanges per transition
# ------------------------------------------------------------
message("\nConverting significant peaks to GRanges per transition...")

da_granges <- list()
for (tname in names(transitions)) {
  for (dir in c("opening", "closing")) {
    sig_peaks <- da_results[[tname]] |>
      dplyr::filter(significant == TRUE, direction == dir) |>
      dplyr::pull(peak)
    
    key <- paste(tname, dir, sep = "_")
    
    if (length(sig_peaks) == 0) {
      message("  ", key, ": 0 peaks. Skipping.")
      next
    }
    
    parts <- strsplit(sig_peaks, "-")
    gr <- GRanges(
      seqnames = sapply(parts, `[`, 1),
      ranges = IRanges(
        start = as.integer(sapply(parts, `[`, 2)),
        end   = as.integer(sapply(parts, `[`, 3))
      )
    )
    gr$transition <- tname
    gr$direction  <- dir
    gr$peak_id    <- sig_peaks
    
    rows <- da_results[[tname]][sig_peaks, ]
    gr$log2FC <- rows$avg_log2FC
    gr$padj   <- rows$p_val_adj
    
    da_granges[[key]] <- gr
    message("  ", key, ": ", length(gr), " peaks")
  }
}

saveRDS(da_granges, file.path(out_dir, "atac_da_peaks_per_transition.rds"))

# ------------------------------------------------------------
# 5. Volcano plots per transition
# ------------------------------------------------------------
message("\nGenerating volcano plots...")

da_full$log10p <- -log10(da_full$p_val_adj + 1e-300)
da_full$log10p[da_full$log10p > 50] <- 50

volcano_list <- list()
for (tname in names(transitions)) {
  t_info <- transitions[[tname]]
  df <- da_full |> dplyr::filter(transition == tname)
  
  n_open  <- sum(df$significant == TRUE & df$direction == "opening")
  n_close <- sum(df$significant == TRUE & df$direction == "closing")
  
  df$plot_color <- "gray80"
  df$plot_color[df$significant & df$direction == "opening"] <-
    stage_colors[t_info$to]
  df$plot_color[df$significant & df$direction == "closing"] <-
    stage_colors[t_info$from]
  df$plot_alpha <- ifelse(df$significant, 0.8, 0.3)
  
  p <- ggplot(df, aes(x = avg_log2FC, y = log10p)) +
    geom_point(aes(color = plot_color, alpha = plot_alpha), size = 0.6) +
    scale_color_identity() +
    scale_alpha_identity() +
    geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed",
               color = "gray40", linewidth = 0.3) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed",
               color = "gray40", linewidth = 0.3) +
    theme_portfolio() +
    labs(title    = paste0(t_info$from, " -> ", t_info$to),
         subtitle = paste0(n_open, " opening, ", n_close,
                           " closing (FDR < 0.05, |log2FC| > 0.25)"),
         x = paste0("log2FC (", t_info$to, " vs ", t_info$from, ")"),
         y = "-log10(p adj)")
  volcano_list[[tname]] <- p
}

volcano_combined <- wrap_plots(volcano_list, ncol = 3) +
  plot_annotation(
    title = "Chromatin accessibility transitions across the differentiation trajectory",
    subtitle = "Sequential pairwise comparisons. Positive log2FC = peaks gaining accessibility moving forward."
  )

ggsave(file.path(fig_dir, "35_atac_da_volcano_per_transition.pdf"),
       volcano_combined, width = 14, height = 5.5)
ggsave(file.path(fig_dir, "35_atac_da_volcano_per_transition.png"),
       volcano_combined, width = 14, height = 5.5, dpi = 200)

# ------------------------------------------------------------
# 6. Transition summary bar chart (opening vs closing)
# ------------------------------------------------------------
message("\nGenerating transition summary figure...")

transition_long <- summary_table |>
  dplyr::select(transition, opening, closing) |>
  tidyr::pivot_longer(c(opening, closing),
                      names_to = "direction", values_to = "n_peaks") |>
  dplyr::mutate(n_peaks_signed = ifelse(direction == "closing",
                                        -n_peaks, n_peaks))

transition_long$transition_label <- factor(transition_long$transition,
                                           levels = names(transitions),
                                           labels = c("MBC -> prePB", "prePB -> PB", "PB -> PC"))

p_summary <- ggplot(transition_long,
                    aes(x = transition_label, y = n_peaks_signed,
                        fill = direction)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = scales::comma(n_peaks),
                vjust = ifelse(direction == "opening", -0.4, 1.4)),
            size = 3.5, color = "#2a2622") +
  geom_hline(yintercept = 0, color = "#2a2622", linewidth = 0.4) +
  scale_fill_manual(values = c("opening" = "#b15835",
                               "closing" = "#7d8c6e")) +
  scale_y_continuous(labels = function(x) scales::comma(abs(x))) +
  theme_portfolio() +
  labs(title    = "Number of peaks changing per developmental transition",
       subtitle = "Up = chromatin opening (gained accessibility). Down = chromatin closing.",
       x = NULL, y = "Number of significant peaks") +
  theme(legend.title = element_blank())

ggsave(file.path(fig_dir, "36_atac_da_transition_summary.pdf"), p_summary,
       width = 9, height = 6)
ggsave(file.path(fig_dir, "36_atac_da_transition_summary.png"), p_summary,
       width = 9, height = 6, dpi = 200)

# ------------------------------------------------------------
# 7. Peak annotation (robust version)
# ------------------------------------------------------------
# Uses genes() function from EnsDb directly (not the generic
# annotations) and forces UCSC seqlevels so chr1, chr2, ... match
# between peak GRanges and gene GRanges.
# ------------------------------------------------------------
message("\nAnnotating peaks (promoter / gene body / intergenic)...")

genes_gr <- genes(EnsDb.Hsapiens.v86)
seqlevelsStyle(genes_gr) <- "UCSC"
genome(genes_gr) <- "hg19"
message("  Genes loaded: ", length(genes_gr))

# Restrict to standard chromosomes
standard_chr <- paste0("chr", c(1:22, "X", "Y"))
genes_gr <- genes_gr[as.character(seqnames(genes_gr)) %in% standard_chr]

proms_gr <- promoters(genes_gr, upstream = 2000, downstream = 200)

annot_summary_list <- list()
for (key in names(da_granges)) {
  gr <- da_granges[[key]]
  
  # Match seqlevels (only standard chromosomes)
  gr_filt <- gr[as.character(seqnames(gr)) %in% standard_chr]
  
  is_promoter   <- overlapsAny(gr_filt, proms_gr)
  is_gene_body  <- overlapsAny(gr_filt, genes_gr) & !is_promoter
  is_intergenic <- !is_promoter & !is_gene_body
  
  key_parts <- strsplit(key, "_(?=opening|closing)", perl = TRUE)[[1]]
  annot_summary_list[[key]] <- data.frame(
    transition = key_parts[1],
    direction  = key_parts[2],
    promoter   = sum(is_promoter),
    gene_body  = sum(is_gene_body),
    intergenic = sum(is_intergenic),
    total      = length(gr_filt)
  )
  message("  ", key, ": promoter=", sum(is_promoter),
          "  gene_body=", sum(is_gene_body),
          "  intergenic=", sum(is_intergenic))
}

annot_combined <- do.call(rbind, annot_summary_list)
write.csv(annot_combined,
          file.path(tbl_dir, "13_atac_da_peak_annotation.csv"),
          row.names = FALSE)
print(annot_combined)

# Annotation figure
annot_long <- annot_combined |>
  tidyr::pivot_longer(c(promoter, gene_body, intergenic),
                      names_to = "feature", values_to = "n") |>
  dplyr::group_by(transition, direction) |>
  dplyr::mutate(pct = 100 * n / sum(n)) |>
  dplyr::ungroup() |>
  dplyr::mutate(
    transition_label = factor(transition,
                              levels = names(transitions),
                              labels = c("MBC -> prePB", "prePB -> PB", "PB -> PC")),
    feature = factor(feature,
                     levels = c("promoter", "gene_body", "intergenic")),
    direction = factor(direction, levels = c("opening", "closing")))

p_annot <- ggplot(annot_long, aes(x = direction, y = pct, fill = feature)) +
  geom_col(width = 0.7) +
  facet_wrap(~transition_label, ncol = 3) +
  scale_fill_manual(values = c("promoter"   = "#b15835",
                               "gene_body"  = "#b88a3e",
                               "intergenic" = "#7d8c6e"),
                    name = "Peak location") +
  theme_portfolio() +
  labs(title    = "Genomic distribution of transition peaks",
       subtitle = "Promoter = within 2 kb of TSS. Distal/intergenic = candidate enhancer regions.",
       x = NULL, y = "% of significant peaks")

ggsave(file.path(fig_dir, "37_atac_da_peak_annotation.pdf"), p_annot,
       width = 11, height = 5.5)
ggsave(file.path(fig_dir, "37_atac_da_peak_annotation.png"), p_annot,
       width = 11, height = 5.5, dpi = 200)

# ------------------------------------------------------------
# Done
# ------------------------------------------------------------
message("\n========================================")
message("Phase 3.3 complete!")
message("\nSummary:")
print(summary_table)
message("\nFigures:")
message("  35_atac_da_volcano_per_transition  (3-panel volcano)")
message("  36_atac_da_transition_summary      (opening/closing bar chart)")
message("  37_atac_da_peak_annotation         (genomic distribution)")
message("\nTables:")
message("  13_atac_da_transitions_full.csv")
message("  13_atac_da_top_per_transition.csv")
message("  13_atac_da_summary.csv")
message("  13_atac_da_peak_annotation.csv")
message("\nCheckpoints saved (for re-runs):")
message("  atac_da_results_raw.rds       (skip FindMarkers if re-running)")
message("  atac_da_peaks_per_transition.rds  (for motif analysis)")
message("\nNext: 15_atac_motif_analysis.R (requires BSgenome.Hsapiens.UCSC.hg19)")
message("========================================")