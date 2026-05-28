# ============================================================
# 15_atac_motif_analysis.R
#
# Phase 3.4: Motif enrichment analysis on transition-opening peaks.
#
# Question: which TFs drive each developmental transition?
#
# Strategy:
#   1. Add JASPAR2020 motifs to the Signac object
#   2. Run motif enrichment per transition's opening peaks vs
#      a background of all unified peaks
#   3. Generate heatmaps and volcanos showing motif enrichment
#   4. Cross-reference with hdWGCNA module hub genes to identify
#      cross-modal TF candidates (TF in chromatin + TF in RNA module)
#
# Method: motifmatchr scans each peak's underlying DNA sequence
# against JASPAR PWMs. Then we test whether each motif is enriched
# in opening peaks vs the genomic background.
#
# Inputs:
#   - data/processed/atac_signac_clustered.rds
#   - data/processed/atac_da_peaks_per_transition.rds
#   - data/processed/seurat_hdWGCNA_modules.rds (for cross-reference)
#
# Outputs:
#   - results/tables/14_motif_enrichment_per_transition.csv
#   - results/tables/14_motif_module_crossref.csv
#   - results/figures/38_motif_heatmap_top_per_transition.{pdf,png}
#   - results/figures/39_motif_volcano_per_transition.{pdf,png}
#   - results/figures/40_motif_module_crossref.{pdf,png}
#   - data/processed/atac_signac_with_motifs.rds
# ============================================================

suppressPackageStartupMessages({
  library(Signac)
  library(Seurat)
  library(GenomicRanges)
  library(JASPAR2020)
  library(TFBSTools)
  library(motifmatchr)
  library(BSgenome.Hsapiens.UCSC.hg19)
  library(EnsDb.Hsapiens.v86)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

# Namespace conflict resolution
if (requireNamespace("conflicted", quietly = TRUE)) {
  conflicted::conflicts_prefer(base::intersect, .quiet = TRUE)
  conflicted::conflicts_prefer(base::union, .quiet = TRUE)
  conflicted::conflicts_prefer(base::setdiff, .quiet = TRUE)
  conflicted::conflicts_prefer(dplyr::filter, .quiet = TRUE)
  conflicted::conflicts_prefer(dplyr::select, .quiet = TRUE)
  conflicted::conflicts_prefer(dplyr::rename, .quiet = TRUE)
  conflicted::conflicts_prefer(dplyr::mutate, .quiet = TRUE)
  conflicted::conflicts_prefer(dplyr::desc, .quiet = TRUE)
}

set.seed(42)

out_dir <- "data/processed"
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
# 1. Load inputs
# ------------------------------------------------------------
message("Loading inputs...")
seu <- readRDS(file.path(out_dir, "atac_signac_clustered.rds"))
DefaultAssay(seu) <- "peaks"

da_granges <- readRDS(file.path(out_dir, "atac_da_peaks_per_transition.rds"))

message("Cells: ", ncol(seu), "  Peaks: ", nrow(seu))
message("Transition peak sets loaded: ", length(da_granges))

# ------------------------------------------------------------
# 2. Load JASPAR2020 vertebrate motifs
# ------------------------------------------------------------
# We use the JASPAR2020 CORE collection for vertebrates,
# which contains ~750 high-confidence TF binding motifs.
# ------------------------------------------------------------
message("\nLoading JASPAR2020 vertebrate motifs...")
opts <- list()
opts[["collection"]] <- "CORE"
opts[["tax_group"]]  <- "vertebrates"
opts[["all_versions"]] <- FALSE

pwm_list <- getMatrixSet(JASPAR2020, opts = opts)
message("Loaded ", length(pwm_list), " motifs")

# Get motif names mapped to TF gene symbols
motif_names <- sapply(pwm_list, function(x) name(x))
message("Sample motifs: ",
        paste(head(motif_names, 5), collapse = ", "))

# ------------------------------------------------------------
# 3. Add motif information to Signac object
# ------------------------------------------------------------
# AddMotifs scans each peak's DNA sequence and creates a
# binary motif x peak matrix indicating which motifs are
# present in which peaks.
# Takes ~5-10 minutes for ~300k peaks.
# ------------------------------------------------------------
message("\nAdding motifs to Signac object (takes ~5-10 min)...")
message("This scans every unified peak for ", length(pwm_list), " motifs.")

seu <- AddMotifs(
  object  = seu,
  genome  = BSgenome.Hsapiens.UCSC.hg19,
  pfm     = pwm_list,
  verbose = TRUE
)

message("Motif annotation complete.")

# Save the object so we don't need to repeat this expensive step
saveRDS(seu, file.path(out_dir, "atac_signac_with_motifs.rds"))
message("Saved: atac_signac_with_motifs.rds")

# ------------------------------------------------------------
# 4. Motif enrichment per transition
# ------------------------------------------------------------
# For each opening peak set, FindMotifs tests whether each
# motif is enriched in the opening peaks compared to a
# background set of similar peaks (matched for sequence
# composition).
# ------------------------------------------------------------
message("\nRunning motif enrichment per transition...")

# Background: all unified peaks (Signac will match for composition)
all_peak_ids <- rownames(seu)
message("Background peak universe: ", length(all_peak_ids), " peaks")

enrichment_results <- list()

for (tname in names(transitions)) {
  for (dir in c("opening", "closing")) {
    key <- paste(tname, dir, sep = "_")
    
    if (is.null(da_granges[[key]])) {
      message("  ", key, ": no peaks. Skipping.")
      next
    }
    
    target_peaks <- da_granges[[key]]$peak_id
    # Filter to peaks that survived Signac filtering (in the object)
    target_peaks <- intersect(target_peaks, all_peak_ids)
    
    if (length(target_peaks) < 50) {
      message("  ", key, ": only ", length(target_peaks),
              " peaks survived. Skipping.")
      next
    }
    
    message("\n--- ", key, " (", length(target_peaks), " peaks) ---")
    
    enrich <- FindMotifs(
      object   = seu,
      features = target_peaks
    )
    
    enrich$transition <- tname
    enrich$direction  <- dir
    enrich$key        <- key
    
    enrichment_results[[key]] <- enrich
    
    # Show top 10
    message("Top 10 enriched motifs:")
    print(head(enrich[, c("motif.name", "observed", "percent.observed",
                          "fold.enrichment", "pvalue")], 10))
  }
}

# ------------------------------------------------------------
# 5. Combine and save full results
# ------------------------------------------------------------
enrichment_full <- do.call(rbind, enrichment_results)

# Add adjusted p-value across each transition (BH)
enrichment_full <- enrichment_full |>
  dplyr::group_by(key) |>
  dplyr::mutate(padj = p.adjust(pvalue, method = "BH")) |>
  dplyr::ungroup() |>
  dplyr::mutate(
    significant = padj < 0.05 & fold.enrichment > 1.5,
    log2_fold_enrichment = log2(fold.enrichment + 1e-6)
  )

write.csv(enrichment_full,
          file.path(tbl_dir, "14_motif_enrichment_per_transition.csv"),
          row.names = FALSE)

# Per-transition summary
enrichment_summary <- enrichment_full |>
  dplyr::group_by(key) |>
  dplyr::summarize(
    motifs_tested = dplyr::n(),
    significant   = sum(significant),
    .groups = "drop"
  )
message("\nEnrichment summary:")
print(enrichment_summary)

# ------------------------------------------------------------
# 6. Heatmap of top motifs per transition (opening only)
# ------------------------------------------------------------
message("\nGenerating motif heatmap (top per transition)...")

# Pick top 20 enriched motifs per opening transition by p-value
top_motifs <- enrichment_full |>
  dplyr::filter(direction == "opening", significant) |>
  dplyr::group_by(transition) |>
  dplyr::arrange(pvalue, dplyr::desc(fold.enrichment)) |>
  dplyr::slice_head(n = 20) |>
  dplyr::ungroup()

top_motif_set <- unique(top_motifs$motif.name)
message("Top motif set size: ", length(top_motif_set))

# Build matrix: rows = motifs, cols = transitions, values = log2FE
heatmap_data <- enrichment_full |>
  dplyr::filter(direction == "opening",
                motif.name %in% top_motif_set) |>
  dplyr::select(transition, motif.name, log2_fold_enrichment) |>
  tidyr::pivot_wider(names_from = transition,
                     values_from = log2_fold_enrichment,
                     values_fill = 0)

heatmap_mat <- as.matrix(heatmap_data[, -1])
rownames(heatmap_mat) <- heatmap_data$motif.name

# Order rows by which transition shows max enrichment, then by value
row_assign <- apply(heatmap_mat, 1, which.max)
row_order  <- order(row_assign, -apply(heatmap_mat, 1, max))
heatmap_mat <- heatmap_mat[row_order, , drop = FALSE]

# Long format for ggplot
heatmap_df <- as.data.frame(heatmap_mat) |>
  tibble::rownames_to_column("motif") |>
  tidyr::pivot_longer(-motif, names_to = "transition", values_to = "log2FE")

heatmap_df$motif <- factor(heatmap_df$motif, levels = rev(rownames(heatmap_mat)))
heatmap_df$transition <- factor(heatmap_df$transition,
                                levels = names(transitions),
                                labels = c("MBC -> prePB", "prePB -> PB", "PB -> PC"))

p_heatmap <- ggplot(heatmap_df,
                    aes(x = transition, y = motif, fill = log2FE)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.1f", log2FE)),
            color = "#2a2622", size = 2.7) +
  scale_fill_gradient2(low = "#3b6286", mid = "#faf7f2", high = "#b15835",
                       midpoint = 0, name = "log2 FE") +
  theme_minimal(base_size = 9) +
  theme(
    plot.title    = element_text(face = "bold", size = 12, color = "#2a2622"),
    plot.subtitle = element_text(size = 10, color = "#5c544b"),
    axis.text.y   = element_text(size = 8),
    axis.text.x   = element_text(angle = 0, hjust = 0.5),
    panel.grid    = element_blank()
  ) +
  labs(title    = "Top enriched motifs in opening peaks per transition",
       subtitle = "log2 fold enrichment (motif occurrence in opening peaks vs background)",
       x = NULL, y = NULL)

ggsave(file.path(fig_dir, "38_motif_heatmap_top_per_transition.pdf"),
       p_heatmap, width = 8, height = max(7, length(top_motif_set) * 0.18))
ggsave(file.path(fig_dir, "38_motif_heatmap_top_per_transition.png"),
       p_heatmap, width = 8, height = max(7, length(top_motif_set) * 0.18),
       dpi = 200)

# ------------------------------------------------------------
# 7. Volcano per transition (opening direction only)
# ------------------------------------------------------------
message("\nGenerating motif volcanos per transition...")

volcano_list <- list()
for (tname in names(transitions)) {
  key <- paste(tname, "opening", sep = "_")
  if (is.null(enrichment_results[[key]])) next
  
  df <- enrichment_full |>
    dplyr::filter(key == !!key) |>
    dplyr::mutate(neg_log10p = -log10(pvalue + 1e-300))
  df$neg_log10p[df$neg_log10p > 50] <- 50
  
  # Label top 10 most significant enriched motifs
  top_labels <- df |>
    dplyr::filter(significant) |>
    dplyr::arrange(pvalue) |>
    dplyr::slice_head(n = 10)
  
  p <- ggplot(df, aes(x = log2_fold_enrichment, y = neg_log10p)) +
    geom_point(aes(color = significant, alpha = significant),
               size = 1.2) +
    scale_color_manual(values = c("FALSE" = "gray70", "TRUE" = "#b15835")) +
    scale_alpha_manual(values = c("FALSE" = 0.4, "TRUE" = 0.85)) +
    geom_text(data = top_labels,
              aes(label = motif.name),
              size = 2.8, color = "#2a2622",
              hjust = -0.1, vjust = -0.2,
              check_overlap = TRUE) +
    geom_vline(xintercept = log2(1.5), linetype = "dashed",
               color = "gray40", linewidth = 0.3) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed",
               color = "gray40", linewidth = 0.3) +
    theme_portfolio() +
    labs(title = paste0(transitions[[tname]]$from, " -> ",
                        transitions[[tname]]$to),
         subtitle = paste0(sum(df$significant), " significant motifs"),
         x = "log2 fold enrichment", y = "-log10(p-value)") +
    theme(legend.position = "none")
  volcano_list[[tname]] <- p
}

volcano_combined <- wrap_plots(volcano_list, ncol = 3) +
  plot_annotation(
    title = "TF motif enrichment in opening peaks per transition",
    subtitle = "Each point is a JASPAR2020 motif. Right side = enriched in opening peaks."
  )

ggsave(file.path(fig_dir, "39_motif_volcano_per_transition.pdf"),
       volcano_combined, width = 14, height = 5.5)
ggsave(file.path(fig_dir, "39_motif_volcano_per_transition.png"),
       volcano_combined, width = 14, height = 5.5, dpi = 200)

# ------------------------------------------------------------
# 8. Cross-reference with hdWGCNA modules
# ------------------------------------------------------------
# Identify motifs whose corresponding TF gene is a hub or
# member of one of our hdWGCNA modules. This finds the
# TFs that are evidenced from BOTH chromatin (motif in opening
# peaks) AND RNA (hub gene in a co-expression module).
# ------------------------------------------------------------
message("\nCross-referencing motifs with hdWGCNA modules...")

hdwgcna_path <- file.path(out_dir, "seurat_hdWGCNA_modules.rds")
if (!file.exists(hdwgcna_path)) {
  message("hdWGCNA object not found; skipping cross-reference.")
  modules_df <- NULL
} else {
  seu_rna <- readRDS(hdwgcna_path)
  suppressPackageStartupMessages(library(hdWGCNA))
  modules_df <- GetModules(seu_rna)
  message("hdWGCNA modules loaded: ",
          length(unique(modules_df$module)), " modules.")
  rm(seu_rna); gc()
}

if (!is.null(modules_df)) {
  # Join motif TF name to module membership
  # Motif names from JASPAR include compound names like "FOS::JUND"
  # We split these to also match the components
  
  # Build a TF -> module lookup
  module_lookup <- modules_df |>
    dplyr::select(gene_name, module) |>
    dplyr::filter(module != "grey")
  
  # Find significant motifs in opening transitions
  enrich_open <- enrichment_full |>
    dplyr::filter(direction == "opening", significant)
  
  # For each motif, split compound names and check membership
  enrich_open$tf_components <- strsplit(enrich_open$motif.name, "::")
  enrich_open$matched_modules <- sapply(enrich_open$tf_components, function(tfs) {
    tfs_upper <- toupper(tfs)
    found <- module_lookup$module[module_lookup$gene_name %in% tfs_upper]
    if (length(found) == 0) return(NA_character_)
    paste(unique(found), collapse = ",")
  })
  
  # Filter to those with module matches
  crossref <- enrich_open |>
    dplyr::filter(!is.na(matched_modules)) |>
    dplyr::select(transition, motif.name, fold.enrichment,
                  pvalue, padj, matched_modules) |>
    dplyr::arrange(transition, pvalue)
  
  write.csv(crossref,
            file.path(tbl_dir, "14_motif_module_crossref.csv"),
            row.names = FALSE)
  
  message("\nMotifs matching hdWGCNA modules:")
  print(crossref)
  
  # Visualization: motif by module by transition
  if (nrow(crossref) > 0) {
    crossref$transition_label <- factor(crossref$transition,
                                        levels = names(transitions),
                                        labels = c("MBC -> prePB", "prePB -> PB", "PB -> PC"))
    
    p_crossref <- ggplot(crossref,
                         aes(x = transition_label, y = motif.name,
                             size = -log10(pvalue + 1e-300),
                             color = matched_modules)) +
      geom_point() +
      scale_size_continuous(name = "-log10(p)", range = c(2, 8)) +
      scale_color_manual(values = c("blue" = "#4472c4",
                                    "yellow" = "#e8b800",
                                    "turquoise" = "#4cb5b0",
                                    "brown" = "#964b00"),
                         name = "hdWGCNA\nmodule") +
      theme_portfolio() +
      labs(title    = "TF motifs evidenced from chromatin AND RNA co-expression",
           subtitle = "Motifs significantly enriched in opening peaks, whose TF is a member of an hdWGCNA module",
           x = NULL, y = NULL)
    
    ggsave(file.path(fig_dir, "40_motif_module_crossref.pdf"),
           p_crossref, width = 10,
           height = max(5, length(unique(crossref$motif.name)) * 0.25))
    ggsave(file.path(fig_dir, "40_motif_module_crossref.png"),
           p_crossref, width = 10,
           height = max(5, length(unique(crossref$motif.name)) * 0.25),
           dpi = 200)
  } else {
    message("No motif-module matches found.")
  }
}

# ------------------------------------------------------------
# Done
# ------------------------------------------------------------
message("\n========================================")
message("Phase 3.4 complete!")
message("\nMotifs tested:        ", length(pwm_list))
message("Enrichment tables:    14_motif_enrichment_per_transition.csv")
message("Cross-ref table:      14_motif_module_crossref.csv")
message("\nFigures:")
message("  38_motif_heatmap_top_per_transition  (top motifs by transition)")
message("  39_motif_volcano_per_transition      (enrichment volcanos)")
message("  40_motif_module_crossref             (chromatin + RNA validated TFs)")
message("\nObject saved (for downstream chromVAR if installed):")
message("  atac_signac_with_motifs.rds")
message("\nNext: 16_atac_peak_to_gene_linkage.R")
message("       or Phase 4: SCENIC TF-target inference (Python)")
message("========================================")