# ============================================================
# 19_hdwgcna_tf_regulatory_network.R   (final, all fixes)
#
# Builds a TF regulatory network using hdWGCNA's ConstructTFNetwork
# pipeline:
#   1. MotifScan: scan each gene's promoter for JASPAR motifs
#   2. ConstructTFNetwork: xgboost TF -> target prediction
#   3. AssignTFRegulons: produce non-overlapping regulons
#   4. Visualization: regulon bar plots, TF x module heatmap
#
# Adapted from hdWGCNA tutorial:
#   https://smorabit.github.io/hdWGCNA/articles/regulatory_network.html
#
# Critical corrections for this dataset:
#   - hg19 reference (Alaterre et al. used hg19)
#   - EnsDb.Hsapiens.v75 (last Ensembl release on GRCh37/hg19)
#   - JASPAR2020 (consistent with Phase 3.4 ATAC motif analysis)
#   - Chromosome-boundary safety filter (avoid MotifScan errors)
#   - TF list derived from JASPAR PFMs (not from motif_df summary,
#     which only contains per-motif aggregate counts)
#   - 'all_cells' column added to BOTH Seurat and metacell objects
#     (SetDatExpr operates on the metacell internal object)
#   - MotifScan checkpoint: skip if seu already has motif data
#
# Inputs:
#   - data/processed/seurat_hdWGCNA_modules.rds
#
# Outputs:
#   - data/processed/seurat_hdwgcna_tf_network.rds        (final object)
#   - data/processed/seurat_hdwgcna_post_motifscan.rds    (post-scan checkpoint)
#   - results/tables/17_motif_scan_summary.csv
#   - results/tables/17_tf_network_results.csv
#   - results/tables/17_regulon_assignments.csv
#   - results/tables/17_regulon_module_overlap.csv
#   - results/figures/46_tf_count_summary.{pdf,png}
#   - results/figures/47_top_regulons_barplots.{pdf,png}
#   - results/figures/48_regulon_module_overlap.{pdf,png}
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(hdWGCNA)
  library(JASPAR2020)
  library(motifmatchr)
  library(TFBSTools)
  library(EnsDb.Hsapiens.v75)
  library(BSgenome.Hsapiens.UCSC.hg19)
  library(GenomicRanges)
  library(xgboost)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
})

# Namespace conflict resolution
if (requireNamespace("conflicted", quietly = TRUE)) {
  conflicted::conflicts_prefer(base::intersect, .quiet = TRUE)
  conflicted::conflicts_prefer(base::setdiff, .quiet = TRUE)
  conflicted::conflicts_prefer(base::union, .quiet = TRUE)
  conflicted::conflicts_prefer(dplyr::filter, .quiet = TRUE)
  conflicted::conflicts_prefer(dplyr::select, .quiet = TRUE)
  conflicted::conflicts_prefer(dplyr::mutate, .quiet = TRUE)
  conflicted::conflicts_prefer(dplyr::desc, .quiet = TRUE)
  conflicted::conflicts_prefer(dplyr::slice, .quiet = TRUE)
}

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

module_colors <- c(
  "yellow"    = "#e8b800",
  "brown"     = "#964b00",
  "turquoise" = "#4cb5b0",
  "blue"      = "#4472c4"
)

# ------------------------------------------------------------
# 1. Load hdWGCNA object (or post-MotifScan checkpoint)
# ------------------------------------------------------------
motifscan_checkpoint <- file.path(out_dir,
                                  "seurat_hdwgcna_post_motifscan.rds")

if (file.exists(motifscan_checkpoint)) {
  message("Found MotifScan checkpoint. Loading...")
  message("(Delete the file to recompute from scratch.)")
  seu <- readRDS(motifscan_checkpoint)
  skip_motifscan <- TRUE
} else {
  message("Loading hdWGCNA object...")
  seu <- readRDS(file.path(out_dir, "seurat_hdWGCNA_modules.rds"))
  skip_motifscan <- FALSE
}

DefaultAssay(seu) <- "RNA"
message("Cells: ", ncol(seu))
print(table(seu$stage))

# ------------------------------------------------------------
# 2. Load JASPAR2020 vertebrate motifs
# ------------------------------------------------------------
message("\nLoading JASPAR2020 vertebrate motifs...")
opts <- list()
opts[["collection"]]   <- "CORE"
opts[["tax_group"]]    <- "vertebrates"
opts[["all_versions"]] <- FALSE

pfm_core <- getMatrixSet(JASPAR2020, opts = opts)
message("Loaded ", length(pfm_core), " motifs")

# Derive the canonical TF gene-symbol list from JASPAR PFM names.
# Compound names (FOS::JUN) are split, variant suffixes ((var.2)) are
# stripped, and symbols are normalized to uppercase for matching.
jaspar_tf_raw <- sapply(pfm_core, function(x) name(x))
jaspar_tfs    <- unique(unlist(strsplit(jaspar_tf_raw, "::")))
jaspar_tfs    <- unique(gsub("\\(.*\\)$", "", jaspar_tfs))
jaspar_tfs_upper <- toupper(jaspar_tfs)
message("Unique JASPAR TF symbols: ", length(jaspar_tfs_upper))

# ------------------------------------------------------------
# 3. Build chromosome-safety gene filter
# ------------------------------------------------------------
# Some genes sit close to chromosome ends. A 2 kb upstream
# promoter window can extend past the chromosome boundary,
# causing MotifScan to crash. Pre-filter to genes whose
# promoter is fully within bounds. With v75 + hg19 this
# typically drops zero genes (the assemblies match exactly).
# Belt-and-suspenders insurance.
# ------------------------------------------------------------
message("\nBuilding chromosome-safety gene filter...")

hg19_genome <- BSgenome.Hsapiens.UCSC.hg19
chr_lengths <- seqlengths(hg19_genome)

all_genes_gr <- genes(EnsDb.Hsapiens.v75)
seqlevelsStyle(all_genes_gr) <- "UCSC"
genome(all_genes_gr) <- "hg19"

standard_chr <- paste0("chr", c(1:22, "X", "Y"))
all_genes_gr <- all_genes_gr[
  as.character(seqnames(all_genes_gr)) %in% standard_chr
]

PROMOTER_UP   <- 2000
PROMOTER_DOWN <- 500
proms <- promoters(all_genes_gr,
                   upstream   = PROMOTER_UP,
                   downstream = PROMOTER_DOWN)
prom_chr_lengths <- chr_lengths[as.character(seqnames(proms))]
keep_safe <- start(proms) >= 1 & end(proms) <= prom_chr_lengths

safe_gene_names <- all_genes_gr$gene_name[keep_safe]
safe_gene_names <- unique(safe_gene_names[
  !is.na(safe_gene_names) & safe_gene_names != ""
])

message("  v75 genes (standard chr): ", length(all_genes_gr))
message("  Safe-promoter genes: ", sum(keep_safe))
message("  Unique safe symbols: ", length(safe_gene_names))

# ------------------------------------------------------------
# 4. Module assignment and gene universe
# ------------------------------------------------------------
modules <- GetModules(seu)
nongrey_genes <- subset(modules, module != "grey")$gene_name
all_data_genes <- rownames(seu)

# TFs present in this dataset (intersect JASPAR TFs with our genes)
tfs_in_data_upper <- intersect(toupper(all_data_genes), jaspar_tfs_upper)
tfs_in_data <- all_data_genes[toupper(all_data_genes) %in% tfs_in_data_upper]

# Final universe: TFs + non-grey module genes, intersected with
# safe-promoter genes
expanded_universe <- intersect(
  unique(c(tfs_in_data, nongrey_genes)),
  intersect(safe_gene_names, all_data_genes)
)

message("\nGene universe summary:")
message("  Non-grey module genes:        ", length(nongrey_genes))
message("  JASPAR TFs present in data:   ", length(tfs_in_data))
message("  Total universe (safe-filtered):", length(expanded_universe))

# ------------------------------------------------------------
# 5. MotifScan (skip if checkpoint exists)
# ------------------------------------------------------------
if (!skip_motifscan) {
  message("\nRunning MotifScan (hg19, EnsDb v75)...")
  message("(10-30 minutes; will save a checkpoint when done)")
  
  seu <- SetWGCNAGenes(seu, expanded_universe)
  
  seu <- MotifScan(
    seu,
    species_genome = "hg19",
    pfm            = pfm_core,
    EnsDb          = EnsDb.Hsapiens.v75
  )
  
  # Save the checkpoint so we never have to re-run MotifScan
  saveRDS(seu, motifscan_checkpoint)
  message("MotifScan complete. Checkpoint saved.")
} else {
  message("\nSkipping MotifScan (using checkpoint)")
  # Make sure the universe is still set on the loaded object
  seu <- SetWGCNAGenes(seu, expanded_universe)
}

active_wgcna <- seu@misc$active_wgcna

# ---- 1. Build motif_name -> gene_name mapping ----
# motif_name examples: "Arnt", "Ahr::Arnt", "MZF1(var.2)", "FOS::JUN"
# We split compound names, strip variant suffixes, uppercase for matching
# Then map each compound motif to the FIRST gene symbol that is present
# in the data (or keep multiple rows if both are present).

motif_df <- GetMotifs(seu)
message("Original motif_df columns: ", paste(colnames(motif_df), collapse = ", "))
message("Rows: ", nrow(motif_df))

# Genes present in the Seurat data, uppercased for matching
data_genes <- rownames(seu)
data_genes_upper <- toupper(data_genes)
upper_to_data <- setNames(data_genes, data_genes_upper)

# Expand each motif row into one row per gene symbol it contains
motif_df_expanded <- motif_df %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    tf_symbols = list({
      raw <- motif_name
      parts <- unlist(strsplit(raw, "::"))
      parts <- gsub("\\(.*\\)$", "", parts)
      toupper(parts)
    })
  ) %>%
  tidyr::unnest(tf_symbols) %>%
  dplyr::mutate(gene_name = upper_to_data[tf_symbols]) %>%
  dplyr::filter(!is.na(gene_name)) %>%
  dplyr::ungroup() %>%
  as.data.frame()

message("Expanded motif_df rows: ", nrow(motif_df_expanded))
message("Unique gene_name (TFs mapped to data): ",
        length(unique(motif_df_expanded$gene_name)))

# Quick sanity check
print(head(motif_df_expanded[, c("motif_name", "motif_ID",
                                 "tf_symbols", "gene_name")], 10))

# ---- 2. Write the expanded motif_df back to the Seurat object ----
# hdWGCNA stores it via SetMotifs (some versions) or directly in misc.
# Try SetMotifs first, fall back to direct assignment.

if (exists("SetMotifs") && is.function(SetMotifs)) {
  seu <- SetMotifs(seu, motif_df_expanded)
  message("Updated via SetMotifs()")
} else {
  seu@misc[[active_wgcna]]$motif_info <- motif_df_expanded
  message("Updated via direct slot assignment")
}

# Verify
check <- GetMotifs(seu)
message("After update: ", nrow(check), " rows, columns: ",
        paste(colnames(check), collapse = ", "))
message("'gene_name' column present: ", "gene_name" %in% colnames(check))

# ------------------------------------------------------------
# 6. SetDatExpr on the full trajectory
# ------------------------------------------------------------
# SetDatExpr operates on the internal metacell Seurat object,
# which has its own metadata slot. We need to add 'all_cells'
# to BOTH the outer object and the metacell object.
# ------------------------------------------------------------
message("\nSetting up trajectory-wide expression pool...")


# Update the metacell object too


seu <- SetDatExpr(
  seu,
  group_name  = c("MBC", "prePB", "PB", "PC"),
  group.by= "stage",
  assay      = "RNA"
)
message("Expression matrix set")

# ------------------------------------------------------------
# 7. Construct the TF network with xgboost
# ------------------------------------------------------------
message("\nConstructing TF network with xgboost (30-60 min)...")

n_cores <-8
message("Using ", n_cores, " threads")

model_params <- list(
  objective = "reg:squarederror",
  max_depth = 1,
  eta       = 0.1,
  nthread   = n_cores,
  alpha     = 0.5
)

seu <- ConstructTFNetwork2(seu, model_params = model_params)

tf_net_results <- GetTFNetwork(seu)
message("\nTF network results:")
message("  Total TF -> target predictions: ", nrow(tf_net_results))
print(head(tf_net_results))

write.csv(tf_net_results,
          file.path(tbl_dir, "17_tf_network_results.csv"),
          row.names = FALSE)

# Save intermediate state (so AssignTFRegulons can be re-tuned
# without re-running the expensive xgboost step)
saveRDS(seu,
        file.path(out_dir, "seurat_hdwgcna_tf_network_pre_regulons.rds"))

# ------------------------------------------------------------
# 8. AssignTFRegulons (strategy C, non-overlapping)
# ------------------------------------------------------------
message("\nAssigning regulons (strategy C, threshold 0.01)...")

seu <- AssignTFRegulons(
  seu,
  strategy   = "C",
  reg_thresh = 0.01,
  n_tfs      = 10
)

regulons <- GetTFRegulons(seu)
message("\nRegulon assignments:")
message("  Total regulon edges: ", nrow(regulons))
message("  Unique regulator TFs: ", length(unique(regulons$tf)))
message("  Unique target genes: ", length(unique(regulons$gene_name)))
print(head(regulons))

write.csv(regulons,
          file.path(tbl_dir, "17_regulon_assignments.csv"),
          row.names = FALSE)

regulon_summary <- regulons |>
  dplyr::group_by(tf) |>
  dplyr::summarize(
    n_targets       = dplyr::n(),
    mean_importance = mean(score, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::arrange(dplyr::desc(n_targets))

message("\nTop regulators by regulon size:")
print(head(regulon_summary, 20))

# ------------------------------------------------------------
# 9. Save the final enriched Seurat object
# ------------------------------------------------------------
saveRDS(seu, file.path(out_dir, "seurat_hdwgcna_tf_network.rds"))
message("\nObject saved: seurat_hdwgcna_tf_network.rds")

# ------------------------------------------------------------
# 10. Figure 46: top 30 TFs by regulon size
# ------------------------------------------------------------
message("\nGenerating top-TF summary figure...")

p_tf_count <- ggplot(head(regulon_summary, 30),
                     aes(x = reorder(tf, n_targets), y = n_targets)) +
  geom_col(fill = "#b15835", width = 0.75) +
  coord_flip() +
  theme_portfolio() +
  labs(
    title    = "Top 30 transcription factors by regulon size",
    subtitle = "Number of target genes assigned to each TF (xgboost importance > 0.01)",
    x = NULL, y = "Number of target genes"
  )

ggsave(file.path(fig_dir, "46_tf_count_summary.pdf"), p_tf_count,
       width = 8, height = 9)
ggsave(file.path(fig_dir, "46_tf_count_summary.png"), p_tf_count,
       width = 8, height = 9, dpi = 200)

# ------------------------------------------------------------
# 11. Figure 47: regulon bar plots for key TFs
# ------------------------------------------------------------
message("\nGenerating regulon bar plots for key TFs...")

candidate_tfs <- c("BATF3", "IRF3", "STAT2", "POU2F2", "XBP1",
                   "JUN", "FOS", "JUNB", "EGR1", "PAX5", "SPIB",
                   "KLF6", "ATF7", "BHLHE40", "RELB")

available_tfs <- intersect(candidate_tfs, unique(regulons$tf))
message("Candidate TFs with regulons in this data: ",
        paste(available_tfs, collapse = ", "))

if (length(available_tfs) >= 2) {
  display_tfs <- head(available_tfs, 6)
  bar_plots <- list()
  
  for (tf_name in display_tfs) {
    tryCatch({
      p <- RegulonBarPlot(seu, selected_tf = tf_name) +
        theme_portfolio() +
        theme(plot.title = element_text(size = 10))
      bar_plots[[tf_name]] <- p
    }, error = function(e) {
      message("  Skipping ", tf_name, ": ", e$message)
    })
  }
  
  if (length(bar_plots) > 0) {
    n_plots  <- length(bar_plots)
    n_rows   <- ceiling(n_plots / 2)
    combined <- wrap_plots(bar_plots, ncol = 2) +
      plot_annotation(
        title    = "Top regulons of key transcription factors",
        subtitle = "Target gene xgboost importance scores per TF"
      )
    
    ggsave(file.path(fig_dir, "47_top_regulons_barplots.pdf"),
           combined, width = 12, height = 4 * n_rows)
    ggsave(file.path(fig_dir, "47_top_regulons_barplots.png"),
           combined, width = 12, height = 4 * n_rows, dpi = 200)
  }
}

# ------------------------------------------------------------
# 12. Figure 48: TF regulons x hdWGCNA module heatmap
# ------------------------------------------------------------
message("\nComputing regulon-module overlap (headline figure)...")

regulons_with_module <- regulons |>
  dplyr::left_join(
    modules |> dplyr::select(gene_name, module),
    by = "gene_name"
  ) |>
  dplyr::filter(!is.na(module), module != "grey")

overlap_summary <- regulons_with_module |>
  dplyr::group_by(tf, module) |>
  dplyr::summarize(n_targets = dplyr::n(), .groups = "drop")

tfs_with_module_signal <- overlap_summary |>
  dplyr::group_by(tf) |>
  dplyr::summarize(max_targets = max(n_targets), .groups = "drop") |>
  dplyr::filter(max_targets >= 5) |>
  dplyr::pull(tf)

overlap_filtered <- overlap_summary |>
  dplyr::filter(tf %in% tfs_with_module_signal) |>
  dplyr::arrange(tf, dplyr::desc(n_targets))

message("TFs with >= 5 targets in any module: ",
        length(tfs_with_module_signal))

write.csv(overlap_filtered,
          file.path(tbl_dir, "17_regulon_module_overlap.csv"),
          row.names = FALSE)

if (nrow(overlap_filtered) > 0) {
  top_tfs_for_heatmap <- overlap_filtered |>
    dplyr::group_by(tf) |>
    dplyr::summarize(total = sum(n_targets), .groups = "drop") |>
    dplyr::arrange(dplyr::desc(total)) |>
    dplyr::slice_head(n = 25) |>
    dplyr::pull(tf)
  
  heatmap_df <- overlap_filtered |>
    dplyr::filter(tf %in% top_tfs_for_heatmap) |>
    dplyr::mutate(
      module = factor(module,
                      levels = c("yellow", "brown", "turquoise", "blue"))
    )
  
  tf_dominant_module <- heatmap_df |>
    dplyr::group_by(tf) |>
    dplyr::slice_max(n_targets, n = 1, with_ties = FALSE) |>
    dplyr::arrange(module, dplyr::desc(n_targets)) |>
    dplyr::pull(tf)
  
  heatmap_df$tf <- factor(heatmap_df$tf, levels = rev(tf_dominant_module))
  
  p_heat <- ggplot(heatmap_df,
                   aes(x = module, y = tf, fill = n_targets)) +
    geom_tile(color = "white", linewidth = 0.4) +
    geom_text(aes(label = n_targets),
              color = "#2a2622", size = 3) +
    scale_fill_gradient(low = "#faf7f2", high = "#b15835",
                        name = "Targets in\nmodule") +
    theme_minimal(base_size = 9) +
    theme(
      plot.title    = element_text(face = "bold", size = 12,
                                   color = "#2a2622"),
      plot.subtitle = element_text(size = 10, color = "#5c544b"),
      axis.text.y   = element_text(size = 9),
      axis.text.x   = element_text(angle = 0, hjust = 0.5),
      panel.grid    = element_blank()
    ) +
    labs(
      title    = "TF regulons mapped to hdWGCNA modules",
      subtitle = "Each cell: number of a TF's targets belonging to each module",
      x = NULL, y = NULL
    )
  
  ggsave(file.path(fig_dir, "48_regulon_module_overlap.pdf"),
         p_heat, width = 7,
         height = max(7, length(top_tfs_for_heatmap) * 0.3))
  ggsave(file.path(fig_dir, "48_regulon_module_overlap.png"),
         p_heat, width = 7,
         height = max(7, length(top_tfs_for_heatmap) * 0.3),
         dpi = 200)
}

# ------------------------------------------------------------
# Done
# ------------------------------------------------------------
message("\n========================================")
message("hdWGCNA TF regulatory network complete!")
message("\nReference: hg19 + EnsDb.Hsapiens.v75")
message("\nTables:")
message("  17_motif_scan_summary.csv    (", nrow(motif_df),
        " motifs scanned)")
message("  17_tf_network_results.csv    (", nrow(tf_net_results),
        " predictions)")
message("  17_regulon_assignments.csv   (", nrow(regulons), " edges)")
message("  17_regulon_module_overlap.csv (", nrow(overlap_filtered),
        " TF-module overlaps)")
message("\nFigures:")
message("  46_tf_count_summary           (top regulators)")
message("  47_top_regulons_barplots      (key TFs from ATAC + RNA)")
message("  48_regulon_module_overlap     (TF x module heatmap)")
message("\nObjects saved:")
message("  seurat_hdwgcna_post_motifscan.rds (MotifScan checkpoint)")
message("  seurat_hdwgcna_tf_network_pre_regulons.rds (xgboost done)")
message("  seurat_hdwgcna_tf_network.rds (final)")
message("========================================")