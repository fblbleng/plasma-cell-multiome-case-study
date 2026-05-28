# ============================================================
# 17_atac_peak_to_gene_linkage.R
#
# Phase 3.5: Proximity-based peak-to-gene linkage.
#
# Approach: assign each transition-opening peak to nearby genes
# by genomic distance to the gene's TSS. Classify by distance
# tier (promoter / proximal / distal). Then cross-reference the
# linked genes against:
#   1. hdWGCNA module membership (RNA co-expression layer)
#   2. The motif found in the peak (chromatin TF layer)
#
# The payoff: a table of "peak -> gene" links where the gene is
# also an hdWGCNA hub, AND the peak contains a motif for a TF
# active in that transition. These are triple-evidenced
# regulatory relationships (chromatin accessibility + TF motif
# + RNA co-expression).
#
# Why proximity (not correlation): the ATAC and RNA come from
# separate cell populations (GSE242324 vs GSE242330), so per-cell
# correlation is impossible. Proximity identifies candidate
# targets; cross-referencing against the RNA layer provides the
# supporting evidence.
#
# Inputs:
#   - data/processed/atac_da_peaks_per_transition.rds
#   - data/processed/atac_signac_with_motifs.rds (for motif-in-peak)
#   - data/processed/seurat_hdWGCNA_modules.rds
#
# Outputs:
#   - results/tables/16_peak_gene_links_all.csv
#   - results/tables/16_peak_gene_links_hdwgcna.csv
#   - results/tables/16_triple_evidence_links.csv
#   - results/figures/42_peak_gene_distance_distribution.{pdf,png}
#   - results/figures/43_linked_genes_per_module.{pdf,png}
# ============================================================

suppressPackageStartupMessages({
  library(Signac)
  library(Seurat)
  library(GenomicRanges)
  library(EnsDb.Hsapiens.v86)
  library(hdWGCNA)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

if (requireNamespace("conflicted", quietly = TRUE)) {
  conflicted::conflicts_prefer(base::intersect, .quiet = TRUE)
  conflicted::conflicts_prefer(base::setdiff, .quiet = TRUE)
  conflicted::conflicts_prefer(dplyr::filter, .quiet = TRUE)
  conflicted::conflicts_prefer(dplyr::select, .quiet = TRUE)
  conflicted::conflicts_prefer(dplyr::mutate, .quiet = TRUE)
  conflicted::conflicts_prefer(dplyr::desc, .quiet = TRUE)
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

transition_labels <- c(
  "T1_MBC_to_prePB" = "MBC -> prePB",
  "T2_prePB_to_PB"  = "prePB -> PB",
  "T3_PB_to_PC"     = "PB -> PC"
)

# TFs per transition (must match script 16)
tfs_per_transition <- list(
  T1_MBC_to_prePB = c("JUN::JUNB", "FOSL2::JUND(var.2)", "FOSL2::JUN(var.2)",
                      "FOSB::JUN", "ATF7", "EGR1", "CREM",
                      "CTCF", "CREB3L4(var.2)", "RFX5", "RFX7"),
  T2_prePB_to_PB  = c("POU2F2", "YY1", "KLF2", "KLF6", "KLF16",
                      "EGR3", "XBP1", "BHLHE40", "CREB3L1",
                      "SP3", "SP8", "RFX1"),
  T3_PB_to_PC     = c("IRF3", "STAT2", "RELB", "TFAP2A(var.3)",
                      "TFAP2B", "TFAP2C", "NFATC3", "HOXC9",
                      "HOXC10", "Alx4", "PLAG1", "OVOL2")
)

# ------------------------------------------------------------
# 1. Load inputs
# ------------------------------------------------------------
message("Loading inputs...")
da_granges <- readRDS(file.path(out_dir, "atac_da_peaks_per_transition.rds"))

# Gene annotation with TSS
genes_gr <- genes(EnsDb.Hsapiens.v86)
seqlevelsStyle(genes_gr) <- "UCSC"
genome(genes_gr) <- "hg19"
standard_chr <- paste0("chr", c(1:22, "X", "Y"))
genes_gr <- genes_gr[as.character(seqnames(genes_gr)) %in% standard_chr]
# Keep protein-coding + lncRNA biotypes
genes_gr <- genes_gr[genes_gr$gene_biotype %in%
                       c("protein_coding", "lincRNA", "antisense")]
message("Genes for linkage: ", length(genes_gr))

# TSS positions (single base at start of gene, strand-aware)
tss_gr <- resize(genes_gr, width = 1, fix = "start")

# hdWGCNA modules
seu_rna <- readRDS(file.path(out_dir, "seurat_hdWGCNA_modules.rds"))
modules_df <- GetModules(seu_rna)
hub_genes_table <- GetHubGenes(seu_rna, n_hubs = 25)
rm(seu_rna); gc()

module_lookup <- modules_df |>
  dplyr::filter(module != "grey") |>
  dplyr::select(gene_name, module) |>
  dplyr::distinct()
message("Genes with module assignment: ", nrow(module_lookup))

# ------------------------------------------------------------
# 2. Proximity-based linkage per transition (opening peaks)
# ------------------------------------------------------------
# For each opening peak, find genes whose TSS is within a
# distance window. Classify into distance tiers:
#   promoter  : peak overlaps TSS +/- 2kb
#   proximal  : TSS within 2-50 kb
#   distal    : TSS within 50-500 kb (candidate enhancer)
# ------------------------------------------------------------
message("\nLinking peaks to genes by proximity...")

DISTAL_WINDOW <- 500000  # 500 kb maximum

link_peaks_to_genes <- function(peaks_gr, transition_name) {
  # Find all gene TSS within DISTAL_WINDOW of each peak
  hits <- findOverlaps(
    peaks_gr,
    resize(tss_gr, width = 2 * DISTAL_WINDOW, fix = "center"),
    ignore.strand = TRUE
  )
  
  if (length(hits) == 0) return(NULL)
  
  peak_idx <- queryHits(hits)
  gene_idx <- subjectHits(hits)
  
  # Distance from peak center to TSS
  peak_centers <- start(peaks_gr) + width(peaks_gr) / 2
  tss_pos <- start(tss_gr)
  
  dist <- abs(peak_centers[peak_idx] - tss_pos[gene_idx])
  
  df <- data.frame(
    transition = transition_name,
    peak_id    = peaks_gr$peak_id[peak_idx],
    peak_log2FC= peaks_gr$log2FC[peak_idx],
    gene_name  = genes_gr$gene_name[gene_idx],
    gene_biotype = genes_gr$gene_biotype[gene_idx],
    distance_bp = dist,
    stringsAsFactors = FALSE
  )
  
  # Distance tier
  df$tier <- dplyr::case_when(
    df$distance_bp <= 2000   ~ "promoter",
    df$distance_bp <= 50000  ~ "proximal",
    TRUE                     ~ "distal"
  )
  
  # Remove rows with NA gene names
  df <- df[!is.na(df$gene_name) & df$gene_name != "", ]
  
  df
}

all_links <- list()
for (tname in names(transition_labels)) {
  key <- paste0(tname, "_opening")
  if (is.null(da_granges[[key]])) {
    message("  ", key, ": no peaks, skipping.")
    next
  }
  links <- link_peaks_to_genes(da_granges[[key]], tname)
  if (!is.null(links)) {
    all_links[[tname]] <- links
    message("  ", tname, ": ", nrow(links), " peak-gene links from ",
            length(unique(links$peak_id)), " peaks")
  }
}

links_all <- do.call(rbind, all_links)
write.csv(links_all, file.path(tbl_dir, "16_peak_gene_links_all.csv"),
          row.names = FALSE)
message("\nTotal peak-gene links: ", nrow(links_all))
message("Distance tier distribution:")
print(table(links_all$tier))

# ------------------------------------------------------------
# 3. Cross-reference with hdWGCNA modules
# ------------------------------------------------------------
message("\nCross-referencing linked genes with hdWGCNA modules...")

links_module <- links_all |>
  dplyr::inner_join(module_lookup, by = "gene_name")

write.csv(links_module,
          file.path(tbl_dir, "16_peak_gene_links_hdwgcna.csv"),
          row.names = FALSE)
message("Links where target gene is in an hdWGCNA module: ",
        nrow(links_module))
message("  Unique genes: ", length(unique(links_module$gene_name)))

# Per-transition x module summary
link_module_summary <- links_module |>
  dplyr::group_by(transition, module) |>
  dplyr::summarize(
    n_links = dplyr::n(),
    n_genes = dplyr::n_distinct(gene_name),
    n_peaks = dplyr::n_distinct(peak_id),
    .groups = "drop"
  )
message("\nLinks by transition and module:")
print(as.data.frame(link_module_summary))

# ------------------------------------------------------------
# 4. Triple-evidence links (the headline result)
# ------------------------------------------------------------
# A triple-evidenced link is one where:
#   1. A peak opens at a transition (chromatin accessibility)
#   2. The peak is near a gene that is an hdWGCNA hub (RNA)
#   3. The peak contains a motif for a TF active in that
#      transition (chromatin TF binding)
#
# We approximate (3) by checking whether the peak contains ANY
# motif from the transition's TF list, using the motif matrix
# in the with-motifs Signac object.
# ------------------------------------------------------------
message("\nIdentifying triple-evidence links...")

# Restrict to hub genes (top 25 per module)
hub_gene_set <- unique(hub_genes_table$gene_name)
links_hub <- links_module |>
  dplyr::filter(gene_name %in% hub_gene_set)
message("Links to hub genes: ", nrow(links_hub))

# Load motif object to check motif-in-peak
motif_obj_path <- file.path(out_dir, "atac_signac_with_motifs.rds")
if (file.exists(motif_obj_path)) {
  message("Loading motif object to verify motif-in-peak...")
  seu_motif <- readRDS(motif_obj_path)
  
  # Motif presence matrix: peaks x motifs (binary)
  motif_matrix <- GetMotifData(seu_motif, slot = "data")
  # Map motif IDs to TF names
  motif_name_map <- seu_motif@assays$peaks@motifs@motif.names
  rm(seu_motif); gc()
  
  # For each transition, check which hub-gene-linked peaks
  # contain a motif for that transition's TFs
  triple_links <- list()
  for (tname in names(tfs_per_transition)) {
    tfs <- tfs_per_transition[[tname]]
    
    # Find motif IDs whose name matches the transition's TFs
    matching_motif_ids <- names(motif_name_map)[
      sapply(motif_name_map, function(nm) {
        any(sapply(tfs, function(tf) {
          # Match TF name to motif name (handle compound names)
          grepl(gsub("\\(.*\\)", "", tf), nm, fixed = TRUE) ||
            nm == tf
        }))
      })
    ]
    
    if (length(matching_motif_ids) == 0) {
      message("  ", tname, ": no matching motif IDs found")
      next
    }
    
    # Hub-linked peaks for this transition
    t_links <- links_hub |> dplyr::filter(transition == tname)
    if (nrow(t_links) == 0) next
    
    # Check motif presence in each peak
    peaks_to_check <- intersect(unique(t_links$peak_id),
                                rownames(motif_matrix))
    
    if (length(peaks_to_check) == 0) next
    
    motif_submat <- motif_matrix[peaks_to_check, matching_motif_ids,
                                 drop = FALSE]
    peaks_with_motif <- peaks_to_check[Matrix::rowSums(motif_submat) > 0]
    
    # Triple-evidence links: hub-linked peaks that also contain
    # a transition-TF motif
    t_triple <- t_links |>
      dplyr::filter(peak_id %in% peaks_with_motif)
    
    if (nrow(t_triple) > 0) {
      triple_links[[tname]] <- t_triple
      message("  ", tname, ": ", nrow(t_triple),
              " triple-evidence links (",
              length(unique(t_triple$gene_name)), " genes)")
    }
  }
  
  if (length(triple_links) > 0) {
    triple_df <- do.call(rbind, triple_links)
    # Add module info and clean up
    triple_df <- triple_df |>
      dplyr::arrange(transition, module, distance_bp) |>
      dplyr::select(transition, gene_name, module, peak_id, tier,
                    distance_bp, peak_log2FC, gene_biotype)
    write.csv(triple_df,
              file.path(tbl_dir, "16_triple_evidence_links.csv"),
              row.names = FALSE)
    message("\nTriple-evidence links saved: ", nrow(triple_df))
    message("\nTop examples:")
    print(head(triple_df, 20))
  } else {
    message("No triple-evidence links found.")
    triple_df <- NULL
  }
} else {
  message("Motif object not found; skipping triple-evidence step.")
  triple_df <- NULL
}

# ------------------------------------------------------------
# 5. Figure: distance distribution
# ------------------------------------------------------------
message("\nGenerating figures...")

links_all$transition_label <- factor(
  transition_labels[links_all$transition],
  levels = transition_labels
)

p_dist <- ggplot(links_all,
                 aes(x = distance_bp / 1000, fill = tier)) +
  geom_histogram(bins = 50, color = NA) +
  facet_wrap(~transition_label, ncol = 3) +
  scale_fill_manual(values = c("promoter"  = "#b15835",
                               "proximal"  = "#b88a3e",
                               "distal"    = "#7d8c6e"),
                    name = "Distance tier") +
  scale_x_continuous(labels = function(x) paste0(x, " kb")) +
  theme_portfolio() +
  labs(title    = "Peak-to-gene distance distribution per transition",
       subtitle = "Distance from each opening peak center to the nearest gene TSS",
       x = "Distance to TSS (kb)", y = "Number of peak-gene links")

ggsave(file.path(fig_dir, "42_peak_gene_distance_distribution.pdf"),
       p_dist, width = 12, height = 5)
ggsave(file.path(fig_dir, "42_peak_gene_distance_distribution.png"),
       p_dist, width = 12, height = 5, dpi = 200)

# ------------------------------------------------------------
# 6. Figure: linked genes per module
# ------------------------------------------------------------
module_colors <- c("yellow" = "#e8b800", "brown" = "#964b00",
                   "turquoise" = "#4cb5b0", "blue" = "#4472c4")

p_mod <- ggplot(link_module_summary,
                aes(x = factor(transition,
                               levels = names(transition_labels),
                               labels = transition_labels),
                    y = n_genes, fill = module)) +
  geom_col(position = "dodge", width = 0.7) +
  scale_fill_manual(values = module_colors, name = "hdWGCNA module") +
  theme_portfolio() +
  labs(title    = "Linked module-member genes per transition",
       subtitle = "Genes near transition-opening peaks that are also hdWGCNA module members",
       x = NULL, y = "Number of unique linked genes")

ggsave(file.path(fig_dir, "43_linked_genes_per_module.pdf"),
       p_mod, width = 9, height = 5.5)
ggsave(file.path(fig_dir, "43_linked_genes_per_module.png"),
       p_mod, width = 9, height = 5.5, dpi = 200)

# ------------------------------------------------------------
# Done
# ------------------------------------------------------------
message("\n========================================")
message("Phase 3.5 complete!")
message("\nTables:")
message("  16_peak_gene_links_all.csv      (", nrow(links_all), " links)")
message("  16_peak_gene_links_hdwgcna.csv  (", nrow(links_module), " module-gene links)")
if (!is.null(triple_df)) {
  message("  16_triple_evidence_links.csv    (", nrow(triple_df), " triple-evidenced)")
}
message("\nFigures:")
message("  42_peak_gene_distance_distribution")
message("  43_linked_genes_per_module")
message("\nNext: Phase 4 SCENIC TF-target inference (Python)")
message("========================================")