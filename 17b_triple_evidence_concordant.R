# ============================================================
# 17b_triple_evidence_concordant.R
#
# Extracts the high-confidence subset of triple-evidence links:
#   1. CONCORDANT: transition's opening peaks link to genes in
#      the module that biologically matches that transition
#   2. PROXIMAL DISTANCE: only promoter (<=2kb) and proximal
#      (2-50kb) links are kept. Distal links (50-500kb) are
#      excluded from the headline because at that range, nearby
#      genes are largely a function of genomic gene density
#      rather than specific regulation, and without per-cell
#      ATAC-RNA correlation they cannot be confidently assigned.
#
# Expected transition -> module pairings:
#   T1 (MBC -> prePB):  yellow (activation / class switching)
#   T2 (prePB -> PB):   brown + turquoise (memory + proliferation)
#   T3 (PB -> PC):      blue (terminal IgG-PC / IFN response)
#
# Reuses triple_df from script 17 (or loads the saved CSV).
#
# Outputs:
#   - results/tables/16b_triple_evidence_concordant.csv
#   - results/figures/44_concordant_links_summary.{pdf,png}
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
})

if (requireNamespace("conflicted", quietly = TRUE)) {
  conflicted::conflicts_prefer(dplyr::filter, .quiet = TRUE)
  conflicted::conflicts_prefer(dplyr::select, .quiet = TRUE)
  conflicted::conflicts_prefer(dplyr::mutate, .quiet = TRUE)
  conflicted::conflicts_prefer(dplyr::desc, .quiet = TRUE)
}

tbl_dir <- "results/tables"
fig_dir <- "results/figures"

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
# 1. Load triple_df
# ------------------------------------------------------------
if (!exists("triple_df") || is.null(triple_df)) {
  triple_df <- read.csv(
    file.path(tbl_dir, "16_triple_evidence_links.csv"),
    stringsAsFactors = FALSE
  )
  message("Loaded triple_df from CSV: ", nrow(triple_df), " links")
} else {
  message("Using triple_df from session: ", nrow(triple_df), " links")
}

# ------------------------------------------------------------
# 2. Define expected transition -> module concordance
# ------------------------------------------------------------
concordant_modules <- list(
  "T1_MBC_to_prePB" = c("yellow"),
  "T2_prePB_to_PB"  = c("brown", "turquoise"),
  "T3_PB_to_PC"     = c("blue")
)

triple_df <- triple_df |>
  dplyr::rowwise() |>
  dplyr::mutate(
    concordant = module %in% concordant_modules[[transition]]
  ) |>
  dplyr::ungroup()

message("\nConcordance breakdown (all distance tiers):")
print(table(triple_df$transition, triple_df$concordant,
            dnn = c("transition", "concordant")))

# ------------------------------------------------------------
# 3. Apply BOTH filters: concordant AND promoter/proximal
# ------------------------------------------------------------
message("\nDistance tier breakdown of concordant links:")
concordant_all_tiers <- triple_df |> dplyr::filter(concordant)
print(table(concordant_all_tiers$transition, concordant_all_tiers$tier,
            dnn = c("transition", "tier")))

concordant_df <- triple_df |>
  dplyr::filter(concordant, tier %in% c("promoter", "proximal")) |>
  dplyr::arrange(transition, module, distance_bp) |>
  dplyr::select(transition, gene_name, module, peak_id, tier,
                distance_bp, peak_log2FC, gene_biotype)

write.csv(concordant_df,
          file.path(tbl_dir, "16b_triple_evidence_concordant.csv"),
          row.names = FALSE)

message("\n========================================")
message("High-confidence concordant links")
message("(concordant module + promoter/proximal distance): ",
        nrow(concordant_df))
message("========================================")

# Per transition summary
concordant_summary <- concordant_df |>
  dplyr::group_by(transition, module) |>
  dplyr::summarize(
    n_links = dplyr::n(),
    n_genes = dplyr::n_distinct(gene_name),
    genes   = paste(sort(unique(gene_name)), collapse = ", "),
    .groups = "drop"
  )

message("\nHigh-confidence links per transition:")
for (i in seq_len(nrow(concordant_summary))) {
  row <- concordant_summary[i, ]
  message("\n  ", row$transition, " (", row$module, " module): ",
          row$n_genes, " genes, ", row$n_links, " links")
  message("    Genes: ", row$genes)
}

# ------------------------------------------------------------
# 4. Figure: concordant links summary
# ------------------------------------------------------------
module_colors <- c("yellow" = "#e8b800", "brown" = "#964b00",
                   "turquoise" = "#4cb5b0", "blue" = "#4472c4")

transition_labels <- c(
  "T1_MBC_to_prePB" = "MBC -> prePB\n(AP-1 -> yellow)",
  "T2_prePB_to_PB"  = "prePB -> PB\n(POU/KLF/XBP1 ->\nbrown/turquoise)",
  "T3_PB_to_PC"     = "PB -> PC\n(IRF/STAT -> blue)"
)

plot_df <- concordant_df |>
  dplyr::group_by(transition, module, tier) |>
  dplyr::summarize(n_genes = dplyr::n_distinct(gene_name),
                   .groups = "drop") |>
  dplyr::mutate(
    transition_label = factor(transition,
                              levels = names(transition_labels),
                              labels = transition_labels),
    tier = factor(tier, levels = c("promoter", "proximal"))
  )

p_concordant <- ggplot(plot_df,
                       aes(x = transition_label, y = n_genes,
                           fill = module, alpha = tier)) +
  geom_col(width = 0.65) +
  scale_fill_manual(values = module_colors, name = "hdWGCNA module") +
  scale_alpha_manual(values = c("promoter" = 1.0, "proximal" = 0.55),
                     name = "Distance tier") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  theme_portfolio() +
  labs(
    title    = "High-confidence concordant regulatory links",
    subtitle = "Promoter/proximal links where chromatin opening, TF motif, and RNA module all agree",
    x = NULL, y = "Unique target genes"
  )

ggsave(file.path(fig_dir, "44_concordant_links_summary.pdf"),
       p_concordant, width = 9, height = 6)
ggsave(file.path(fig_dir, "44_concordant_links_summary.png"),
       p_concordant, width = 9, height = 6, dpi = 200)

# ------------------------------------------------------------
# 5. Headline examples per transition (closest, deduplicated)
# ------------------------------------------------------------
message("\n========================================")
message("HEADLINE HIGH-CONFIDENCE LINKS (closest per transition):")
message("========================================")

for (tname in names(concordant_modules)) {
  sub <- concordant_df |>
    dplyr::filter(transition == tname) |>
    dplyr::distinct(gene_name, .keep_all = TRUE) |>
    dplyr::arrange(distance_bp) |>
    dplyr::slice_head(n = 12)
  if (nrow(sub) > 0) {
    message("\n", tname, ":")
    for (j in seq_len(nrow(sub))) {
      message("  ", sub$gene_name[j], " (", sub$module[j], ", ",
              sub$tier[j], ", ", round(sub$distance_bp[j] / 1000, 1), " kb",
              ", ", sub$gene_biotype[j], ")")
    }
  } else {
    message("\n", tname, ": no high-confidence concordant links")
  }
}

message("\n========================================")
message("Phase 3.5b complete!")
message("\nTable:  16b_triple_evidence_concordant.csv (",
        nrow(concordant_df), " high-confidence links)")
message("Figure: 44_concordant_links_summary")
message("\nNote: distal links (>50kb) excluded from this headline subset.")
message("Full triple-evidence table (all tiers) remains in")
message("16_triple_evidence_links.csv for completeness.")
message("========================================")