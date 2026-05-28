# ============================================================
# 12_atac_data_ingest.R   (v2 - count matrices only)
#
# Phase 3.1: Load scATAC-seq data using the precomputed peak x cell
# count matrices (the fragments files have coordinate corruption,
# see methodology note in case study writeup).
#
# Strategy:
#   1. Load count matrix + peak list per stage
#   2. Build a unified peak set across stages (GRanges reduce)
#   3. Map each stage's peaks -> unified peaks via overlap,
#      sum counts of overlapping source peaks into unified peaks
#   4. Merge into single Signac/Seurat object
#   5. Basic QC metrics (peak count per cell, fragments proxy)
#   6. Save unified object for downstream clustering
#
# Inputs (data/raw/GSE242324/):
#   - GSM7758{111,112,113,114}_{stage}_count_matrix.csv.gz
#   - GSM7758{111,112,113,114}_{stage}_peak_names_out.csv.gz
#
# Outputs:
#   - data/processed/atac_unified_peaks.rds
#   - data/processed/atac_signac_merged.rds
#   - results/figures/29_atac_peak_overlap.{pdf,png}
#   - results/figures/30_atac_qc_dashboard.{pdf,png}
#   - results/tables/11_atac_peak_counts.csv
#   - results/tables/11_atac_qc_per_stage.csv
#
# Reference: hg19
# Note: TSS enrichment and nucleosome signal QC not computed,
#       fragments-derived metrics unavailable due to source data issues.
# ============================================================

suppressPackageStartupMessages({
  library(Signac)
  library(Seurat)
  library(GenomicRanges)
  library(EnsDb.Hsapiens.v86)
  library(Matrix)
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

set.seed(42)

raw_dir <- "data/raw/GSE242324"
out_dir <- "data/processed"
fig_dir <- "results/figures"
tbl_dir <- "results/tables"

stage_colors <- c("MBC"   = "#7d8c6e",
                  "prePB" = "#b88a3e",
                  "PB"    = "#b15835",
                  "PC"    = "#7c5c6b")

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
# Helper: parse "chr1:713973-714223" strings into GRanges
# ------------------------------------------------------------
parse_peaks_to_granges <- function(peak_strings) {
  parts <- strsplit(peak_strings, "[:-]")
  GRanges(
    seqnames = sapply(parts, `[`, 1),
    ranges = IRanges(
      start = as.integer(sapply(parts, `[`, 2)),
      end   = as.integer(sapply(parts, `[`, 3))
    )
  )
}

# ------------------------------------------------------------
# Sample manifest
# ------------------------------------------------------------
samples <- data.frame(
  stage = c("MBC", "prePB", "PB", "PC"),
  gsm   = c("GSM7758111", "GSM7758112", "GSM7758113", "GSM7758114"),
  stringsAsFactors = FALSE
)
samples$peaks_file <- file.path(raw_dir,
  paste0(samples$gsm, "_", samples$stage, "_peak_names_out.csv.gz"))
samples$counts_file <- file.path(raw_dir,
  paste0(samples$gsm, "_", samples$stage, "_count_matrix.csv.gz"))

# Verify files exist
for (i in seq_len(nrow(samples))) {
  for (col in c("peaks_file", "counts_file")) {
    if (!file.exists(samples[[col]][i])) {
      stop("Missing: ", samples[[col]][i])
    }
  }
}
message("All input files found.")

# ------------------------------------------------------------
# 1. Load peaks + counts per stage
# ------------------------------------------------------------
message("\n[1/6] Loading peaks and count matrices per stage...")

stage_data <- list()
for (i in seq_len(nrow(samples))) {
  stage <- samples$stage[i]
  message("\n  Loading ", stage, "...")

  # Peaks (chr:start-end strings)
  peak_strings <- fread(samples$peaks_file[i], header = FALSE)$V1
  gr <- parse_peaks_to_granges(peak_strings)
  message("    Peaks: ", length(gr))

  # Count matrix (peaks x cells, comma-separated, header = cell barcodes)
  message("    Reading count matrix...")
  cm <- fread(samples$counts_file[i], header = TRUE, sep = ",")
  cm_mat <- as.matrix(cm)
  rownames(cm_mat) <- peak_strings
  cm_sparse <- as(cm_mat, "CsparseMatrix")
  rm(cm, cm_mat); gc()
  message("    Matrix: ", nrow(cm_sparse), " peaks x ",
          ncol(cm_sparse), " cells")

  stage_data[[stage]] <- list(
    peaks    = gr,
    counts   = cm_sparse,
    peak_ids = peak_strings
  )
}

# ------------------------------------------------------------
# 2. Build unified peak set
# ------------------------------------------------------------
message("\n[2/6] Building unified peak set...")

all_peaks <- Reduce(c, lapply(stage_data, `[[`, "peaks"))
unified_peaks <- reduce(all_peaks)

# Keep only standard chromosomes
standard_chr <- paste0("chr", c(1:22, "X", "Y"))
unified_peaks <- unified_peaks[seqnames(unified_peaks) %in% standard_chr]
seqlevels(unified_peaks) <- standard_chr

# Filter very wide or very narrow peaks
peak_widths <- width(unified_peaks)
unified_peaks <- unified_peaks[peak_widths < 10000 & peak_widths > 20]

# Name unified peaks as chr-start-end (Signac convention uses hyphen)
unified_peak_names <- paste0(
  as.character(seqnames(unified_peaks)), "-",
  start(unified_peaks), "-",
  end(unified_peaks)
)
names(unified_peaks) <- unified_peak_names

message("Unified peaks: ", length(unified_peaks))
message("  Median width: ", median(width(unified_peaks)), " bp")

saveRDS(unified_peaks, file.path(out_dir, "atac_unified_peaks.rds"))

# Peak counts table
peak_counts_table <- data.frame(
  Stage  = c(samples$stage, "Unified"),
  Peaks  = c(sapply(stage_data, function(x) length(x$peaks)),
             length(unified_peaks))
)
write.csv(peak_counts_table,
          file.path(tbl_dir, "11_atac_peak_counts.csv"),
          row.names = FALSE)
print(peak_counts_table)

# ------------------------------------------------------------
# 3. Re-project each stage's counts onto unified peaks
# ------------------------------------------------------------
# Strategy: for each stage, find which original peaks overlap each
# unified peak. Sum the counts of overlapping source peaks into
# the unified peak. This aggregates correctly: if 2 narrow source
# peaks overlap one wider unified peak, the cell's count for that
# unified peak is the sum.
# ------------------------------------------------------------
message("\n[3/6] Re-projecting counts onto unified peaks...")

stage_objects <- list()

for (stage in names(stage_data)) {
  message("\n  Re-projecting ", stage, "...")

  src_peaks  <- stage_data[[stage]]$peaks
  src_counts <- stage_data[[stage]]$counts
  src_peak_ids <- stage_data[[stage]]$peak_ids

  # Find overlaps: which source peak overlaps which unified peak
  hits <- findOverlaps(src_peaks, unified_peaks)
  src_idx <- queryHits(hits)
  uni_idx <- subjectHits(hits)

  # Build sparse projection matrix: rows = unified peaks, cols = source peaks
  # M[i, j] = 1 if source peak j overlaps unified peak i
  proj <- sparseMatrix(
    i = uni_idx,
    j = src_idx,
    x = 1,
    dims = c(length(unified_peaks), length(src_peaks)),
    dimnames = list(unified_peak_names, src_peak_ids)
  )

  # New counts = projection %*% source counts
  # (rows are unified peaks, cols are cells)
  new_counts <- proj %*% src_counts

  message("    Re-projected: ", nrow(new_counts), " unified peaks x ",
          ncol(new_counts), " cells")

  # Build Signac ChromatinAssay (without fragments)
  chrom_assay <- CreateChromatinAssay(
    counts       = new_counts,
    sep          = c("-", "-"),
    genome       = "hg19",
    min.cells    = 5,
    min.features = 200
  )

  # Build Seurat object
  # Important: ChromatinAssay drops cells below min.features, so use
  # its colnames (not the raw matrix colnames) when building metadata.
  surviving_cells <- colnames(chrom_assay)
  meta_df <- data.frame(
    row.names = surviving_cells,
    stage     = rep(stage, length(surviving_cells)),
    gsm       = rep(samples$gsm[samples$stage == stage],
                    length(surviving_cells))
  )

  seu_stage <- CreateSeuratObject(
    counts   = chrom_assay,
    assay    = "peaks",
    meta.data = meta_df
  )

  stage_objects[[stage]] <- seu_stage
  message("    Object: ", ncol(seu_stage), " cells, ",
          nrow(seu_stage), " peaks (after Signac min.cells/min.features)")

  # Free memory: don't need the raw stage data anymore
  rm(new_counts, proj, src_counts); gc()
}

# Free the raw count matrices (we keep the peak ranges via unified_peaks)
rm(stage_data); gc()

# ------------------------------------------------------------
# 4. Merge stages into one Signac object
# ------------------------------------------------------------
message("\n[4/6] Merging stages...")

# Prefix cell names with stage so they're unique across the merge
for (stage in names(stage_objects)) {
  stage_objects[[stage]] <- RenameCells(stage_objects[[stage]],
                                         add.cell.id = stage)
}

seu_atac <- merge(
  x = stage_objects[[1]],
  y = stage_objects[-1]
)
seu_atac$stage <- factor(seu_atac$stage,
                          levels = c("MBC", "prePB", "PB", "PC"))

message("Merged object: ", ncol(seu_atac), " cells, ",
        nrow(seu_atac), " peaks")
print(table(seu_atac$stage))

# Free per-stage objects
rm(stage_objects); gc()

# ------------------------------------------------------------
# 5. Add hg19 gene annotations
# ------------------------------------------------------------
message("\n[5/6] Adding hg19 gene annotations from EnsDb.Hsapiens.v86...")
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- "UCSC"
genome(annotations) <- "hg19"
Annotation(seu_atac) <- annotations

# ------------------------------------------------------------
# 6. Basic QC metrics (without fragments)
# ------------------------------------------------------------
message("\n[6/6] Computing basic QC metrics...")

# Per-cell stats already in metadata via Seurat:
#   nCount_peaks    = total counts per cell (fragments-in-peaks proxy)
#   nFeature_peaks  = unique peaks detected per cell

qc_table <- seu_atac@meta.data %>%
  group_by(stage) %>%
  summarize(
    n_cells          = n(),
    median_peaks     = median(nFeature_peaks),
    median_counts    = median(nCount_peaks),
    iqr_counts_low   = quantile(nCount_peaks, 0.25),
    iqr_counts_high  = quantile(nCount_peaks, 0.75),
    .groups = "drop"
  )
write.csv(qc_table,
          file.path(tbl_dir, "11_atac_qc_per_stage.csv"),
          row.names = FALSE)
print(qc_table)

# ------------------------------------------------------------
# QC dashboard figure
# ------------------------------------------------------------
p1 <- ggplot(seu_atac@meta.data,
             aes(x = stage, y = log10(nCount_peaks + 1), fill = stage)) +
  geom_violin(alpha = 0.7, color = NA) +
  geom_boxplot(width = 0.15, fill = "white", color = "#2a2622",
               outlier.shape = NA, alpha = 0.9) +
  scale_fill_manual(values = stage_colors) +
  theme_portfolio() +
  labs(title = "Total counts per cell (log10)", x = NULL,
       y = "log10(nCount_peaks + 1)") +
  theme(legend.position = "none")

p2 <- ggplot(seu_atac@meta.data,
             aes(x = stage, y = nFeature_peaks, fill = stage)) +
  geom_violin(alpha = 0.7, color = NA) +
  geom_boxplot(width = 0.15, fill = "white", color = "#2a2622",
               outlier.shape = NA, alpha = 0.9) +
  scale_fill_manual(values = stage_colors) +
  theme_portfolio() +
  labs(title = "Unique peaks per cell", x = NULL,
       y = "Peaks detected") +
  theme(legend.position = "none")

p3 <- ggplot(seu_atac@meta.data,
             aes(x = nCount_peaks, y = nFeature_peaks, color = stage)) +
  geom_point(alpha = 0.4, size = 0.4) +
  scale_color_manual(values = stage_colors) +
  scale_x_log10(labels = scales::comma) +
  scale_y_log10(labels = scales::comma) +
  theme_portfolio() +
  labs(title = "Counts vs unique peaks (per cell)",
       x = "log10 counts", y = "log10 unique peaks") +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))

# Cell counts per stage
cell_counts <- as.data.frame(table(seu_atac$stage))
colnames(cell_counts) <- c("stage", "n_cells")

p4 <- ggplot(cell_counts, aes(x = stage, y = n_cells, fill = stage)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = stage_colors) +
  geom_text(aes(label = scales::comma(n_cells)),
            vjust = -0.4, size = 3.5, color = "#2a2622") +
  scale_y_continuous(labels = scales::comma,
                     expand = expansion(mult = c(0, 0.15))) +
  theme_portfolio() +
  labs(title = "Cells per stage after QC", x = NULL,
       y = "Cells") +
  theme(legend.position = "none")

qc_dashboard <- (p1 | p2) / (p3 | p4) +
  plot_annotation(
    title    = "Phase 3.1: scATAC-seq QC across differentiation stages",
    subtitle = paste0(ncol(seu_atac), " cells, ", nrow(seu_atac),
                      " peaks (unified)")
  )

ggsave(file.path(fig_dir, "30_atac_qc_dashboard.pdf"),
       qc_dashboard, width = 12, height = 8)
ggsave(file.path(fig_dir, "30_atac_qc_dashboard.png"),
       qc_dashboard, width = 12, height = 8, dpi = 200)

# Peak set overlap figure
peak_overlap <- data.frame(
  Stage = factor(peak_counts_table$Stage,
                 levels = c("MBC", "prePB", "PB", "PC", "Unified")),
  Peaks = peak_counts_table$Peaks
)

p_overlap <- ggplot(peak_overlap, aes(x = Stage, y = Peaks, fill = Stage)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = c(stage_colors, "Unified" = "#2a2622")) +
  geom_text(aes(label = scales::comma(Peaks)),
            vjust = -0.5, size = 3.5, color = "#2a2622") +
  scale_y_continuous(labels = scales::comma,
                     expand = expansion(mult = c(0, 0.15))) +
  theme_portfolio() +
  labs(title = "Peaks per stage and unified peak set",
       subtitle = "Peaks were called separately per stage, then reduced to a unified set",
       x = NULL, y = "Number of peaks") +
  theme(legend.position = "none")

ggsave(file.path(fig_dir, "29_atac_peak_overlap.pdf"), p_overlap,
       width = 8, height = 5)
ggsave(file.path(fig_dir, "29_atac_peak_overlap.png"), p_overlap,
       width = 8, height = 5, dpi = 200)

# ------------------------------------------------------------
# 7. Save final merged object
# ------------------------------------------------------------
message("\nSaving final object...")
saveRDS(seu_atac, file.path(out_dir, "atac_signac_merged.rds"))

message("\n========================================")
message("Phase 3.1 complete!")
message("\nNumbers:")
message("  Cells per stage:")
print(table(seu_atac$stage))
message("  Total peaks (unified): ", nrow(seu_atac))
message("\nFigures:")
message("  29_atac_peak_overlap")
message("  30_atac_qc_dashboard")
message("\nTables:")
message("  11_atac_peak_counts.csv")
message("  11_atac_qc_per_stage.csv")
message("\nObject saved:")
message("  atac_signac_merged.rds (~", round(ncol(seu_atac)/1000, 1),
        "k cells, ", round(nrow(seu_atac)/1000), "k peaks)")
message("\nNote: TSS enrichment and nucleosome signal were not computed")
message("      due to coordinate corruption in the deposited fragments files.")
message("      Cell-level QC is based on peak count metrics from the count matrix.")
message("\nNext: 13_atac_dimreduction.R (TF-IDF, LSI, UMAP, clustering)")
message("========================================")
