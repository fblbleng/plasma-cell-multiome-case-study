# ============================================================
# 08_lncRNA_annotation.R   (v2 - fast GTF parser)
#
# Same purpose as v1, but uses data.table::fread to read the
# GTF directly as a TSV. This skips the GenomicRanges
# construction overhead and is ~10x faster than rtracklayer
# for large GENCODE files.
#
# Annotates lncRNAs (lincRNA + antisense biotypes) from
# GENCODE v19 (matches the hg19 alignment used by BD Rhapsody
# in Alaterre et al., Blood 2024).
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(Matrix)
})

set.seed(42)

# ------------------------------------------------------------
# 1. Paths
# ------------------------------------------------------------
ann_dir   <- "data/annotation"
out_dir   <- "data/processed"
fig_dir   <- "results/figures"
tbl_dir   <- "results/tables"

dir.create(ann_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tbl_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

gtf_url    <- "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz"
gtf_local  <- file.path(ann_dir, "gencode.v19.annotation.gtf.gz")

# ------------------------------------------------------------
# 2. Download if needed
# ------------------------------------------------------------
if (!file.exists(gtf_local)) {
  message("Downloading GENCODE v19 GTF (~30 MB)...")
  options(timeout = 600)  # 10 minutes for slow connections
  download.file(gtf_url, gtf_local, mode = "wb")
}
message("GTF: ", gtf_local)

# ------------------------------------------------------------
# 3. FAST GTF parsing with data.table::fread
# ------------------------------------------------------------
# GTF format is TSV with 9 columns. fread reads it natively
# and is much faster than rtracklayer.
# ------------------------------------------------------------
message("\nParsing GTF (this should take ~10-20 seconds)...")

gtf <- fread(gtf_local,
             sep = "\t",
             header = FALSE,
             skip = "chr",
             col.names = c("chromosome", "source", "feature_type",
                           "start", "end", "score", "strand",
                           "frame", "attributes"),
             showProgress = FALSE)

message("Total features in GTF: ", nrow(gtf))

# Keep only gene-level rows
gtf_genes <- gtf[feature_type == "gene"]
message("Gene-level records: ", nrow(gtf_genes))

# ------------------------------------------------------------
# 4. Parse attributes (biotype + symbol) with vectorized regex
# ------------------------------------------------------------
extract_attr_vec <- function(attrs, key) {
  pattern <- paste0('.*', key, ' "([^"]+)".*')
  ifelse(grepl(paste0(key, ' "'), attrs),
         sub(pattern, "\\1", attrs),
         NA_character_)
}

message("\nExtracting biotype and gene symbol from attributes...")
gtf_genes[, gene_id     := extract_attr_vec(attributes, "gene_id")]
gtf_genes[, gene_type   := extract_attr_vec(attributes, "gene_type")]
gtf_genes[, gene_symbol := extract_attr_vec(attributes, "gene_name")]
gtf_genes[, attributes := NULL]

biotype_table <- sort(table(gtf_genes$gene_type), decreasing = TRUE)
message("\nTop 10 biotypes in GENCODE v19:")
print(head(biotype_table, 10))

# ------------------------------------------------------------
# 5. Filter to focused lncRNA panel
# ------------------------------------------------------------
target_biotypes <- c("lincRNA", "antisense")
lncRNA_dt <- gtf_genes[gene_type %in% target_biotypes]
message("\nlncRNAs (lincRNA + antisense): ", nrow(lncRNA_dt))
message("  lincRNA:   ", nrow(lncRNA_dt[gene_type == "lincRNA"]))
message("  antisense: ", nrow(lncRNA_dt[gene_type == "antisense"]))

lncRNA_dt[, gene_id_clean := sub("\\..*$", "", gene_id)]

lncRNA_annotation <- as.data.frame(lncRNA_dt[, .(
  gene_id, gene_id_clean, gene_symbol, biotype = gene_type,
  chromosome, start, end, strand, width = end - start + 1
)])

message("\nFirst few rows of the annotation:")
print(head(lncRNA_annotation, 5))

# ------------------------------------------------------------
# 6. Match to the Seurat object
# ------------------------------------------------------------
seu <- readRDS(file.path(out_dir, "seurat_integrated_annotated.rds"))
message("\nSeurat object: ", ncol(seu), " cells, ", nrow(seu), " genes")

genes_in_data <- rownames(seu)
detected_lncRNAs <- intersect(lncRNA_annotation$gene_symbol, genes_in_data)
message("\nlncRNAs detected in the data: ", length(detected_lncRNAs),
        " of ", nrow(lncRNA_annotation), " annotated")

lncRNA_annotation$detected_in_data <- lncRNA_annotation$gene_symbol %in% genes_in_data

detection_by_biotype <- lncRNA_annotation %>%
  group_by(biotype) %>%
  summarize(
    n_annotated  = n(),
    n_detected   = sum(detected_in_data),
    pct_detected = round(100 * n_detected / n_annotated, 1),
    .groups = "drop"
  )
message("\nDetection rate by biotype:")
print(detection_by_biotype)

# ------------------------------------------------------------
# 7. Per-stage expression statistics
# ------------------------------------------------------------
DefaultAssay(seu) <- "RNA"
seu$stage <- factor(seu$stage, levels = c("MBC", "prePB", "PB", "PC"))

expr_data <- GetAssayData(seu, layer = "data")

if (length(detected_lncRNAs) > 0) {
  lncRNA_expr <- expr_data[detected_lncRNAs, , drop = FALSE]

  pct_expressing <- Matrix::rowMeans(lncRNA_expr > 0) * 100
  mean_expr      <- Matrix::rowMeans(lncRNA_expr)

  stage_means <- sapply(levels(seu$stage), function(s) {
    cells_in_stage <- colnames(seu)[seu$stage == s]
    Matrix::rowMeans(lncRNA_expr[, cells_in_stage, drop = FALSE])
  })

  expression_stats <- data.frame(
    gene_symbol     = detected_lncRNAs,
    pct_expressing  = pct_expressing,
    mean_expression = mean_expr,
    mean_MBC        = stage_means[, "MBC"],
    mean_prePB      = stage_means[, "prePB"],
    mean_PB         = stage_means[, "PB"],
    mean_PC         = stage_means[, "PC"],
    stringsAsFactors = FALSE
  )

  lncRNA_annotation_full <- merge(
    lncRNA_annotation, expression_stats,
    by = "gene_symbol", all.x = TRUE
  )

  message("\nExpression statistics computed for ", nrow(expression_stats), " detected lncRNAs")
  message("\nTop 10 most broadly expressed lncRNAs:")
  print(
    expression_stats %>%
      arrange(desc(pct_expressing)) %>%
      slice_head(n = 10) %>%
      select(gene_symbol, pct_expressing, mean_expression)
  )
} else {
  warning("No detected lncRNAs - check gene symbol matching")
  lncRNA_annotation_full <- lncRNA_annotation
}

# ------------------------------------------------------------
# 8. Save outputs
# ------------------------------------------------------------
write.csv(lncRNA_annotation_full,
          file.path(tbl_dir, "07_lncRNA_annotation.csv"),
          row.names = FALSE)

detected_only <- lncRNA_annotation_full[lncRNA_annotation_full$detected_in_data, ]
write.csv(detected_only,
          file.path(tbl_dir, "07_lncRNA_detected.csv"),
          row.names = FALSE)

saveRDS(detected_lncRNAs,
        file.path(out_dir, "lncRNA_detected_symbols.rds"))

message("\nSaved:")
message("  results/tables/07_lncRNA_annotation.csv")
message("  results/tables/07_lncRNA_detected.csv")
message("  data/processed/lncRNA_detected_symbols.rds")

# ------------------------------------------------------------
# 9. Diagnostic figures
# ------------------------------------------------------------
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

p_detection <- ggplot(detection_by_biotype,
                      aes(x = biotype, y = pct_detected)) +
  geom_col(fill = "#b15835", width = 0.5) +
  geom_text(aes(label = paste0(n_detected, " / ", n_annotated)),
            vjust = -0.5, color = "#2a2622", size = 3.5) +
  scale_y_continuous(limits = c(0, max(detection_by_biotype$pct_detected) * 1.2)) +
  theme_portfolio() +
  labs(
    title    = "lncRNA detection rate in plasma cell scRNA-seq",
    subtitle = "GENCODE v19 annotation, mapped to BD Rhapsody hg19-aligned data",
    x = "Biotype",
    y = "% of annotated genes detected"
  )

ggsave(file.path(fig_dir, "18_lncRNA_detection.pdf"),
       p_detection, width = 7, height = 5)
ggsave(file.path(fig_dir, "18_lncRNA_detection.png"),
       p_detection, width = 7, height = 5, dpi = 200)

if (length(detected_lncRNAs) > 0) {
  p_breadth <- ggplot(expression_stats, aes(x = pct_expressing)) +
    geom_histogram(bins = 50, fill = "#7d8c6e", color = "white", linewidth = 0.2) +
    theme_portfolio() +
    labs(
      title    = "Expression breadth of detected lncRNAs",
      subtitle = paste0(length(detected_lncRNAs), " lncRNAs across ", ncol(seu), " cells"),
      x = "% of cells expressing the lncRNA",
      y = "Number of lncRNAs"
    )

  ggsave(file.path(fig_dir, "19_lncRNA_expression_breadth.pdf"),
         p_breadth, width = 7, height = 5)
  ggsave(file.path(fig_dir, "19_lncRNA_expression_breadth.png"),
         p_breadth, width = 7, height = 5, dpi = 200)
}

message("\n========================================")
message("Script 08 (lncRNA annotation) complete!")
message("\nNumbers:")
message("  Annotated lncRNAs:    ", nrow(lncRNA_annotation))
message("  Detected in data:     ", length(detected_lncRNAs))
message("  Detection rate (all): ",
        round(100 * length(detected_lncRNAs) / nrow(lncRNA_annotation), 1), "%")
message("\nFigures:")
message("  18_lncRNA_detection")
message("  19_lncRNA_expression_breadth")
message("\nNext: 09_lncRNA_differential.R")
message("========================================")
