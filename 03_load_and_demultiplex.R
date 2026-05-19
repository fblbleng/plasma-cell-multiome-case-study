# ============================================================
# 03_load_and_demultiplex.R
# Load BD Rhapsody scRNA-seq matrices (GSE242330) and assign
# each cell to its differentiation stage + replicate based on
# the CellID prefix.
#
# Expected input files (in data/raw/GSE242330/):
#   - GSM7758185_NPCD_rep1_RSEC_ReadsPerCell.csv.gz
#   - GSM7758186_NPCD_rep2_RSEC_ReadsPerCell.csv.gz
#
# Sample tag → stage mapping (from GSE242330_SampleTag.xlsx):
#   Rep 1: 1_1_=MBC | 2_3_=prePB | 3_5_=PB | 4_7_=PC
#   Rep 2: 1_2_=MBC | 2_4_=prePB | 3_6_=PB | 4_8_=PC
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(Matrix)
  library(ggplot2)
  library(patchwork)
})

set.seed(42)
raw_dir <- "data/raw/GSE242330"
out_dir <- "data/processed"
fig_dir <- "results/figures"
tbl_dir <- "results/tables"
for (d in c(out_dir, fig_dir, tbl_dir)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------
# 1. Sample tag → stage mapping
# ------------------------------------------------------------
# Each row: prefix used in the CellID column of the matrix
tag_map <- data.frame(
  prefix    = c("1_1_", "2_3_", "3_5_", "4_7_",
                "1_2_", "2_4_", "3_6_", "4_8_"),
  stage     = c("MBC",  "prePB", "PB",  "PC",
                "MBC",  "prePB", "PB",  "PC"),
  replicate = c(1, 1, 1, 1, 2, 2, 2, 2),
  stringsAsFactors = FALSE
)
tag_map$sample_id <- paste0(tag_map$stage, "_rep", tag_map$replicate)
print(tag_map)

# ------------------------------------------------------------
# 2. Helper: BD Rhapsody RSEC matrix loader
# ------------------------------------------------------------
# BD Rhapsody RSEC_ReadsPerCell.csv format (typical):
#   - Header lines starting with "#" (skip these)
#   - First column: Cell_Index (numeric cell barcode)
#   - Remaining columns: gene names with counts
#   - Rows = cells, Columns = genes (TRANSPOSED relative to Seurat!)
#
# We'll auto-detect the header lines and parse robustly.
# ------------------------------------------------------------
read_bd_rsec <- function(filepath, sample_label) {
  message("\n--- Reading ", basename(filepath), " ---")

  # Peek first lines to find where the data starts
  con <- gzfile(filepath, "r")
  preview <- readLines(con, n = 30)
  close(con)

  # Count comment lines
  n_skip <- sum(grepl("^#", preview))
  message("Skipping ", n_skip, " header/comment lines")

  # Read full table
  df <- read.csv(filepath, skip = n_skip, header = TRUE,
                 stringsAsFactors = FALSE, check.names = FALSE)
  message("Dimensions (rows x cols): ", nrow(df), " x ", ncol(df))
  message("First 3 column names: ", paste(head(colnames(df), 3), collapse = " | "))

  # The first column is typically Cell_Index — confirm
  message("First column class: ", class(df[[1]]))
  message("First 3 values of column 1: ", paste(head(df[[1]], 3), collapse = ", "))

  return(df)
}

# ------------------------------------------------------------
# 3. Load both replicates
# ------------------------------------------------------------
rep1_path <- file.path(raw_dir, "GSM7758185_NPCD_rep1_RSEC_ReadsPerCell.csv.gz")
rep2_path <- file.path(raw_dir, "GSM7758186_NPCD_rep2_RSEC_ReadsPerCell.csv.gz")

stopifnot(file.exists(rep1_path), file.exists(rep2_path))

rep1_df <- read_bd_rsec(rep1_path, "rep1")
rep2_df <- read_bd_rsec(rep2_path, "rep2")

# ------------------------------------------------------------
# 4. Inspect the CellID format
# ------------------------------------------------------------
# We expect each cell to have a CellID like "1_1_<barcode>" (MBC rep1)
# or "2_3_<barcode>" (prePB rep1), etc.
# The first column of the matrix should hold these.
# ------------------------------------------------------------
message("\n=== Cell ID format check ===")
message("\nRep1 - first 10 cell IDs:")
print(head(rep1_df[[1]], 10))
message("\nRep2 - first 10 cell IDs:")
print(head(rep2_df[[1]], 10))

# ------------------------------------------------------------
# >>>>>>>>>>>>>>>>  STOP HERE ON FIRST RUN  <<<<<<<<<<<<<<<<
# Inspect the console output above. If cell IDs LOOK LIKE
# "1_1_AAAGTC...", "2_3_GGCATC..." etc., we're good and the
# rest of the script will run correctly.
#
# If cell IDs are pure integers (e.g., 1, 2, 3, ...) or have
# a different format, we'll need to adjust. Send me the first
# 10 cell IDs from each replicate before proceeding.
# ------------------------------------------------------------
# Uncomment the line below ONLY after confirming format:
# proceed <- TRUE

if (!exists("proceed")) {
  stop("Pause: confirm CellID format above, then set 'proceed <- TRUE' and re-run.")
}

# ------------------------------------------------------------
# 5. Assign stage + replicate from CellID prefix
# ------------------------------------------------------------
assign_stage <- function(df, replicate_num) {
  cell_ids <- as.character(df[[1]])

  # Extract the "N_N_" prefix
  prefix <- sub("^(\\d+_\\d+_).*", "\\1", cell_ids)

  # Subset tag_map to this replicate
  tm <- tag_map[tag_map$replicate == replicate_num, ]
  idx <- match(prefix, tm$prefix)

  df$stage      <- tm$stage[idx]
  df$replicate  <- replicate_num
  df$sample_id  <- tm$sample_id[idx]

  # Report stats
  message("\nReplicate ", replicate_num, " — cell counts by stage:")
  print(table(df$stage, useNA = "ifany"))

  # Drop cells with no recognized prefix
  n_drop <- sum(is.na(df$stage))
  if (n_drop > 0) {
    message("Dropping ", n_drop, " cells with unrecognized prefix")
    df <- df[!is.na(df$stage), ]
  }
  return(df)
}

rep1_df <- assign_stage(rep1_df, 1)
rep2_df <- assign_stage(rep2_df, 2)

# ------------------------------------------------------------
# 6. Convert to Seurat objects (per replicate, all stages together)
# ------------------------------------------------------------
df_to_seurat <- function(df, project_name) {
  # Separate metadata from counts
  meta_cols <- c("stage", "replicate", "sample_id")
  meta <- df[, c(colnames(df)[1], meta_cols), drop = FALSE]
  colnames(meta)[1] <- "cell_barcode"

  # Build cell names like "rep1_MBC_<cell_index>"
  cell_names <- paste0("rep", df$replicate[1], "_", df$cell_barcode <- df[[1]])
  cell_names <- paste(df$sample_id, df[[1]], sep = "_")
  rownames(meta) <- cell_names

  # Counts: rows = genes, cols = cells. BD gives us cells × genes,
  # so we transpose.
  count_cols <- setdiff(colnames(df), c(colnames(df)[1], meta_cols))
  counts <- as.matrix(df[, count_cols])
  rownames(counts) <- cell_names

  counts_t <- t(counts)                       # genes × cells
  counts_t <- Matrix::Matrix(counts_t, sparse = TRUE)

  seurat_obj <- CreateSeuratObject(
    counts = counts_t,
    project = project_name,
    meta.data = meta[, meta_cols, drop = FALSE],
    min.cells = 3,
    min.features = 200
  )
  message("\nSeurat object for ", project_name, ":")
  print(seurat_obj)
  return(seurat_obj)
}

seu_rep1 <- df_to_seurat(rep1_df, "NPCD_rep1")
seu_rep2 <- df_to_seurat(rep2_df, "NPCD_rep2")

# ------------------------------------------------------------
# 7. Save objects
# ------------------------------------------------------------
saveRDS(seu_rep1, file.path(out_dir, "seurat_rep1_raw.rds"))
saveRDS(seu_rep2, file.path(out_dir, "seurat_rep2_raw.rds"))

# ------------------------------------------------------------
# 8. Quick sanity figure: cells per stage per replicate
# ------------------------------------------------------------
sanity_tbl <- rbind(
  data.frame(replicate = "rep1", stage = seu_rep1$stage),
  data.frame(replicate = "rep2", stage = seu_rep2$stage)
)
sanity_tbl$stage <- factor(sanity_tbl$stage, levels = c("MBC", "prePB", "PB", "PC"))

p <- ggplot(sanity_tbl, aes(x = stage, fill = replicate)) +
  geom_bar(position = position_dodge()) +
  scale_fill_manual(values = c("rep1" = "#b15835", "rep2" = "#7d8c6e")) +
  theme_minimal(base_size = 13) +
  labs(title = "Cells recovered per stage and replicate",
       subtitle = "GSE242330 — Moreaux et al., Blood 2024",
       x = NULL, y = "Number of cells") +
  theme(plot.title = element_text(face = "bold"))

ggsave(file.path(fig_dir, "01_cells_per_stage.pdf"), p, width = 6, height = 4)
ggsave(file.path(fig_dir, "01_cells_per_stage.png"), p, width = 6, height = 4, dpi = 200)
print(p)

write.csv(table(sanity_tbl$replicate, sanity_tbl$stage),
          file.path(tbl_dir, "01_cells_per_stage.csv"))

message("\n=== Stage 2.1 complete ===")
message("Seurat objects saved to: ", out_dir)
message("First figure saved to: ", fig_dir, "/01_cells_per_stage.pdf")
message("\nNext: run 04_qc_and_filter.R")
