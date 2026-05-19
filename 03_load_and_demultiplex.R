# ============================================================
# 03_load_and_demultiplex.R   (REVISED)
# Load BD Rhapsody scRNA-seq matrices (GSE242330) and assign
# each cell to its differentiation stage using the per-cell
# Sample_Tag metadata from GSE242330_SampleTag.xlsx.
#
# Excel file structure:
#   - Sheet 1: Legend (tag -> stage mapping)
#   - Sheet 2: Sample_tag_rep1 (Cell_Index, Sample_Tag, Sample_Name)
#   - Sheet 3: Sample_tag_rep2 (Cell_Index, Sample_Tag, Sample_Name)
#
# Input files in data/raw/GSE242330/:
#   - GSM7758185_NPCD_rep1_RSEC_ReadsPerCell.csv.gz
#   - GSM7758186_NPCD_rep2_RSEC_ReadsPerCell.csv.gz
#   - GSE242330_SampleTag.xlsx
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(Matrix)
  library(ggplot2)
  library(patchwork)
  library(readxl)
})

set.seed(42)
raw_dir <- "data/raw/GSE242330"
out_dir <- "data/processed"
fig_dir <- "results/figures"
tbl_dir <- "results/tables"
for (d in c(out_dir, fig_dir, tbl_dir)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------
# 1. Read sample tag metadata from Excel
# ------------------------------------------------------------
xlsx_path <- file.path(raw_dir, "GSE242330_SampleTag.xlsx")
stopifnot(file.exists(xlsx_path))

message("\n--- Sheets in the Excel file ---")
sheets <- excel_sheets(xlsx_path)
print(sheets)

# The metadata sheets have a multi-line header block before the actual table.
# Find the row where "Cell_Index" appears and read from there.
read_tag_sheet <- function(path, sheet) {
  raw <- read_excel(path, sheet = sheet, col_names = FALSE,
                    .name_repair = "minimal")
  header_row <- which(raw[[1]] == "Cell_Index")[1]
  if (is.na(header_row)) stop("Couldn't find 'Cell_Index' row in sheet: ", sheet)
  message("Sheet '", sheet, "': data starts at row ", header_row)

  df <- read_excel(path, sheet = sheet, skip = header_row - 1)
  df <- as.data.frame(df)
  df$Cell_Index <- as.integer(df$Cell_Index)
  return(df)
}

# Find the per-replicate sheets (names may vary slightly)
rep1_meta_sheet <- grep("rep1", sheets, value = TRUE, ignore.case = TRUE)[1]
rep2_meta_sheet <- grep("rep2", sheets, value = TRUE, ignore.case = TRUE)[1]
message("Using sheet for rep1: ", rep1_meta_sheet)
message("Using sheet for rep2: ", rep2_meta_sheet)

rep1_meta <- read_tag_sheet(xlsx_path, rep1_meta_sheet)
rep2_meta <- read_tag_sheet(xlsx_path, rep2_meta_sheet)

message("\n--- Rep1 metadata: cells per Sample_Name ---")
print(table(rep1_meta$Sample_Name, useNA = "ifany"))
message("\n--- Rep2 metadata: cells per Sample_Name ---")
print(table(rep2_meta$Sample_Name, useNA = "ifany"))

# ------------------------------------------------------------
# 2. BD Rhapsody RSEC matrix loader
# ------------------------------------------------------------
read_bd_rsec <- function(filepath) {
  message("\n--- Reading ", basename(filepath), " ---")
  con <- gzfile(filepath, "r")
  preview <- readLines(con, n = 30)
  close(con)
  n_skip <- sum(grepl("^#", preview))

  df <- read.csv(filepath, skip = n_skip, header = TRUE,
                 stringsAsFactors = FALSE, check.names = FALSE)
  df$Cell_Index <- as.integer(df$Cell_Index)
  message("Dim: ", nrow(df), " cells x ", ncol(df) - 1, " genes")
  return(df)
}

rep1_df <- read_bd_rsec(file.path(raw_dir, "GSM7758185_NPCD_rep1_RSEC_ReadsPerCell.csv.gz"))
rep2_df <- read_bd_rsec(file.path(raw_dir, "GSM7758186_NPCD_rep2_RSEC_ReadsPerCell.csv.gz"))

# ------------------------------------------------------------
# 3. Join metadata onto count matrices
# ------------------------------------------------------------
attach_meta <- function(df, meta, replicate_num) {
  before <- nrow(df)
  df <- df %>%
    left_join(meta[, c("Cell_Index", "Sample_Tag", "Sample_Name")],
              by = "Cell_Index")

  n_matched   <- sum(!is.na(df$Sample_Name))
  n_unmatched <- sum(is.na(df$Sample_Name))
  message("\nRep ", replicate_num, " join: ",
          n_matched, " matched / ", n_unmatched, " unmatched / ",
          before, " total")

  df <- df[!is.na(df$Sample_Name), ]
  message("Stage breakdown after join:")
  print(table(df$Sample_Name))
  df$replicate <- replicate_num
  return(df)
}

rep1_df <- attach_meta(rep1_df, rep1_meta, 1)
rep2_df <- attach_meta(rep2_df, rep2_meta, 2)

# ------------------------------------------------------------
# 4. Drop Multiplets and Undetermined
# ------------------------------------------------------------
valid_stages <- c("MBC", "prePB", "PB", "PC")

clean_stages <- function(df, replicate_num) {
  before <- nrow(df)
  df <- df[df$Sample_Name %in% valid_stages, ]
  after <- nrow(df)
  message("Rep ", replicate_num, ": dropped ", before - after,
          " Multiplets/Undetermined; kept ", after, " cells")
  return(df)
}

rep1_df <- clean_stages(rep1_df, 1)
rep2_df <- clean_stages(rep2_df, 2)

# ------------------------------------------------------------
# 5. Build Seurat objects
# ------------------------------------------------------------
df_to_seurat <- function(df, project_name) {
  meta_cols <- c("Cell_Index", "Sample_Tag", "Sample_Name", "replicate")
  meta <- df[, meta_cols, drop = FALSE]
  colnames(meta)[colnames(meta) == "Sample_Name"] <- "stage"

  # Unique cell names: stage_repN_cellindex
  cell_names <- paste0(meta$stage, "_rep", meta$replicate, "_", meta$Cell_Index)
  rownames(meta) <- cell_names

  # Build counts matrix (genes x cells)
  count_cols <- setdiff(colnames(df), meta_cols)
  counts <- as.matrix(df[, count_cols])
  rownames(counts) <- cell_names
  counts_t <- t(counts)
  counts_t <- Matrix::Matrix(counts_t, sparse = TRUE)

  meta$stage <- factor(meta$stage, levels = c("MBC", "prePB", "PB", "PC"))
  meta$replicate <- factor(meta$replicate)

  obj <- CreateSeuratObject(
    counts = counts_t,
    project = project_name,
    meta.data = meta[, c("stage", "replicate", "Sample_Tag", "Cell_Index")],
    min.cells = 3,
    min.features = 200
  )

  message("\n", project_name, " Seurat object:")
  print(obj)
  message("Cells per stage:")
  print(table(obj$stage))
  return(obj)
}

seu_rep1 <- df_to_seurat(rep1_df, "NPCD_rep1")
seu_rep2 <- df_to_seurat(rep2_df, "NPCD_rep2")

# ------------------------------------------------------------
# 6. Save objects
# ------------------------------------------------------------
saveRDS(seu_rep1, file.path(out_dir, "seurat_rep1_raw.rds"))
saveRDS(seu_rep2, file.path(out_dir, "seurat_rep2_raw.rds"))
message("\nSaved raw Seurat objects to ", out_dir)

# ------------------------------------------------------------
# 7. First case study figure: cells per stage per replicate
# ------------------------------------------------------------
counts_df <- rbind(
  data.frame(replicate = "rep1", stage = seu_rep1$stage),
  data.frame(replicate = "rep2", stage = seu_rep2$stage)
)

p <- ggplot(counts_df, aes(x = stage, fill = replicate)) +
  geom_bar(position = position_dodge(width = 0.8), width = 0.7,
           color = "white", linewidth = 0.3) +
  geom_text(stat = "count", aes(label = after_stat(count)),
            position = position_dodge(width = 0.8),
            vjust = -0.4, size = 3.2, color = "#5c544b") +
  scale_fill_manual(values = c("rep1" = "#b15835", "rep2" = "#7d8c6e"),
                    name = "Replicate") +
  theme_minimal(base_size = 13) +
  labs(
    title = "Cells recovered per differentiation stage",
    subtitle = "GSE242330 - Moreaux et al., Blood 2024 | BD Rhapsody WTA",
    x = NULL,
    y = "Number of cells",
    caption = "After sample tag demultiplexing; Multiplets and Undetermined excluded"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 14, color = "#2a2622"),
    plot.subtitle = element_text(size = 11, color = "#5c544b"),
    plot.caption = element_text(size = 9, color = "#8a7f73", hjust = 0),
    axis.text = element_text(color = "#2a2622"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position = "top"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12)))

ggsave(file.path(fig_dir, "01_cells_per_stage.pdf"), p, width = 7, height = 4.5)
ggsave(file.path(fig_dir, "01_cells_per_stage.png"), p, width = 7, height = 4.5, dpi = 200)
print(p)

write.csv(as.data.frame.matrix(table(counts_df$stage, counts_df$replicate)),
          file.path(tbl_dir, "01_cells_per_stage.csv"))

message("\n========================================")
message("Stage 2.1 complete!")
message("- Cells in rep1: ", ncol(seu_rep1))
message("- Cells in rep2: ", ncol(seu_rep2))
message("- Figure: ", fig_dir, "/01_cells_per_stage.pdf")
message("Next: 04_qc_and_filter.R")
message("========================================")
