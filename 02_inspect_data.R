# ============================================================
# 02_inspect_data.R
# Characterize what we actually downloaded so we know
# how to proceed in Stage 2 (scRNA-seq analysis).
# ============================================================

library(R.utils)

raw_dir <- "data/raw"
scrna_acc <- "GSE242330"
scatac_acc <- "GSE242324"

# ---- Helper: report file info ----
report_files <- function(path) {
  files <- list.files(path, recursive = TRUE, full.names = TRUE)
  if (length(files) == 0) {
    message("  (no files)")
    return(invisible(NULL))
  }
  info <- file.info(files)
  df <- data.frame(
    file = basename(files),
    size_MB = round(info$size / 1024^2, 2),
    stringsAsFactors = FALSE
  )
  print(df[order(-df$size_MB), ])
}

message("\n========== scRNA-seq (", scrna_acc, ") ==========")
report_files(file.path(raw_dir, scrna_acc))

message("\n========== scATAC-seq (", scatac_acc, ") ==========")
report_files(file.path(raw_dir, scatac_acc))

# ============================================================
# Peek inside any tab-delimited / matrix files
# This helps you understand the format:
#   - Are these dense count matrices?  (genes × cells)
#   - Cell Ranger triplets (matrix.mtx, barcodes.tsv, features.tsv)?
#   - Seurat .rds objects?
#   - Fragment files (.tsv.gz) for ATAC?
# ============================================================

peek_file <- function(filepath, n = 5) {
  message("\n--- ", basename(filepath), " ---")
  ext <- tools::file_ext(filepath)

  tryCatch({
    if (grepl("\\.gz$", filepath) && !grepl("\\.tar\\.gz$", filepath)) {
      # gzipped text — read first few lines
      con <- gzfile(filepath, "r")
      lines <- readLines(con, n = n)
      close(con)
      cat(lines, sep = "\n")
    } else if (ext %in% c("txt", "csv", "tsv")) {
      lines <- readLines(filepath, n = n)
      cat(lines, sep = "\n")
    } else if (ext == "rds") {
      obj <- readRDS(filepath)
      cat("Class: ", class(obj)[1], "\n")
      print(obj)
    } else if (ext %in% c("h5", "hdf5")) {
      cat("HDF5 file — inspect with rhdf5::h5ls() or Seurat::Read10X_h5()\n")
    } else if (ext == "tar") {
      cat("Tar archive — list contents with: untar('", filepath, "', list=TRUE)\n", sep = "")
    } else {
      cat("Unknown format — manual inspection needed.\n")
    }
  }, error = function(e) {
    message("Could not peek: ", e$message)
  })
}

# Auto-peek every file
all_files <- c(
  list.files(file.path(raw_dir, scrna_acc),  recursive = TRUE, full.names = TRUE),
  list.files(file.path(raw_dir, scatac_acc), recursive = TRUE, full.names = TRUE)
)

for (f in all_files) {
  peek_file(f)
}

message("\n=== Inspection complete ===")
message("Based on what you see above, we'll choose the right loader:")
message("  - matrix.mtx + barcodes.tsv + features.tsv  -> Seurat::Read10X()")
message("  - filtered_feature_bc_matrix.h5             -> Seurat::Read10X_h5()")
message("  - .rds file                                 -> readRDS()")
message("  - dense .txt/.csv matrix                    -> read.table() then CreateSeuratObject()")
