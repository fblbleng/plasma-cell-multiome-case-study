# ============================================================
# fix_joinlayers.R   (one-off utility)
#
# In Seurat v5, merged objects have split counts layers
# (counts.NPCD_rep1, counts.NPCD_rep2). FindAllMarkers and
# similar v5 tools require these consolidated into one layer.
# This script rewrites all saved Seurat objects with joined
# layers so downstream scripts work without modification.
#
# Run this ONCE, then re-run scripts 07 (and any script 06
# that didn't have the inline fix).
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
})

out_dir <- "data/processed"

# All Seurat objects that might have split layers
files <- c(
  "seurat_with_cellcycle.rds",
  "seurat_A_full_regression.rds",
  "seurat_B_CCdiff_regression.rds",
  "seurat_C_no_regression.rds",
  "seurat_integrated_clustered.rds"
)

for (f in files) {
  fp <- file.path(out_dir, f)
  if (!file.exists(fp)) {
    message("Skipping (does not exist): ", f)
    next
  }
  message("\nProcessing: ", f)
  obj <- readRDS(fp)

  if ("RNA" %in% Assays(obj)) {
    DefaultAssay(obj) <- "RNA"
    layers_before <- Layers(obj, assay = "RNA")
    message("  Layers before: ", paste(layers_before, collapse = ", "))

    # Only join if there's more than one counts layer
    if (sum(grepl("^counts", layers_before)) > 1 ||
        sum(grepl("^data", layers_before)) > 1) {
      obj <- JoinLayers(obj, assay = "RNA")
      layers_after <- Layers(obj, assay = "RNA")
      message("  Layers after:  ", paste(layers_after, collapse = ", "))
      saveRDS(obj, fp)
      message("  Saved.")
    } else {
      message("  Already joined - skipping save.")
    }
  } else {
    message("  No RNA assay found, skipping.")
  }
}

message("\n=== JoinLayers fix complete ===")
message("You can now re-run scripts 06 and 07 without errors.")
