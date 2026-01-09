suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(data.table)
  library(jsonlite)
})
source("3_analysis/utils/utils_io.R")
cfg <- read_config("config/config.yaml")
set.seed(cfg$analysis$seed)

raw_dir <- cfg$paths$raw
int_dir <- cfg$paths$intermediate
res_dir <- cfg$paths$results

# 01_build_objects.R
# Load GSE167232 TB BAL scRNA-seq data (Pisu et al., Nature Immunology 2021)
message("Loading GSE167232 TB BAL integrated dataset...")

# Read the pre-processed Seurat object from GEO (decompressed)
gse_path <- file.path(raw_dir, "GSE167232_bal_integrated.RDS")
if (!file.exists(gse_path)) {
  stop("GSE167232 dataset not found. Please download from GEO and decompress.")
}

seu <- readRDS(gse_path)
# Update Seurat object to current version for compatibility
seu <- UpdateSeuratObject(seu)
message("Loaded Seurat object with ", ncol(seu), " cells and ", nrow(seu), " genes")

# Inspect the object structure
print(head(seu@meta.data))
message("Available metadata columns: ", paste(colnames(seu@meta.data), collapse=", "))

# Standardize condition column if needed
# The GSE167232 dataset has TB infection status in metadata
if ("infection" %in% colnames(seu@meta.data)) {
  seu$condition <- seu$infection
} else if ("sample" %in% colnames(seu@meta.data)) {
  seu$condition <- seu$sample
} else if ("orig.ident" %in% colnames(seu@meta.data)) {
  seu$condition <- seu$orig.ident
}

# Create sample_id if not present
if (!"sample_id" %in% colnames(seu@meta.data)) {
  seu$sample_id <- seu$orig.ident
}

# Save the initial object
out_path <- file.path(int_dir, "seurat_objects", "tb_scrna_raw.rds")
ensure_dir(dirname(out_path))
saveRDS(seu, out_path)
message("Saved: ", out_path)

# Generate summary stats
summary_df <- data.frame(
  total_cells = ncol(seu),
  total_genes = nrow(seu),
  conditions = paste(unique(seu$condition), collapse=", "),
  samples = length(unique(seu$sample_id))
)
write_table(summary_df, file.path(res_dir, "tables", "Table0_data_summary.csv"))
message("Data summary saved.")
