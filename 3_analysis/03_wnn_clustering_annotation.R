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

# 03_wnn_clustering_annotation.R
# Annotate cell types using marker genes or existing annotations

seu <- readRDS(file.path(int_dir, "seurat_objects", "tb_scrna_integrated.rds"))
message("Loaded object with ", ncol(seu), " cells")

# Check if cell type annotations already exist in metadata
cell_type_cols <- c("cell_type", "celltype", "CellType", "cell.type", "cluster_annotation", "predicted.celltype")
existing_ct <- intersect(cell_type_cols, colnames(seu@meta.data))

if (length(existing_ct) > 0) {
  # Use existing cell type annotation
  seu$cell_type <- seu@meta.data[[existing_ct[1]]]
  message("Using existing cell type annotation from column: ", existing_ct[1])
} else {
  # Annotate based on canonical markers (for macrophage/monocyte enriched BAL data)
  message("No existing cell type annotation found. Using marker-based annotation...")
  
  # For BAL macrophage data, identify major subsets
  markers <- list(
    "Alveolar_Macrophage" = c("FABP4", "MARCO", "MCEMP1"),
    "Interstitial_Macrophage" = c("C1QA", "C1QB", "APOE"),
    "Monocyte" = c("CD14", "S100A8", "S100A9", "VCAN"),
    "T_cell" = c("CD3D", "CD3E", "CD3G"),
    "B_cell" = c("CD79A", "MS4A1"),
    "NK_cell" = c("GNLY", "NKG7", "KLRD1"),
    "Dendritic_cell" = c("CD1C", "FCER1A", "CLEC10A"),
    "Neutrophil" = c("FCGR3B", "CSF3R", "CXCR2")
  )
  
  # Score each cell type
  for (ct in names(markers)) {
    valid_markers <- intersect(markers[[ct]], rownames(seu))
    if (length(valid_markers) > 0) {
      seu <- AddModuleScore(seu, features = list(valid_markers), name = paste0("score_", ct))
    }
  }
  
  # Assign cell type based on highest score
  score_cols <- grep("^score_", colnames(seu@meta.data), value = TRUE)
  if (length(score_cols) > 0) {
    score_mat <- as.matrix(seu@meta.data[, score_cols])
    seu$cell_type <- gsub("^score_|1$", "", score_cols[apply(score_mat, 1, which.max)])
  } else {
    seu$cell_type <- paste0("Cluster_", Idents(seu))
  }
}

# Summary of cell types
ct_summary <- as.data.frame(table(seu$cell_type))
colnames(ct_summary) <- c("cell_type", "n_cells")
ct_summary$proportion <- ct_summary$n_cells / sum(ct_summary$n_cells)
write_table(ct_summary, file.path(res_dir, "tables", "Table_celltype_summary.csv"))

saveRDS(seu, file.path(int_dir, "seurat_objects", "tb_scrna_annotated.rds"))
message("Saved annotated object")

# Generate cell type UMAP
umap_name <- if ("umap" %in% names(seu@reductions)) "umap" else "UMAP"
p <- DimPlot(seu, reduction = umap_name, group.by = "cell_type") + 
  ggtitle("TB scRNA UMAP (cell types)")
ggsave(file.path(res_dir, "figures", "Fig2A_umap_celltypes.png"), p, width=9, height=6, dpi=300)
message("Cell type UMAP saved.")
