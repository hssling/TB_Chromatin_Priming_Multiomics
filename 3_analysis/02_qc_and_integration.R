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

# 02_qc_and_integration.R
# The GSE167232 data is already QC'd and integrated, so we just need to ensure standard processing

seu <- readRDS(file.path(int_dir, "seurat_objects", "tb_scrna_raw.rds"))
message("Loaded object with ", ncol(seu), " cells")

# Check if already normalized and processed
if (!"nCount_RNA" %in% colnames(seu@meta.data)) {
  seu$nCount_RNA <- Matrix::colSums(seu@assays$RNA@counts)
  seu$nFeature_RNA <- Matrix::colSums(seu@assays$RNA@counts > 0)
}

# Ensure standard processing is done
if (!"PCA" %in% names(seu@reductions)) {
  message("Running standard Seurat processing...")
  seu <- NormalizeData(seu, verbose = FALSE)
  seu <- FindVariableFeatures(seu, nfeatures = cfg$analysis$hvgs, verbose = FALSE)
  seu <- ScaleData(seu, verbose = FALSE)
  seu <- RunPCA(seu, verbose = FALSE)
}

if (!"umap" %in% names(seu@reductions) && !"UMAP" %in% names(seu@reductions)) {
  message("Running UMAP...")
  seu <- FindNeighbors(seu, dims = 1:30, verbose = FALSE)
  seu <- FindClusters(seu, resolution = cfg$analysis$clustering_resolution, verbose = FALSE)
  seu <- RunUMAP(seu, dims = 1:30, verbose = FALSE)
}

out_path <- file.path(int_dir, "seurat_objects", "tb_scrna_integrated.rds")
saveRDS(seu, out_path)
message("Saved: ", out_path)

# Generate UMAP plot
umap_name <- if ("umap" %in% names(seu@reductions)) "umap" else "UMAP"
p <- DimPlot(seu, reduction = umap_name, group.by = "seurat_clusters") + 
  ggtitle("TB scRNA UMAP (clusters)")
ensure_dir(file.path(res_dir, "figures"))
ggsave(file.path(res_dir, "figures", "Fig2A_umap_clusters.png"), p, width=7, height=5, dpi=300)
message("UMAP plot saved.")

# Also plot by condition if available
if ("condition" %in% colnames(seu@meta.data)) {
  p2 <- DimPlot(seu, reduction = umap_name, group.by = "condition") +
    ggtitle("TB scRNA UMAP (by condition)")
  ggsave(file.path(res_dir, "figures", "Fig2A_umap_condition.png"), p2, width=7, height=5, dpi=300)
}
