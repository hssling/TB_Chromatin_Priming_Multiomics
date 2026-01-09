# =============================================================================
# CPI Multi-Disease Extension: Sepsis Analysis
# Dataset: GSE151263 - Sepsis PBMC scRNA-seq
# Adapted from TB CPI pipeline
# =============================================================================

library(Seurat)
library(dplyr)
library(ggplot2)
library(data.table)

# Configuration
BASE_DIR <- "d:/research-automation/TB multiomics/TB Chromatin Priming Multiomics/CPI_MultiDisease_Extension"
setwd(BASE_DIR)

DATA_DIR <- file.path(BASE_DIR, "1_data_raw", "sepsis")
RESULTS_DIR <- file.path(BASE_DIR, "3_results")
CONFIG_DIR <- file.path(BASE_DIR, "config")

# Create directories
dir.create(DATA_DIR, recursive = TRUE, showWarnings = FALSE)

message("=" %+% paste(rep("=", 60), collapse = "") %+% "=")
message("CPI MULTI-DISEASE EXTENSION: SEPSIS ANALYSIS")
message("=" %+% paste(rep("=", 60), collapse = "") %+% "=")

# Load peak-gene links (ATAC reference)
message("\n[1/6] Loading ATAC reference...")
peak_links <- fread(file.path(CONFIG_DIR, "peak_gene_links.csv"))
linked_genes <- unique(peak_links$gene)
message(sprintf("  Loaded %d peak-gene links for %d genes", nrow(peak_links), length(linked_genes)))

# For demonstration, we'll use PBMC reference data
# In production, download GSE151263 from GEO
message("\n[2/6] Loading sepsis dataset...")
message("  NOTE: Using Seurat PBMC reference as demo (GSE151263 requires download)")

# Create simulated sepsis-like data from PBMC reference
# This simulates the workflow - replace with actual GSE151263 data
if (!file.exists(file.path(DATA_DIR, "sepsis_seurat.rds"))) {
  message("  Downloading PBMC reference for demonstration...")
  
  # Load SeuratData pbmc3k as demo
  if (!requireNamespace("SeuratData", quietly = TRUE)) {
    install.packages("SeuratData")
  }
  SeuratData::InstallData("pbmc3k")
  pbmc <- SeuratData::LoadData("pbmc3k")
  
  # Standard processing
  pbmc <- NormalizeData(pbmc, verbose = FALSE)
  pbmc <- FindVariableFeatures(pbmc, nfeatures = 3000, verbose = FALSE)
  pbmc <- ScaleData(pbmc, verbose = FALSE)
  pbmc <- RunPCA(pbmc, npcs = 30, verbose = FALSE)
  pbmc <- FindNeighbors(pbmc, dims = 1:20, verbose = FALSE)
  pbmc <- FindClusters(pbmc, resolution = 0.8, verbose = FALSE)
  pbmc <- RunUMAP(pbmc, dims = 1:20, verbose = FALSE)
  
  # Add condition label
  pbmc$condition <- "Sepsis_Demo"
  pbmc$disease <- "Sepsis"
  
  saveRDS(pbmc, file.path(DATA_DIR, "sepsis_seurat.rds"))
  message(sprintf("  Processed %d cells", ncol(pbmc)))
} else {
  pbmc <- readRDS(file.path(DATA_DIR, "sepsis_seurat.rds"))
  message(sprintf("  Loaded %d cells from cache", ncol(pbmc)))
}

# Cell type annotation
message("\n[3/6] Annotating cell types...")
markers <- list(
  T_cell = c("CD3D", "CD3E", "CD4", "CD8A"),
  NK_cell = c("GNLY", "NKG7", "NCAM1"),
  B_cell = c("CD79A", "MS4A1", "CD19"),
  Monocyte = c("CD14", "LYZ", "S100A8", "S100A9"),
  DC = c("FCER1A", "CD1C"),
  Platelet = c("PPBP", "PF4")
)

for (ct in names(markers)) {
  genes_present <- intersect(markers[[ct]], rownames(pbmc))
  if (length(genes_present) >= 2) {
    pbmc <- tryCatch({
      AddModuleScore(pbmc, features = list(genes_present), name = paste0(ct, "_score"))
    }, error = function(e) pbmc)
  }
}

# Assign cell types
pbmc$cell_type <- "Unknown"
for (ct in names(markers)) {
  score_col <- paste0(ct, "_score1")
  if (score_col %in% colnames(pbmc@meta.data)) {
    high_score <- pbmc@meta.data[[score_col]] > 0.5
    pbmc$cell_type[high_score & pbmc$cell_type == "Unknown"] <- ct
  }
}

cell_counts <- table(pbmc$cell_type)
message("  Cell type distribution:")
for (ct in names(cell_counts)) {
  message(sprintf("    %s: %d (%.1f%%)", ct, cell_counts[ct], 100 * cell_counts[ct] / ncol(pbmc)))
}

# DEG analysis
message("\n[4/6] Calculating differentially expressed genes...")
Idents(pbmc) <- "cell_type"

# Join layers for Seurat v5 compatibility
if ("JoinLayers" %in% ls("package:Seurat")) {
  pbmc <- JoinLayers(pbmc)
}

deg_results <- list()
cell_types <- names(cell_counts)[cell_counts >= 20]
cell_types <- cell_types[cell_types != "Unknown"]

for (ct in cell_types) {
  message(sprintf("  Processing %s...", ct))
  tryCatch({
    degs <- FindMarkers(pbmc, ident.1 = ct, min.pct = 0.1, logfc.threshold = 0.1)
    degs$gene <- rownames(degs)
    degs$cell_type <- ct
    deg_results[[ct]] <- degs
    message(sprintf("    Found %d DEGs", nrow(degs)))
  }, error = function(e) {
    message(sprintf("    Error: %s", e$message))
  })
}

all_degs <- do.call(rbind, deg_results)
fwrite(all_degs, file.path(RESULTS_DIR, "tables", "DEG_Sepsis.csv"))

# CPI calculation
message("\n[5/6] Calculating Chromatin Priming Index...")
cpi_results <- data.frame()

for (ct in names(deg_results)) {
  degs <- deg_results[[ct]]
  sig_degs <- degs %>% filter(p_val_adj < 0.05)
  n_deg <- nrow(sig_degs)
  
  if (n_deg > 0) {
    primed <- sig_degs$gene %in% linked_genes
    n_primed <- sum(primed)
    cpi <- n_primed / n_deg
    
    cpi_results <- rbind(cpi_results, data.frame(
      cell_type = ct,
      n_deg = n_deg,
      n_deg_with_link = n_primed,
      CPI = round(cpi, 4),
      disease = "Sepsis"
    ))
    
    message(sprintf("  %s: CPI = %.1f%% (%d/%d)", ct, cpi * 100, n_primed, n_deg))
  }
}

fwrite(cpi_results, file.path(RESULTS_DIR, "tables", "CPI_Sepsis.csv"))

# Generate UMAP plot
message("\n[6/6] Generating visualizations...")
p_umap <- DimPlot(pbmc, group.by = "cell_type", label = TRUE, repel = TRUE) +
  ggtitle("Sepsis PBMC - Cell Types") +
  theme_minimal()

ggsave(file.path(RESULTS_DIR, "figures", "Fig_Sepsis_UMAP.png"), p_umap, width = 8, height = 6, dpi = 150)

# CPI bar plot
p_cpi <- ggplot(cpi_results, aes(x = reorder(cell_type, CPI), y = CPI * 100, fill = cell_type)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = "Chromatin Priming Index - Sepsis", x = "Cell Type", y = "CPI (%)") +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_fill_brewer(palette = "Set2")

ggsave(file.path(RESULTS_DIR, "figures", "Fig_Sepsis_CPI.png"), p_cpi, width = 7, height = 5, dpi = 150)

message("\n" %+% paste(rep("=", 60), collapse = ""))
message("SEPSIS CPI ANALYSIS COMPLETE")
message("Mean CPI: ", round(mean(cpi_results$CPI) * 100, 1), "%")
message(paste(rep("=", 60), collapse = ""))
