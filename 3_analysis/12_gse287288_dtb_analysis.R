suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(data.table)
  library(Matrix)
})
source("3_analysis/utils/utils_io.R")
cfg <- read_config("config/config.yaml")
set.seed(cfg$analysis$seed)

res_dir <- cfg$paths$results

# 12_gse287288_dtb_analysis.R
# Analyze GSE287288 Disseminated TB PBMC data and calculate CPI

message("=== Loading GSE287288 DTB Samples ===")

# Load all DTB (pre-treatment) samples
data_dir <- "1_data_raw/GSE287288"
sample_ids <- sprintf("DTB_%02d", 1:7)

all_objects <- list()

for (i in 1:7) {
  sample_name <- sample_ids[i]
  gsm <- sprintf("GSM8743%d", 122 + i)
  prefix <- paste0(gsm, "_", sample_name, "_")
  
  message(sprintf("Loading %s...", sample_name))
  
  # Read 10x files
  barcodes_file <- file.path(data_dir, paste0(prefix, "barcodes.tsv.gz"))
  features_file <- file.path(data_dir, paste0(prefix, "features.tsv.gz"))
  matrix_file <- file.path(data_dir, paste0(prefix, "matrix.mtx.gz"))
  
  if (!file.exists(matrix_file)) {
    message(sprintf("  Skipping %s - matrix file not found", sample_name))
    next
  }
  
  # Read matrix
  mat <- readMM(matrix_file)
  mat <- as(mat, "CsparseMatrix")  # Convert to dgCMatrix
  barcodes <- read.table(barcodes_file, header = FALSE, stringsAsFactors = FALSE)$V1
  features <- read.table(features_file, header = FALSE, stringsAsFactors = FALSE, sep = "\t")
  
  # Make gene names unique to avoid duplicate rownames error
  gene_names <- make.unique(features$V2)
  rownames(mat) <- gene_names
  colnames(mat) <- paste0(sample_name, "_", barcodes)
  
  # Create Seurat object
  seu <- CreateSeuratObject(counts = mat, project = sample_name, min.cells = 3, min.features = 200)
  seu$sample_id <- sample_name
  seu$patient_id <- sub("DTB_", "", sample_name)
  seu$condition <- "DTB"
  
  # Downsample to max 3000 cells to avoid memory issues
  if (ncol(seu) > 3000) {
    set.seed(42)
    cells_to_keep <- sample(colnames(seu), 3000)
    seu <- subset(seu, cells = cells_to_keep)
    message(sprintf("    Downsampled to %d cells", ncol(seu)))
  }
  
  all_objects[[sample_name]] <- seu
  message(sprintf("  %s: %d cells, %d genes", sample_name, ncol(seu), nrow(seu)))
}

# Merge all samples
message("\nMerging samples...")
seu_merged <- merge(all_objects[[1]], y = all_objects[-1], project = "GSE287288_DTB")
message(sprintf("Total: %d cells, %d genes", ncol(seu_merged), nrow(seu_merged)))

# Standard preprocessing
message("\nPerforming standard preprocessing...")
seu_merged <- NormalizeData(seu_merged, verbose = FALSE)
seu_merged <- FindVariableFeatures(seu_merged, selection.method = "vst", nfeatures = 3000, verbose = FALSE)
seu_merged <- ScaleData(seu_merged, verbose = FALSE)
seu_merged <- RunPCA(seu_merged, npcs = 30, verbose = FALSE)
seu_merged <- FindNeighbors(seu_merged, dims = 1:20, verbose = FALSE)
seu_merged <- FindClusters(seu_merged, resolution = 0.8, verbose = FALSE)
seu_merged <- RunUMAP(seu_merged, dims = 1:20, verbose = FALSE)

# Cell type annotation using markers
message("\nAnnotating cell types...")
markers <- list(
  "T_cell" = c("CD3D", "CD3E", "CD4", "CD8A"),
  "NK_cell" = c("GNLY", "NKG7", "NCAM1"),
  "B_cell" = c("CD79A", "MS4A1", "CD19"),
  "Monocyte" = c("CD14", "LYZ", "S100A8", "S100A9"),
  "DC" = c("CD1C", "FCER1A", "CLEC4C"),
  "Platelet" = c("PPBP", "PF4")
)

for (ct in names(markers)) {
  genes <- intersect(markers[[ct]], rownames(seu_merged))
  if (length(genes) >= 2) {  # Require at least 2 genes
    tryCatch({
      seu_merged <- AddModuleScore(seu_merged, features = list(genes), name = paste0(ct, "_score"), nbin = 24)
    }, error = function(e) {
      message(sprintf("  Skipping %s markers: %s", ct, e$message))
    })
  } else {
    message(sprintf("  Skipping %s: not enough markers found (%d)", ct, length(genes)))
  }
}

# Assign cell types based on highest score
score_cols <- paste0(names(markers), "_score1")
score_cols <- intersect(score_cols, colnames(seu_merged@meta.data))

if (length(score_cols) > 0) {
  scores <- seu_merged@meta.data[, score_cols]
  max_score <- apply(scores, 1, max)
  max_idx <- apply(scores, 1, which.max)
  cell_types <- gsub("_score1", "", score_cols)[max_idx]
  cell_types[max_score < 0] <- "Unknown"
  seu_merged$cell_type <- cell_types
} else {
  seu_merged$cell_type <- "PBMC"
}

message("Cell type distribution:")
print(table(seu_merged$cell_type))

# Save processed object
saveRDS(seu_merged, file.path(cfg$paths$intermediate, "seu_gse287288_dtb.rds"))

# Generate UMAP plot
p_umap <- DimPlot(seu_merged, group.by = "cell_type", label = TRUE) +
  ggtitle("GSE287288: Disseminated TB PBMC Cell Types")
ggsave(file.path(res_dir, "figures", "Fig7_DTB_umap_celltypes.png"), p_umap, width = 8, height = 6, dpi = 300)

# DEG analysis per cell type
message("\n=== Differential Expression Analysis (DTB) ===")
peak_links <- fread("2_data_intermediate/atac_reference/peak_gene_links.csv")
linked_genes <- unique(peak_links$gene)

# Join layers for Seurat v5 compatibility
seu_merged <- JoinLayers(seu_merged)

Idents(seu_merged) <- "cell_type"
all_deg <- list()
cpi_results <- list()

for (ct in unique(seu_merged$cell_type)) {
  if (ct == "Unknown") next
  
  cells_ct <- WhichCells(seu_merged, idents = ct)
  if (length(cells_ct) < 50) {
    message(sprintf("  %s: Too few cells (%d), skipping", ct, length(cells_ct)))
    next
  }
  
  message(sprintf("  Analyzing %s (%d cells)...", ct, length(cells_ct)))
  
  # Find markers for this cell type vs all others
  markers_ct <- FindMarkers(seu_merged, ident.1 = ct, min.pct = 0.1, logfc.threshold = 0.1, verbose = FALSE)
  markers_ct$gene <- rownames(markers_ct)
  markers_ct$cell_type <- ct
  
  all_deg[[ct]] <- markers_ct
  
  # Calculate CPI
  sig_deg <- markers_ct[markers_ct$p_val_adj < 0.05, ]$gene
  deg_with_link <- intersect(sig_deg, linked_genes)
  cpi <- length(deg_with_link) / max(1, length(sig_deg))
  
  cpi_results[[ct]] <- data.frame(
    cell_type = ct,
    n_deg = length(sig_deg),
    n_deg_with_link = length(deg_with_link),
    CPI = cpi
  )
  
  message(sprintf("    %d significant DEGs, CPI = %.1f%%", length(sig_deg), cpi * 100))
}

# Combine results
deg_combined <- rbindlist(all_deg, fill = TRUE)
cpi_combined <- rbindlist(cpi_results)

write_table(deg_combined, file.path(res_dir, "tables", "DEG_GSE287288_DTB.csv"))
write_table(cpi_combined, file.path(res_dir, "tables", "Table_CPI_GSE287288_DTB.csv"))

# Generate CPI comparison plot
message("\n=== Generating CPI Comparison ===")

# Load previous CPI data (BAL)
cpi_bal <- fread(file.path(res_dir, "tables", "Table_CPI_by_celltype.csv"))
cpi_bal$dataset <- "GSE167232_BAL"
cpi_combined$dataset <- "GSE287288_DTB_PBMC"

# Combine for plotting
cpi_all <- rbind(cpi_bal[, .(cell_type, CPI, dataset)], cpi_combined[, .(cell_type, CPI, dataset)])

p_cpi <- ggplot(cpi_all, aes(x = reorder(cell_type, -CPI), y = CPI * 100, fill = dataset)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Chromatin Priming Index: BAL (GSE167232) vs PBMC (GSE287288)",
       x = "Cell Type", y = "CPI (%)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("GSE167232_BAL" = "#E41A1C", "GSE287288_DTB_PBMC" = "#377EB8"))

ggsave(file.path(res_dir, "figures", "Fig8_CPI_BAL_vs_DTB_PBMC.png"), p_cpi, width = 10, height = 6, dpi = 300)

message("\n=== Analysis Complete ===")
message("Outputs:")
message("  - ", file.path(res_dir, "tables", "DEG_GSE287288_DTB.csv"))
message("  - ", file.path(res_dir, "tables", "Table_CPI_GSE287288_DTB.csv"))
message("  - ", file.path(res_dir, "figures", "Fig7_DTB_umap_celltypes.png"))
message("  - ", file.path(res_dir, "figures", "Fig8_CPI_BAL_vs_DTB_PBMC.png"))
