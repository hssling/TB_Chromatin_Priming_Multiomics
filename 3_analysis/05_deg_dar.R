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

# 05_deg_dar.R
# Perform differential expression analysis by cell type and condition

seu <- readRDS(file.path(int_dir, "seurat_objects", "tb_scrna_annotated.rds"))
message("Loaded object with ", ncol(seu), " cells")

# Identify condition comparisons
conditions <- unique(seu$condition)
message("Conditions found: ", paste(conditions, collapse=", "))

deg_list <- list()
cell_types <- unique(seu$cell_type)
message("Cell types: ", paste(cell_types, collapse=", "))

for (ct in cell_types) {
  message("Processing cell type: ", ct)
  seu_ct <- subset(seu, subset = cell_type == ct)
  
  if (ncol(seu_ct) < 20) {
    message("  Skipping - too few cells (", ncol(seu_ct), ")")
    next
  }
  
  # Set identities for comparison
  Idents(seu_ct) <- "condition"
  
  # Try to find TB vs control comparison
  conds <- unique(Idents(seu_ct))
  
  # Find appropriate groups for comparison
  tb_like <- conds[grepl("TB|Mtb|infected|disease", conds, ignore.case=TRUE)]
  ctrl_like <- conds[grepl("control|healthy|uninfected|mock", conds, ignore.case=TRUE)]
  
  if (length(tb_like) > 0 && length(ctrl_like) > 0) {
    message("  Comparing ", tb_like[1], " vs ", ctrl_like[1])
    
    tryCatch({
      deg <- FindMarkers(seu_ct, 
                        ident.1 = tb_like[1], 
                        ident.2 = ctrl_like[1],
                        logfc.threshold = cfg$analysis$deg$logfc_threshold,
                        min.pct = cfg$analysis$deg$min_pct,
                        verbose = FALSE)
      
      deg$gene <- rownames(deg)
      deg$cell_type <- ct
      deg$comparison <- paste0(tb_like[1], "_vs_", ctrl_like[1])
      
      deg_list[[ct]] <- deg
      
      # Save individual cell type results
      write_table(deg, file.path(res_dir, "tables", paste0("DEG_", gsub("/|\\s", "_", ct), ".csv")))
      message("  Found ", sum(deg$p_val_adj < 0.05, na.rm=TRUE), " significant DEGs")
      
    }, error = function(e) {
      message("  Error in DEG analysis: ", e$message)
    })
    
  } else if (length(conds) >= 2) {
    # Compare first two conditions
    message("  Comparing ", conds[1], " vs ", conds[2])
    
    tryCatch({
      deg <- FindMarkers(seu_ct,
                        ident.1 = conds[1],
                        ident.2 = conds[2],
                        logfc.threshold = cfg$analysis$deg$logfc_threshold,
                        min.pct = cfg$analysis$deg$min_pct,
                        verbose = FALSE)
      
      deg$gene <- rownames(deg)
      deg$cell_type <- ct
      deg$comparison <- paste0(conds[1], "_vs_", conds[2])
      
      deg_list[[ct]] <- deg
      write_table(deg, file.path(res_dir, "tables", paste0("DEG_", gsub("/|\\s", "_", ct), ".csv")))
      message("  Found ", sum(deg$p_val_adj < 0.05, na.rm=TRUE), " significant DEGs")
      
    }, error = function(e) {
      message("  Error in DEG analysis: ", e$message)
    })
  } else {
    message("  Skipping - insufficient conditions for comparison")
  }
}

# Combine all DEGs
if (length(deg_list) > 0) {
  deg_all <- rbindlist(deg_list, fill=TRUE)
  
  # Filter significant DEGs
  deg_sig <- deg_all[deg_all$p_val_adj < cfg$analysis$deg$padj, ]
  
  write_table(deg_all, file.path(res_dir, "tables", "DEG_all_celltypes.csv"))
  message("Saved all DEGs: ", nrow(deg_all), " total, ", nrow(deg_sig), " significant")
  
  # Generate volcano plot for top cell type
  if (nrow(deg_all) > 0) {
    top_ct <- names(which.max(sapply(deg_list, nrow)))
    deg_top <- deg_list[[top_ct]]
    deg_top$significant <- deg_top$p_val_adj < 0.05 & abs(deg_top$avg_log2FC) > 0.5
    
    p <- ggplot(deg_top, aes(x=avg_log2FC, y=-log10(p_val_adj), color=significant)) +
      geom_point(alpha=0.5, size=1) +
      scale_color_manual(values=c("grey", "red")) +
      theme_bw() +
      ggtitle(paste0("DEGs in ", top_ct)) +
      xlab("Log2 Fold Change") + ylab("-Log10 Adjusted P-value")
    ggsave(file.path(res_dir, "figures", "Fig3_volcano_DEG.png"), p, width=8, height=6, dpi=300)
  }
} else {
  message("No DEGs computed - creating empty file")
  deg_all <- data.frame(gene=character(), p_val_adj=double(), avg_log2FC=double(), cell_type=character())
  write_table(deg_all, file.path(res_dir, "tables", "DEG_all_celltypes.csv"))
}
