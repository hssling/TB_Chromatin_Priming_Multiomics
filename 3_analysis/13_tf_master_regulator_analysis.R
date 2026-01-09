suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(data.table)
})
source("3_analysis/utils/utils_io.R")
cfg <- read_config("config/config.yaml")
set.seed(cfg$analysis$seed)

res_dir <- cfg$paths$results
int_dir <- cfg$paths$intermediate

# 13_tf_master_regulator_analysis.R
# Identify master regulators among chromatin-primed DEGs

message("=== TF Master Regulator Analysis ===")

# Known immune-related transcription factors
known_tfs <- c(
  # NFkB family
  "NFKB1", "NFKB2", "RELA", "RELB", "REL",
  # AP-1 family
  "FOS", "FOSB", "FOSL1", "FOSL2", "JUN", "JUNB", "JUND",
  # IRF family
  "IRF1", "IRF2", "IRF3", "IRF4", "IRF5", "IRF7", "IRF8", "IRF9",
  # STAT family
  "STAT1", "STAT2", "STAT3", "STAT4", "STAT5A", "STAT5B", "STAT6",
  # C/EBP family
  "CEBPA", "CEBPB", "CEBPD", "CEBPE",
  # Other immune TFs
  "SPI1", "ETS1", "ETS2", "ELF1", "ELF4", "SPIB",
  "BATF", "BATF3", "BCL6", "PRDM1", "TBX21", "GATA3",
  "RORC", "FOXP3", "EOMES", "RUNX1", "RUNX2", "RUNX3",
  "HIF1A", "ARNT", "EPAS1", "ATF3", "ATF4", "XBP1",
  "MYC", "MYCN", "MAX", "E2F1", "E2F2", "E2F3",
  # Myeloid TFs
  "PU.1", "MAFB", "MAF", "NR4A1", "NR4A2", "NR4A3"
)

# Load data
deg_bal <- fread(file.path(res_dir, "tables", "DEG_all_celltypes.csv"))
deg_dtb <- fread(file.path(res_dir, "tables", "DEG_GSE287288_DTB.csv"))
links <- fread(file.path(int_dir, "atac_reference", "peak_gene_links.csv"))
linked_genes <- unique(links$gene)

message("BAL DEGs: ", nrow(deg_bal))
message("DTB PBMC DEGs: ", nrow(deg_dtb))
message("Known immune TFs: ", length(known_tfs))

# Analyze TFs in each dataset
results <- list()

for (dataset in c("BAL", "DTB_PBMC")) {
  if (dataset == "BAL") {
    deg_data <- deg_bal
  } else {
    deg_data <- deg_dtb
  }
  
  # Get significant DEGs
  sig_deg <- unique(deg_data[deg_data$p_val_adj < 0.05, ]$gene)
  primed_deg <- intersect(sig_deg, linked_genes)
  
  # TFs that are DEGs
  tf_deg <- intersect(known_tfs, sig_deg)
  
  # TFs that are chromatin-primed DEGs
  tf_primed <- intersect(known_tfs, primed_deg)
  
  # Get expression info for TF DEGs
  for (tf in tf_deg) {
    tf_info <- deg_data[deg_data$gene == tf, ]
    if (nrow(tf_info) > 0) {
      for (i in 1:nrow(tf_info)) {
        results[[length(results) + 1]] <- data.frame(
          dataset = dataset,
          tf = tf,
          cell_type = tf_info$cell_type[i],
          avg_log2FC = tf_info$avg_log2FC[i],
          p_val_adj = tf_info$p_val_adj[i],
          is_primed = tf %in% primed_deg
        )
      }
    }
  }
}

# Combine results
tf_results <- rbindlist(results)
message("\nTF DEGs found: ", length(unique(tf_results$tf)))

# Rank TFs by number of cell types where they're differentially expressed
tf_summary <- tf_results %>%
  group_by(dataset, tf) %>%
  summarize(
    n_cell_types = n_distinct(cell_type),
    mean_log2FC = mean(avg_log2FC, na.rm = TRUE),
    min_pval = min(p_val_adj, na.rm = TRUE),
    is_primed = first(is_primed),
    .groups = "drop"
  ) %>%
  mutate(score = n_cell_types * abs(mean_log2FC) * ifelse(is_primed, 1.5, 1)) %>%
  arrange(desc(score))

write_table(tf_summary, file.path(res_dir, "tables", "Table_TF_master_regulators.csv"))
message("Saved TF summary")

# Top TFs
message("\nTop Master Regulators:")
print(head(tf_summary, 15))

# Visualization - Top TFs by dataset
p1 <- ggplot(head(tf_summary[tf_summary$dataset == "BAL", ], 15), 
             aes(x = reorder(tf, score), y = score, fill = is_primed)) +
  geom_col() +
  coord_flip() +
  labs(title = "BAL: Top TF Master Regulators", x = "", y = "Composite Score") +
  scale_fill_manual(values = c("TRUE" = "#E41A1C", "FALSE" = "#377EB8"),
                    labels = c("TRUE" = "Chromatin-primed", "FALSE" = "Not primed")) +
  theme_bw()

p2 <- ggplot(head(tf_summary[tf_summary$dataset == "DTB_PBMC", ], 15), 
             aes(x = reorder(tf, score), y = score, fill = is_primed)) +
  geom_col() +
  coord_flip() +
  labs(title = "DTB PBMC: Top TF Master Regulators", x = "", y = "Composite Score") +
  scale_fill_manual(values = c("TRUE" = "#E41A1C", "FALSE" = "#377EB8"),
                    labels = c("TRUE" = "Chromatin-primed", "FALSE" = "Not primed")) +
  theme_bw()

p_combined <- p1 / p2
ggsave(file.path(res_dir, "figures", "Fig9_TF_master_regulators.png"), p_combined, width = 10, height = 10, dpi = 300)
message("Saved TF figure")

# TF overlap between datasets
bal_tfs <- unique(tf_summary[tf_summary$dataset == "BAL", ]$tf)
dtb_tfs <- unique(tf_summary[tf_summary$dataset == "DTB_PBMC", ]$tf)
shared_tfs <- intersect(bal_tfs, dtb_tfs)

message("\nTF Summary:")
message("  BAL unique TFs: ", length(setdiff(bal_tfs, dtb_tfs)))
message("  DTB PBMC unique TFs: ", length(setdiff(dtb_tfs, bal_tfs)))
message("  Shared TFs: ", length(shared_tfs))
message("  Shared: ", paste(shared_tfs, collapse = ", "))

message("\n=== Analysis Complete ===")
