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

# 07_trained_immunity_module.R (NEW)
seu <- readRDS(file.path(int_dir, "seurat_objects", "tb_scrna_annotated.rds"))
genes <- readLines(cfg$modules$trained_immunity$gene_list_file)
genes <- genes[genes != ""]

seu <- AddModuleScore(seu, features = list(genes), name = "trained_immunity_score")

df <- FetchData(seu, vars = c("trained_immunity_score1", "condition", "cell_type"))
colnames(df)[1] <- "trained_immunity_score"
write_table(df, file.path(res_dir, "tables", "Table_trained_immunity_scores.csv"))

p <- ggplot(df, aes(x=condition, y=trained_immunity_score)) +
  geom_boxplot() + theme_bw() +
  facet_wrap(~cell_type, scales="free_y") +
  ggtitle("Trained immunity module score by condition")
ggsave(file.path(res_dir, "figures", "Fig5A_trained_immunity_score.png"), p, width=10, height=6, dpi=300)

saveRDS(seu, file.path(int_dir, "seurat_objects", "tb_scrna_trained_immunity_scored.rds"))
