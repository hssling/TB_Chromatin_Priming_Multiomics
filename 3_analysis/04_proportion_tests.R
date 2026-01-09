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

# 04_proportion_tests.R
seu <- readRDS(file.path(int_dir, "seurat_objects", "tb_scrna_annotated.rds"))

df <- as.data.frame(table(seu$sample_id, seu$cell_type, seu$condition))
colnames(df) <- c("sample_id","cell_type","condition","n_cells")

df <- df %>% group_by(sample_id) %>% mutate(prop = n_cells / sum(n_cells)) %>% ungroup()
df$padj <- NA_real_

write_table(df, file.path(res_dir, "tables", "Table1_cell_proportions.csv"))

p <- ggplot(df, aes(x=cell_type, y=prop)) + geom_boxplot() + theme_bw() + coord_flip() +
  ggtitle("Cell-type proportions by condition (placeholder)")
ggsave(file.path(res_dir, "figures", "Fig2B_cell_proportions.png"), p, width=8, height=5, dpi=300)
