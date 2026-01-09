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

# 08_metabolic_epigenetic_link.R (NEW)
seu <- readRDS(file.path(int_dir, "seurat_objects", "tb_scrna_trained_immunity_scored.rds"))
sets <- jsonlite::fromJSON(cfg$modules$metabolism$gene_list_file)

for (nm in names(sets)) {
  seu <- AddModuleScore(seu, features = list(sets[[nm]]), name = paste0("ms_", nm))
}

vars <- c("condition","cell_type")
for (nm in names(sets)) vars <- c(vars, paste0("ms_", nm, "1"))
df <- FetchData(seu, vars = vars)

write_table(df, file.path(res_dir, "tables", "Table_metabolic_scores.csv"))

if ("ms_glycolysis1" %in% colnames(df)) {
  p <- ggplot(df, aes(x=condition, y=ms_glycolysis1)) + geom_boxplot() + theme_bw() +
    facet_wrap(~cell_type, scales="free_y") + ggtitle("Glycolysis module score by condition")
  ggsave(file.path(res_dir, "figures", "Fig5B_metabolic_scores.png"), p, width=10, height=6, dpi=300)
}

cpi_path <- file.path(res_dir, "tables", "Table_CPI_by_celltype.csv")
if (file.exists(cpi_path)) {
  cpi <- fread(cpi_path)
  summ <- df %>% group_by(cell_type) %>%
    summarise(across(starts_with("ms_"), ~median(.x, na.rm=TRUE))) %>% as.data.frame()
  out <- merge(cpi, summ, by="cell_type", all.x=TRUE)
  write_table(out, file.path(res_dir, "tables", "Table_metabolism_vs_CPI.csv"))
}

saveRDS(seu, file.path(int_dir, "seurat_objects", "tb_scrna_metabolic_scored.rds"))
