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

# 06_chromatin_priming_index.R (NEW)
deg_all <- fread(file.path(res_dir, "tables", "DEG_all_celltypes.csv"))

link_path <- file.path(int_dir, "atac_reference", "peak_gene_links.csv")
if (!file.exists(link_path)) {
  links <- data.frame(gene=character(), peak=character(), cell_type=character())
  write_table(links, link_path)
} else {
  links <- fread(link_path)
}

deg_all$has_peak_link <- deg_all$gene %in% links$gene

cpi <- deg_all %>%
  group_by(cell_type) %>%
  summarise(
    n_deg = n(),
    n_deg_with_link = sum(has_peak_link, na.rm=TRUE),
    CPI = ifelse(n_deg>0, n_deg_with_link / n_deg, NA_real_)
  ) %>% as.data.frame()

write_table(cpi, file.path(res_dir, "tables", "Table_CPI_by_celltype.csv"))

p <- ggplot(cpi, aes(x=reorder(cell_type, CPI), y=CPI)) + geom_col() + coord_flip() + theme_bw() +
  ggtitle("Chromatin Priming Index (CPI) by cell type")
ggsave(file.path(res_dir, "figures", "Fig4A_CPI.png"), p, width=8, height=5, dpi=300)

cats <- deg_all %>%
  mutate(category = ifelse(has_peak_link, "DAR_to_DEG", "DEG_only")) %>%
  count(cell_type, category) %>%
  group_by(cell_type) %>% mutate(prop = n / sum(n)) %>% ungroup()

write_table(cats, file.path(res_dir, "tables", "Table_priming_categories.csv"))

p2 <- ggplot(cats, aes(x=cell_type, y=prop, fill=category)) +
  geom_col() + coord_flip() + theme_bw() +
  ggtitle("Priming categories (placeholder)")
ggsave(file.path(res_dir, "figures", "Fig4B_priming_categories.png"), p2, width=8, height=5, dpi=300)
