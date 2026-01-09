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

# 09_tf_regulatory_networks.R (NEW)
deg <- fread(file.path(res_dir, "tables", "DEG_all_celltypes.csv"))
links <- fread(file.path(int_dir, "atac_reference", "peak_gene_links.csv"))

motif_path <- file.path(int_dir, "atac_reference", "motif_enrichment_homer.csv")
if (!file.exists(motif_path)) {
  motifs <- data.frame(tf=character(), pval=double(), cell_type=character())
  write_table(motifs, motif_path)
} else {
  motifs <- fread(motif_path)
}

deg_tfs <- unique(deg$gene)

mr <- motifs %>%
  mutate(motif_score = ifelse(is.finite(pval) & pval>0, -log10(pval), 0),
         is_deg = tf %in% deg_tfs,
         composite = motif_score + ifelse(is_deg, 1, 0)) %>%
  arrange(desc(composite))

write_table(mr, file.path(res_dir, "tables", "Table_master_regulators.csv"))

if (nrow(mr) > 0) {
  top <- head(mr, 20)
  p <- ggplot(top, aes(x=reorder(tf, composite), y=composite)) + geom_col() + coord_flip() + theme_bw() +
    ggtitle("Master regulator ranking (composite, placeholder)")
  ggsave(file.path(res_dir, "figures", "Fig6B_master_regulator_ranking.png"), p, width=8, height=6, dpi=300)
}
