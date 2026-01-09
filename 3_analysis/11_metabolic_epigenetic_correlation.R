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

int_dir <- cfg$paths$intermediate
res_dir <- cfg$paths$results

# 11_metabolic_epigenetic_correlation.R
# Novel analysis: Are metabolic/trained immunity genes enriched among chromatin-primed DEGs?

message("=== Metabolic-Epigenetic Correlation Analysis ===")

# Load data
peak_links <- fread("2_data_intermediate/atac_reference/peak_gene_links.csv")
deg_all <- fread(file.path(res_dir, "tables", "DEG_all_celltypes.csv"))
metabolic_genes <- fromJSON("config/metabolic_genesets.json")
trained_genes <- readLines("config/trained_immunity_genes.txt")
trained_genes <- trained_genes[nchar(trained_genes) > 0]

message("Peak-gene links: ", nrow(peak_links))
message("DEGs total: ", nrow(deg_all))

# Combine all metabolic pathway genes
all_metabolic <- unique(unlist(metabolic_genes))
message("Metabolic pathway genes: ", length(all_metabolic))
message("Trained immunity genes: ", length(trained_genes))

# Get unique genes with peak links
linked_genes <- unique(peak_links$gene)
message("Genes with chromatin links: ", length(linked_genes))

# Calculate enrichment for each gene category
results <- list()

for (ct in unique(deg_all$cell_type)) {
  deg_ct <- deg_all[deg_all$cell_type == ct, ]
  
  # Significant DEGs
  sig_deg <- deg_ct[deg_ct$p_val_adj < 0.05, ]$gene
  
  # DEGs with chromatin links
  primed_deg <- intersect(sig_deg, linked_genes)
  non_primed_deg <- setdiff(sig_deg, linked_genes)
  
  # Metabolic gene enrichment in primed vs non-primed
  metabolic_in_primed <- sum(primed_deg %in% all_metabolic)
  metabolic_in_nonprimed <- sum(non_primed_deg %in% all_metabolic)
  
  trained_in_primed <- sum(primed_deg %in% trained_genes)
  trained_in_nonprimed <- sum(non_primed_deg %in% trained_genes)
  
  # Calculate proportions
  prop_metabolic_primed <- metabolic_in_primed / length(primed_deg)
  prop_metabolic_nonprimed <- metabolic_in_nonprimed / max(1, length(non_primed_deg))
  
  prop_trained_primed <- trained_in_primed / length(primed_deg)
  prop_trained_nonprimed <- trained_in_nonprimed / max(1, length(non_primed_deg))
  
  # Fisher's exact test for metabolic genes
  fisher_metabolic <- fisher.test(matrix(c(
    metabolic_in_primed, 
    length(primed_deg) - metabolic_in_primed,
    metabolic_in_nonprimed, 
    max(1, length(non_primed_deg)) - metabolic_in_nonprimed
  ), nrow=2))
  
  # Fisher's exact test for trained immunity genes
  fisher_trained <- fisher.test(matrix(c(
    trained_in_primed, 
    length(primed_deg) - trained_in_primed,
    trained_in_nonprimed, 
    max(1, length(non_primed_deg)) - trained_in_nonprimed
  ), nrow=2))
  
  results[[ct]] <- data.frame(
    cell_type = ct,
    n_sig_deg = length(sig_deg),
    n_primed = length(primed_deg),
    n_nonprimed = length(non_primed_deg),
    metabolic_primed = metabolic_in_primed,
    metabolic_nonprimed = metabolic_in_nonprimed,
    prop_metabolic_primed = prop_metabolic_primed,
    prop_metabolic_nonprimed = prop_metabolic_nonprimed,
    metabolic_OR = fisher_metabolic$estimate,
    metabolic_pval = fisher_metabolic$p.value,
    trained_primed = trained_in_primed,
    trained_nonprimed = trained_in_nonprimed,
    prop_trained_primed = prop_trained_primed,
    prop_trained_nonprimed = prop_trained_nonprimed,
    trained_OR = fisher_trained$estimate,
    trained_pval = fisher_trained$p.value
  )
  
  message(sprintf("  %s: Metabolic OR=%.2f (p=%.3f), Trained OR=%.2f (p=%.3f)", 
                  ct, fisher_metabolic$estimate, fisher_metabolic$p.value,
                  fisher_trained$estimate, fisher_trained$p.value))
}

# Combine results
results_df <- rbindlist(results)
write_table(results_df, file.path(res_dir, "tables", "Table_metabolic_enrichment_analysis.csv"))
message("Saved enrichment results")

# Per-gene analysis: accessibility score for metabolic vs non-metabolic DEGs
message("\n=== Per-gene accessibility analysis ===")

# Get mean accessibility scores per gene
gene_accessibility <- peak_links[, .(
  mean_score = mean(score, na.rm = TRUE),
  max_score = max(score, na.rm = TRUE),
  n_links = .N
), by = gene]

# Merge with DEG status
sig_degs <- unique(deg_all[deg_all$p_val_adj < 0.05, ]$gene)
gene_accessibility$is_deg <- gene_accessibility$gene %in% sig_degs
gene_accessibility$is_metabolic <- gene_accessibility$gene %in% all_metabolic
gene_accessibility$is_trained <- gene_accessibility$gene %in% trained_genes
gene_accessibility$category <- ifelse(
  gene_accessibility$is_trained, "Trained Immunity",
  ifelse(gene_accessibility$is_metabolic, "Metabolic",
  ifelse(gene_accessibility$is_deg, "Other DEG", "Background"))
)

write_table(gene_accessibility, file.path(res_dir, "tables", "Table_gene_accessibility_categories.csv"))

# Statistical comparisons
metabolic_scores <- gene_accessibility[is_metabolic == TRUE, mean_score]
trained_scores <- gene_accessibility[is_trained == TRUE, mean_score]
other_deg_scores <- gene_accessibility[is_deg == TRUE & is_metabolic == FALSE & is_trained == FALSE, mean_score]
background_scores <- gene_accessibility[is_deg == FALSE, mean_score]

message(sprintf("Mean accessibility - Trained: %.3f, Metabolic: %.3f, Other DEG: %.3f, Background: %.3f",
                mean(trained_scores, na.rm=T), mean(metabolic_scores, na.rm=T),
                mean(other_deg_scores, na.rm=T), mean(background_scores, na.rm=T)))

# Wilcoxon tests
if (length(metabolic_scores) > 2 && length(background_scores) > 2) {
  wilcox_metabolic <- wilcox.test(metabolic_scores, background_scores)
  message(sprintf("Metabolic vs Background: p = %.2e", wilcox_metabolic$p.value))
}

if (length(trained_scores) > 2 && length(background_scores) > 2) {
  wilcox_trained <- wilcox.test(trained_scores, background_scores)
  message(sprintf("Trained vs Background: p = %.2e", wilcox_trained$p.value))
}

# Generate visualization
p1 <- ggplot(gene_accessibility[category != "Background"], 
             aes(x = reorder(category, -mean_score), y = mean_score, fill = category)) +
  geom_boxplot(outlier.size = 0.5) +
  labs(title = "Chromatin Accessibility by Gene Category",
       x = "", y = "Mean Accessibility Score") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("Trained Immunity" = "#E41A1C", 
                               "Metabolic" = "#377EB8", 
                               "Other DEG" = "#4DAF4A"))

ggsave(file.path(res_dir, "figures", "Fig6_metabolic_accessibility.png"), p1, width = 6, height = 5, dpi = 300)
message("Saved accessibility figure")

# Summary statistics
summary_stats <- data.frame(
  category = c("Trained Immunity", "Metabolic", "Other DEG", "Background"),
  mean_accessibility = c(mean(trained_scores, na.rm=T), mean(metabolic_scores, na.rm=T),
                         mean(other_deg_scores, na.rm=T), mean(background_scores, na.rm=T)),
  n_genes = c(length(trained_scores), length(metabolic_scores), 
              length(other_deg_scores), length(background_scores))
)
write_table(summary_stats, file.path(res_dir, "tables", "Table_accessibility_summary.csv"))

message("\n=== Analysis Complete ===")
