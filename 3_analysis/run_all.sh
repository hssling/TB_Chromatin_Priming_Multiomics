#!/usr/bin/env bash
set -euo pipefail

Rscript 3_analysis/01_build_objects.R
Rscript 3_analysis/02_qc_and_integration.R
Rscript 3_analysis/03_wnn_clustering_annotation.R
Rscript 3_analysis/04_proportion_tests.R
Rscript 3_analysis/05_deg_dar.R

# v2 extensions
Rscript 3_analysis/06_chromatin_priming_index.R
Rscript 3_analysis/07_trained_immunity_module.R
Rscript 3_analysis/08_metabolic_epigenetic_link.R
Rscript 3_analysis/09_tf_regulatory_networks.R

Rscript 3_analysis/10_enrichment_and_disease.R

echo "All done. Results in 4_results/"
