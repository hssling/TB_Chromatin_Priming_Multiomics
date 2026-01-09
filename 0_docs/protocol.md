# Protocol (v2): TB Chromatin Priming Multiomics

## Study objective
Quantify cell-type–specific transcriptional changes and chromatin-priming signatures across TB disease states and identify candidate master regulators.

## Primary contrasts
- Active TB vs LTBI (preferred)
- Active TB vs Healthy (if LTBI not available)
- Optional: severe vs non-severe TB (if metadata supports)

## Key derived endpoints
- Differential expression (DEG) per cell type
- Accessibility linkage proxy (DEG ↔ nearby accessible peaks from reference)
- Chromatin Priming Index (CPI) per cell type
- Trained immunity & metabolic module scores
- TF prioritization & TF→target networks with epigenetic priors

## Reproducibility
All thresholds, paths, and dataset identifiers are defined in `config/config.yaml`.
