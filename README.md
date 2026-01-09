# TB Chromatin Priming Multiomics

[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.XXXXXX-blue)](https://doi.org/10.5281/zenodo.XXXXXX)
[![R](https://img.shields.io/badge/R-%3E%3D4.2.0-blue)](https://www.r-project.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Chromatin Priming Index (CPI): A Quantitative Metric for Integrating Single-Cell Transcriptomic and Chromatin Accessibility Data

### Author
**Dr. Siddalingaiah H S, MD**  
Professor, Department of Community Medicine  
Shridevi Institute of Medical Sciences and Research Hospital  
Tumkur, Karnataka, India  
ORCID: [0000-0002-4771-8285](https://orcid.org/0000-0002-4771-8285)

---

## Overview

This repository contains the analysis pipeline for the **Chromatin Priming Index (CPI)** methodology paper. CPI provides a quantitative framework for assessing chromatin priming from single-cell RNA-seq data using reference chromatin accessibility atlases, without requiring matched ATAC-seq profiling.

### Key Findings

| Dataset | Cells | Mean CPI |
|---------|-------|----------|
| BAL (GSE167232) | 10,357 | 78.8% |
| DTB PBMC (GSE287288) | 21,000 | 84.1% |

- **16 shared master regulator TFs** across datasets: CEBPB, IRF7, STAT1, JUN, etc.
- **Dendritic cells** show highest CPI (89.0%)

---

## Repository Structure

```
TB-Chromatin-Priming-Multiomics/
├── 0_docs/                      # Documentation
│   ├── protocol.md             # Study protocol
│   ├── manuscript_scaffold.md  # Manuscript structure
│   └── figure_plan.md          # Figure specifications
├── 1_data_raw/                  # Raw data (not tracked)
│   ├── GSE167232/              # BAL scRNA-seq
│   └── GSE287288/              # DTB PBMC scRNA-seq
├── 2_data_intermediate/         # Processed data
│   ├── atac_reference/         # Peak-gene links
│   └── seurat_objects/         # Processed Seurat objects
├── 3_analysis/                  # R analysis scripts
│   ├── 01_build_objects.R      # Data loading
│   ├── 02_qc_and_integration.R # QC and normalization
│   ├── 03_wnn_clustering_annotation.R
│   ├── 04_proportion_tests.R
│   ├── 05_deg_dar.R            # Differential expression
│   ├── 06_chromatin_priming_index.R  # CPI calculation
│   ├── 07_trained_immunity_module.R
│   ├── 08_metabolic_epigenetic_link.R
│   ├── 09_tf_regulatory_networks.R
│   ├── 11_metabolic_epigenetic_correlation.R
│   ├── 12_gse287288_dtb_analysis.R  # Second dataset
│   ├── 13_tf_master_regulator_analysis.R
│   └── run_all.ps1             # Pipeline orchestrator
├── 4_results/                   # Output files
│   ├── figures/                # Publication figures
│   ├── tables/                 # Data tables
│   └── submission/             # Manuscript files
├── config/                      # Configuration
│   ├── config.yaml             # Main configuration
│   ├── trained_immunity_genes.txt
│   └── metabolic_genesets.json
├── .github/workflows/           # CI/CD
│   └── r-analysis.yml
├── renv/                        # R package management
├── renv.lock                    # Package versions
├── requirements.txt             # Python dependencies
├── LICENSE
└── README.md
```

---

## Quick Start

### Prerequisites

- R >= 4.2.0
- Python >= 3.8
- Git

### Installation

```bash
# Clone repository
git clone https://github.com/hssling/TB_Chromatin_Priming_Multiomics.git
cd TB_Chromatin_Priming_Multiomics

# Install R dependencies
Rscript -e "install.packages('renv'); renv::restore()"

# Install Python dependencies
pip install -r requirements.txt
```

### Download Data

```bash
# Download GSE167232 (BAL)
# Data is automatically downloaded when running the pipeline

# Download GSE287288 (DTB PBMC)
curl -L -o 1_data_raw/GSE287288_RAW.tar \
  "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE287288&format=file"
tar -xf 1_data_raw/GSE287288_RAW.tar -C 1_data_raw/GSE287288/

# Download ATAC reference (10x PBMC multiome)
curl -L -o 2_data_intermediate/atac_reference/feature_linkage.bedpe.gz \
  "https://cf.10xgenomics.com/samples/cell-arc/2.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_feature_linkage.bedpe.gz"
python convert_bedpe.py
```

### Run Analysis

```bash
# Run full pipeline
cd 3_analysis
./run_all.ps1  # Windows
# OR
for script in *.R; do Rscript "$script"; done  # Linux/Mac
```

---

## Reproducibility

### Environment

This analysis was performed with:
- R 4.5.1
- Seurat 5.x
- Python 3.11

All package versions are locked in `renv.lock`.

### Verification

To verify results:

```bash
# Run tests
Rscript tests/verify_cpi_calculation.R

# Check output checksums
md5sum 4_results/tables/Table_CPI_by_celltype.csv
# Expected: [checksum will be added]
```

---

## Data Sources

| Source | Accession | Description |
|--------|-----------|-------------|
| GEO | [GSE167232](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE167232) | BAL macrophages (Pisu et al. 2021) |
| GEO | [GSE287288](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE287288) | DTB PBMC (Gong et al. 2025) |
| 10x Genomics | PBMC Multiome | Peak-gene links (1.28M) |

---

## Citation

If you use this pipeline or the CPI methodology, please cite:

```bibtex
@article{siddalingaiah2025cpi,
  title={Chromatin Priming Index: A Quantitative Metric for Integrating 
         Single-Cell Transcriptomic and Chromatin Accessibility Data in Tuberculosis},
  author={Siddalingaiah, H S},
  journal={Frontiers in Immunology},
  year={2025},
  doi={10.3389/fimmu.2025.XXXXXX}
}
```

---

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

## Contact

**Dr. Siddalingaiah H S**  
Email: hssling@yahoo.com  
ORCID: [0000-0002-4771-8285](https://orcid.org/0000-0002-4771-8285)
