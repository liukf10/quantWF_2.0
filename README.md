# quantWF_2.0 - Quantitative Omics Analysis Workflow

[![R Version](https://img.shields.io/badge/R-%3E%3D%203.6-blue)](https://www.r-project.org/)
[![License](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Maintenance](https://img.shields.io/badge/Maintained%3F-yes-green.svg)](https://github.com/yourname/quantWF_2.0/graphs/commit-activity)

**quantWF_2.0** is a comprehensive R-based bioinformatics workflow framework designed for quantitative omics data analysis, including RNA-seq, mass spectrometry (MS)-based proteomics, and other quantitative data types.

Kimi code help to finish readme and tutorial md document.

> 🚀 **Version 2.0** | 📧 Maintainer: Kefu Liu (liukefu19@163.com) | 📅 July 2025

---

## 📋 Table of Contents

- [Features](#-features)
- [Installation](#-installation)
- [Quick Start](#-quick-start)
- [Workflow Overview](#-workflow-overview)
- [Documentation](#-documentation)
- [File Structure](#-file-structure)
- [Citation](#-citation)

---

## ✨ Features

quantWF_2.0 provides an **end-to-end analysis pipeline** with five major steps:

| Step | Module | Description | Key Methods |
|:---:|:---:|:---|:---|
| **0** | Data Extraction | Extract quantitative data from MS output files | MaxQuant, Proteome Discoverer parser |
| **1** | Preprocessing | Normalization, outlier detection, filtering, imputation | CPM, TMM, quantile, KNN, IAC/bicor outlier detection |
| **2** | Covariate Adjustment | Batch correction and hidden covariate estimation | SVA, PEER, PCA, Combat, limma::removeBatchEffect |
| **3** | Statistical Analysis | Differential expression and co-expression network | limma/edgeR/DESeq2, WGCNA, t-test/wilcoxon |
| **4** | Functional Enrichment | GO/KEGG and GSEA analysis | clusterProfiler |

### 🎯 Key Highlights

- **Dual Interface Modes**: Each step provides both command-line (`*_cmd`) and R-script (`*_R`) versions for flexibility
- **Automatic Covariate Matching**: Smart alignment between samples and metadata
- **Hidden Covariate Discovery**: Three methods (SVA, PEER, PCA) to capture unknown confounders
- **Comprehensive QC**: Built-in PCA, PVCA, BIC, MDS, and density plots for quality assessment
- **QTL-Ready**: Special modules for QTL analysis with INT transformation and tensorQTL/QTLtools output

---

## 🔧 Installation

### System Requirements

- R >= 3.6.0
- Bioconductor >= 3.10
- Linux/macOS/Windows (PEER method requires Linux)

### Install Dependencies

```r
# Run the package installation script
source("install_dependencies.R")
```

Or install manually:

```r
# CRAN packages
install.packages(c(
  "optparse", "DDPNA", "data.table", "remotes", 
  "ggplot2", "ggpubr", "ggrepel", "doParallel",
  "lme4", "nlme", "foreach", "devtools", 
  "WGCNA", "corrplot"
))

# Bioconductor packages
BiocManager::install(c(
  "Biostrings", "sva", "rtracklayer", "DESeq2",
  "limma", "edgeR", "vsn", "impute", "Biobase",
  "org.Hs.eg.db", "org.Mm.eg.db", "org.Rn.eg.db",
  "clusterProfiler"
))
```

### Clone Repository

```bash
git clone https://github.com/yourname/quantWF_2.0.git
cd quantWF_2.0
```

---

## 🚀 Quick Start

### Example 1: RNA-seq Differential Expression Analysis

```r
setwd("quantWF_2.0")

# Step 1: Preprocess RNA-seq count data
quantfile <- "testdata/RNA_data.RData"
covarfile <- "testdata/countTPM_covar.txt"
output <- "my_analysis"
grpname <- "group"
norm.method <- "TMM"
step <- "234"  # 2=outlier removal, 3=gene filter, 4=normalization
source("step1_rna_R")

# Step 2: Adjust for hidden covariates
quantfile <- "my_analysis_filter_TMM.RData"
covarfile <- "testdata/countTPM_covar.txt"
keepvar <- "group"
method <- "sva"
adjcov <- "obs&hid"
evaluate_method <- "pvca&pca"
source("step2_covariate_adj_R")

# Step 3: Differential expression analysis
input <- "my_analysis_adj_final.RData"
groupA <- "control"
groupB <- "treatment"
log2FCcut <- 1
adjPcut <- 0.05
source("step3_DEG_classic_R")
```

### Example 2: Command Line Usage

```bash
# Step 1: RNA-seq preprocessing
Rscript step1_rna_cmd \
  --quantfile="counts.RData" \
  --covarfile="metadata.RData" \
  --output="rna_analysis" \
  --grpname="group" \
  --step="234" \
  --norm.method="TMM"

# Step 2: Covariate adjustment
Rscript step2_covariate_adj_cmd \
  --quantfile="rna_analysis_filter_TMM.RData" \
  --covarfile="metadata.RData" \
  --output="adj_analysis" \
  --keepvar="diagnose" \
  --method="sva" \
  --adjcov="obs&hid"

# Step 3: Differential expression
Rscript step3_DEG_classic_cmd \
  --input="adj_analysis_adj_final.RData" \
  --output="deg_results" \
  --groupA="control" \
  --groupB="treatment"
```

### Example 3: MS Proteomics Analysis

```r
# Step 0: Extract from raw output
quantfile <- "testdata/phos_data.txt"
IDinfCol <- 1:7
output <- "phos_analysis"
source("step0_MS_R")

# Step 1: Preprocess MS data
quantfile <- "phos_analysis_raw.RData"
step <- "1234"
sdout <- 3
maxNA <- 0.5
source("step1_MS_R")
```

---

## 📊 Workflow Overview

### Data Flow

```
Raw Data → Step 0 (Extract) → Step 1 (Preprocess) → Step 2 (Adjust) → Step 3 (Analyze) → Step 4 (Enrich)
     ↓                              ↓                       ↓                  ↓                ↓
  .txt/.csv                  _raw.RData          _adj_obs.RData      _DEG.RData     enrichment
  (MaxQuant,                → _rmOut.RData       _adj_final.RData    _stat.RData    results
  counts, etc.)             → _filter.RData
                            → _norm.RData
```

### Supported Data Types

| Data Type | Step 1 Script | Input Format | Key Features |
|:---|:---|:---|:---|
| RNA-seq counts | `step1_rna_R` | Count matrix | CPM → filter → TMM/quantile |
| RNA-seq counts + TPM | `step1_countTPM_R` | count + TPM | GTEx-style filtering |
| Expression matrix | `step1_rna_R` | Normalized matrix | Log2 input supported |
| MS Proteomics | `step1_MS_R` | Intensity matrix | Total intensity norm → KNN imputation |

---

## 📚 Documentation

- **Detailed Usage Guide**: See [USAGE.md](USAGE.md) for comprehensive Chinese documentation
- **Test Scripts**: Check `testdata/测试脚本R中运行.R` for working examples
- **Function Reference**: See `preprocess_func.R` for utility functions

### Key Parameters Reference

#### Step 1 Parameters
| Parameter | Description | Default |
|-----------|-------------|---------|
| `step` | Steps to run (e.g., "234") | "234" |
| `sdout` | Outlier SD threshold | 3 |
| `min.value` | Minimum expression threshold | 1 (RNA), 0.1 (countTPM) |
| `norm.method` | Normalization method | "TMM" |

#### Step 2 Parameters
| Parameter | Description | Options |
|-----------|-------------|---------|
| `method` | Hidden covariate method | "sva", "peer", "pca" |
| `adjcov` | Adjustment mode | "obs&hid", "obs", "hid", "obs+hid" |
| `keepvar` | Variables to preserve | comma-separated |

---

## 📁 File Structure

```
quantWF_2.0/
├── preprocess_func.R              # Core utility functions
├── DEG_troika.R                   # limma/edgeR/DESeq2 DEG analysis
├── wgcna_plot.R / wgcna_tuning.R  # WGCNA visualization & tuning
├── step0_MS_cmd / step0_MS_R      # Step 0: MS data extraction
├── step1_countTPM_cmd / _R        # Step 1: RNA-seq with count+TPM
├── step1_MS_cmd / step1_MS_R      # Step 1: MS-based proteomics
├── step1_rna_cmd / step1_rna_R    # Step 1: RNA-seq preprocessing
├── step2_batch_adj_cmd / _R       # Step 2: Batch adjustment
├── step2_cov4QTL_cmd / _R         # Step 2: Covariates for QTL
├── step2_covariate_adj_cmd / _R   # Step 2: General covariate adj
├── step3_DEG_classic_cmd / _R     # Step 3: Classic DEG (t/wilcox)
├── step3_WGCNA                    # Step 3: WGCNA analysis
├── step4_GO&KEGGenrich.R          # Step 4: GO/KEGG enrichment
├── step4_GS_GSenrich.R            # Step 4: Gene set enrichment
├── script/                        # Alternative script versions
├── testdata/                      # Test datasets
└── temp/                          # Temporary files
```

---

## 🤝 Contributing

Contributions are welcome! Please feel free to submit a Pull Request. For major changes, please open an issue first to discuss what you would like to change.

### Development Guidelines

1. Follow the existing code style with Chinese comments
2. Maintain backward compatibility
3. Update documentation for any new features
4. Add test cases to `testdata/` for new functionality

---

## 📖 Citation

If you use quantWF_2.0 in your research, please cite:

```
Liu K. (2025). quantWF_2.0: A comprehensive R-based workflow framework for quantitative omics data analysis. GitHub repository.
```

---

## 📝 License

This project is licensed under the MIT License - see the LICENSE file for details.

---

## 🙏 Acknowledgments

- **SVA** package for hidden covariate estimation (Leek et al.)
- **PEER** method from GTEx project
- **WGCNA** for co-expression network analysis (Langfelder & Horvath)
- **clusterProfiler** for functional enrichment (Yu et al.)
- **edgeR/limma/DESeq2** for differential expression analysis

---

<div align="center">

**Made with ❤️ by Kefu Liu**

</div>
