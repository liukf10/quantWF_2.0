# quantWF_2.0 - Quantitative Omics Analysis Workflow

[![R Version](https://img.shields.io/badge/R-%3E%3D%203.6-blue)](https://www.r-project.org/)
[![License](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Maintenance](https://img.shields.io/badge/Maintained%3F-yes-green.svg)](https://github.com/yourname/quantWF_2.0/graphs/commit-activity)

**quantWF_2.0** is a comprehensive R-based bioinformatics workflow framework designed for quantitative omics data analysis, including RNA-seq, mass spectrometry (MS)-based proteomics, and other quantitative data types.



> 🚀 **Version 2.1** | 📧 Maintainer: Kefu Liu (liukefu19@163.com) | 📅 March 2026

---

## 📋 Table of Contents

- [Features](#-features)
- [Changelog](#-changelog)
- [Installation](#-installation)
- [Quick Start](#-quick-start)
- [Workflow Overview](#-workflow-overview)
- [Parameters](#-parameters)
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

## 📝 Changelog

### v2.1 (March 2026)

**New Features:**
- **Added `outputdir` parameter**: Support for specifying output directory
  - All output files (including plot directory) are now saved to the specified directory
  - Default is current working directory (`.`) when not specified
  - Both R script versions (`outputdir` variable) and command-line versions (`--outputdir` flag) are supported

- **Dual Output Directories for Graphics**: Organized plot outputs into two directories
  - `plot/`: High-resolution PDF files suitable for publication
  - `assets/`: PNG files optimized for Markdown report embedding

- **Automatic Markdown Report Generation**: Each step now generates comprehensive execution reports
  - Report file: `{output}_report.md` in the output directory
  - Includes executed code, parameters, screen output, and loginfo
  - Embeds PNG plots from `assets/` directory
  - Lists all generated files with descriptions
  - Supports appending when chaining steps (Step 1 → Step 2 → Step 3)
  - **PCA plots are now ordered** in reports: raw → cpm → rm → norm → INT

**Improvements:**
- Added `get_plot_dir()` helper function in `preprocess_func.R` for automatic plot directory creation
- Added report generation functions in `preprocess_func.R`:
  - `init_report()`: Initialize report file
  - `write_executed_code()`: Record executed R code
  - `write_parameters()`: Document parameter settings
  - `write_loginfo()`: Format execution loginfo as YAML
  - `pdf_to_png()`: Convert PDF plots to PNG for embedding
  - `insert_plot()`: Insert plot images into report
  - `write_generated_files()`: List all output files
  - `check_existing_report()`: Support for report appending
- Updated all step scripts (Step 0-3) to support the new `outputdir` parameter and report generation
- File paths are now constructed using `file.path()` for better cross-platform compatibility
- Fixed `makeScreePlot()` blank output issue by adding explicit `print()` calls
- Fixed `assets_dir` scope issues in step1 scripts
- Fixed gene names with special characters (e.g., "H3-3A") in BIC evaluation

**Files Modified:**
- `preprocess_func.R`: Added `get_plot_dir()` and report generation functions, fixed `makeScreePlot()`
- All `step*_R` and `step*_cmd` files: Added `outputdir` parameter, report generation, and `assets/` directory support
- `wgcna_plot.R`: Added `outputdir` parameter support
- `DEG_troika.R`: Version bump to 2.1, added `outputdir` parameter

### v2.0 (July 2025)

- Initial release with unified interface
- Five-step workflow: Data Extraction → Preprocessing → Covariate Adjustment → Statistical Analysis → Functional Enrichment
- Dual interface modes: Command-line (`*_cmd`) and R-script (`*_R`) versions
- Support for RNA-seq, MS proteomics, and other quantitative data types

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

## 🔧 Parameters

### Common Parameters (All Steps)

| Parameter | Type | Description | Default |
|-----------|------|-------------|---------|
| `output` | string | Output file name prefix | `"data"` |
| `outputdir` | string | Output directory path | `"."` (current directory) |
| `funcdir` | string | Directory containing `preprocess_func.R` | Current directory |

### Step 0: MS Data Extraction (`step0_MS_R` / `step0_MS_cmd`)

| Parameter | Type | Description | Default |
|-----------|------|-------------|---------|
| `quantfile` | string | Input MS data file (txt/csv/RData) | Required |
| `covarfile` | string | Sample metadata file | `NA` |
| `IDinfCol` | string/vector | Column(s) containing feature IDs | `"Protein.IDs"` |
| `output` | string | Output prefix | `"data"` |
| `outputdir` | string | Output directory | `"."` |

### Step 1: RNA-seq Preprocessing (`step1_rna_R` / `step1_rna_cmd`)

| Parameter | Type | Description | Default |
|-----------|------|-------------|---------|
| `quantfile` | string | Input count matrix file | Required |
| `covarfile` | string | Sample metadata file | `NA` |
| `output` | string | Output prefix | `"data"` |
| `outputdir` | string | Output directory | `"."` |
| `step` | string | Steps to execute (1=CPM, 2=outlier, 3=filter, 4=norm) | `"234"` |
| `sdout` | numeric | Outlier detection SD threshold | `3` |
| `min.value` | numeric | Minimum expression threshold | `1` |
| `min.total` | numeric | Minimum total count threshold | `3*min.value` |
| `grpname` | string | Grouping variable name in covar | `"group"` |
| `norm.method` | string | Normalization method ("TMM" or "quantile") | `"TMM"` |

### Step 1: RNA-seq with TPM (`step1_countTPM_R` / `step1_countTPM_cmd`)

| Parameter | Type | Description | Default |
|-----------|------|-------------|---------|
| `quantfile` | string | Input file with count and TPM matrices | Required |
| `covarfile` | string | Sample metadata file | `NA` |
| `output` | string | Output prefix | `"data"` |
| `outputdir` | string | Output directory | `"."` |
| `step` | string | Steps to execute | `"234"` |
| `sdout` | numeric | Outlier detection SD threshold | `3` |
| `min.value` | numeric | Minimum TPM threshold | `0.1` |
| `min.total` | numeric | Minimum total count threshold | `3*min.value` |
| `grpname` | string | Grouping variable name | `"group"` |
| `norm.method` | string | Normalization method | `"TMM"` |

### Step 1: MS Proteomics (`step1_MS_R` / `step1_MS_cmd`)

| Parameter | Type | Description | Default |
|-----------|------|-------------|---------|
| `quantfile` | string | Input intensity matrix file | Required |
| `covarfile` | string | Sample metadata file | `NA` |
| `output` | string | Output prefix | `"data"` |
| `outputdir` | string | Output directory | `"."` |
| `step` | string | Steps (1=norm, 2=outlier, 3=filter, 4=impute) | `"1234"` |
| `sdout` | numeric | Outlier detection threshold | `3` |
| `maxNA` | numeric | Maximum missing value ratio | `0.5` |
| `grpname` | string | Grouping variable name | `"group"` |

### Step 2: Covariate Adjustment (`step2_covariate_adj_R` / `step2_covariate_adj_cmd`)

| Parameter | Type | Description | Default |
|-----------|------|-------------|---------|
| `quantfile` | string | Input quant file from Step 1 | Required |
| `covarfile` | string | Sample metadata file | Required |
| `output` | string | Output prefix | `"data"` |
| `outputdir` | string | Output directory | `"."` |
| `keepvar` | string | Variables to preserve (comma-separated) | Required |
| `method` | string | Hidden covariate method ("sva", "peer", "pca") | `"sva"` |
| `adjcov` | string | Adjustment mode ("obs&hid", "obs", "hid", "obs+hid") | `"obs&hid"` |
| `svamethod` | string | SVA number estimation ("be", "leek") | `"be"` |
| `nPCAchoose` | string | PC selection method ("BE", "Elbow") | `"BE"` |
| `evaluate_method` | string | Evaluation metrics ("pvca", "pca", "BIC") | `"pvca&pca"` |

### Step 2: Batch Adjustment (`step2_batch_adj_R` / `step2_batch_adj_cmd`)

| Parameter | Type | Description | Default |
|-----------|------|-------------|---------|
| `quantfile` | string | Input quant file | Required |
| `covarfile` | string | Sample metadata file | Required |
| `output` | string | Output prefix | `"data"` |
| `outputdir` | string | Output directory | `"."` |
| `batch` | string | Batch variable name | Required |
| `model` | string | Adjustment model ("combat", "limma") | `"combat"` |
| `evaluate_method` | string | Evaluation metrics | `"pvca&pca"` |

### Step 2: Covariates for QTL (`step2_cov4QTL_R` / `step2_cov4QTL_cmd`)

| Parameter | Type | Description | Default |
|-----------|------|-------------|---------|
| `quantfile` | string | Input quant file | Required |
| `covarfile` | string | Sample metadata file | Required |
| `output` | string | Output prefix | `"data"` |
| `outputdir` | string | Output directory | `"."` |
| `genotype_pc` | string | Genotype PC file | `NA` |
| `ngeno_pc` | numeric | Number of genotype PCs | `3` |
| `method` | string | Hidden covariate method ("sva", "peer", "pca") | `"sva"` |
| `n_hid` | numeric | Number of hidden covariates | `20` |
| `saveannotation` | logical | Save annotation files | `TRUE` |

### Step 3: Differential Expression Classic (`step3_DEG_classic_R` / `step3_DEG_classic_cmd`)

| Parameter | Type | Description | Default |
|-----------|------|-------------|---------|
| `input` | string | Input file from Step 2 | Required |
| `output` | string | Output prefix | `"data"` |
| `outputdir` | string | Output directory | `"."` |
| `groupA` | string | Control group name | Required |
| `groupB` | string | Treatment group name | Required |
| `grpname` | string | Group variable name in covar | `"group"` |
| `log2FCcut` | numeric | Log2 fold change cutoff | `1` |
| `adjPcut` | numeric | Adjusted p-value cutoff | `0.05` |
| `do_ttest` | logical | Also run t-test for large samples | `FALSE` |
| `varequal` | logical | Assume equal variance | `TRUE` |

### Step 3: WGCNA (`step3_WGCNA`)

| Parameter | Type | Description | Default |
|-----------|------|-------------|---------|
| `quantfile` | string | Input expression matrix | Required |
| `covarfile` | string | Sample trait file | `NA` |
| `output` | string | Output prefix | `"data"` |
| `outputdir` | string | Output directory | `"."` |
| `corrType` | string | Correlation type ("bicor", "pearson") | `"bicor"` |
| ` RsquaredCut` | numeric | R-squared cutoff for scale-free topology | `0.8` |
| ` MEDissThres` | numeric | Module dissimilarity threshold | `0.25` |

### DEG Troika (`DEG_troika.R`)

| Parameter | Type | Description | Default |
|-----------|------|-------------|---------|
| `counts` | matrix | Raw count matrix | Required |
| `covar` | data.frame | Sample metadata | Required |
| `output` | string | Output prefix | `"data"` |
| `outputdir` | string | Output directory | `"."` |
| `grpname` | string | Group variable name | `"group"` |
| `groupA` | string | Control group name | Required if >2 groups |
| `groupB` | string | Treatment group name | Required if >2 groups |
| `log2FCcut` | numeric | Log2 fold change cutoff | `1` |
| `adjPcut` | numeric | Adjusted p-value cutoff | `0.05` |

---

## 📚 Documentation

- **Detailed Usage Guide**: See [USAGE.md](USAGE.md) for comprehensive Chinese documentation
- **Test Scripts**: Check `testdata/测试脚本R中运行.R` for working examples
- **Function Reference**: See `preprocess_func.R` for utility functions

### 📊 Report Generation

Each step automatically generates a Markdown report documenting the entire analysis process.

#### Report Location
```
{outputdir}/{output}_report.md
```

#### Report Contents

The report includes:

1. **Executed Code**: The exact R code or command used
2. **Parameters**: Table of all parameter settings
3. **Step-by-Step Log**:
   - Screen output (cat/print/message)
   - Execution loginfo (YAML format)
   - Embedded plots (PNG format from `assets/` directory)
4. **Generated Files Summary**: List of all output files with descriptions

**Output Directory Structure**:
```
{outputdir}/
├── {output}_report.md          # Markdown report
├── *.RData                     # Analysis results
├── plot/                       # High-resolution PDF plots
│   ├── *_step1_pca.pdf
│   ├── *_step1_distribution.pdf
│   ├── *_step2_*_pca.pdf
│   ├── *_step2_*_pvca.pdf
│   └── *_step2_*_bic.pdf
└── assets/                     # PNG plots for report embedding
    ├── *_step1_pca_*.png      # Step 1 PCA plots (ordered: raw→cpm→rm→norm→INT)
    ├── *_step1_box_*.png      # Step 1 boxplots
    ├── *_step1_density_*.png  # Step 1 density plots
    ├── *_step2_*_pca_*.png    # Step 2 PCA plots by covariate
    ├── *_step2_*_pvca.png
    ├── *_step2_*_bic.png
    └── *_step2_box_*.png      # Step 2 distribution boxplots
```

#### Example Report Structure

```markdown
# quantWF Analysis Report

## Step 1: RNA-seq Preprocessing
Generated: 2026-03-25 10:00:00

### Executed Code
```r
quantfile <- "testdata/countTPM.RData"
covarfile <- "testdata/countTPM_covar.txt"
output <- "my_analysis"
source("step1_countTPM_R")
```

### Parameters
| Parameter | Value |
|-----------|-------|
| step | 234 |
| norm.method | TMM |
...

### Step 1.1: CPM Normalization

**Output:**
```
Finish cpm normalized.
```

**loginfo:**
```yaml
record: "CPM normalization completed"
input: testdata/countTPM.RData
output: my_analysis_cpm.RData
running_time: 2026-03-25 10:00:01
```

### Step 1.2: Outlier Detection
...

**Plot:**
![Outlier Detection](assets/my_analysis_step1_pca_raw_group.png)

### Generated Files Summary

| File | Path | Description |
|------|------|-------------|
| **Final Result** | **my_analysis_log_TMM.RData** | **For Step 2 input** |
| CPM Result | my_analysis_cpm.RData | After CPM normalization |
...
```

#### Chaining Steps

When running steps sequentially, reports are automatically appended:

```r
# Step 1 creates the report
source("step1_countTPM_R")

# Step 2 appends to the same report
source("step2_covariate_adj_R")

# Step 3 continues appending
source("step3_DEG_classic_R")

# Final report contains all three steps
```

#### Requirements for Plot Embedding

To embed plots in the report, install one of:
- `pdftools` package: `install.packages("pdftools")`
- `magick` package: `install.packages("magick")`

If neither is installed, plots will still be generated as PDF but won't be embedded in the report.

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
- **Kimi code** help to finish readme and tutorial md document and coding in >V2.1.

---

<div align="center">

**Made with ❤️ by Kefu Liu**

</div>
