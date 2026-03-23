# quantWF_2.0 - Quantitative Omics Analysis Workflow

## Project Overview

**quantWF_2.0** is an R-based bioinformatics workflow framework designed for comprehensive quantitative omics data analysis, including RNA-seq, mass spectrometry (MS)-based proteomics, and other quantitative data types. 

- **Version**: 2.0
- **Maintainer**: Kefu Liu (liukefu19@163.com)
- **Language**: R (with Chinese comments and documentation)
- **License**: Not specified

### Key Capabilities

This workflow provides an end-to-end analysis pipeline:

1. **Step 0 - Data Extraction**: Extract quantitative data from raw MS output files
2. **Step 1 - Data Preprocessing**: Normalization, outlier detection, gene/feature filtering, missing value imputation
3. **Step 2 - Covariate Adjustment**: Batch correction and hidden covariate estimation (SVA, PEER, PCA)
4. **Step 3 - Statistical Analysis**: Differential expression analysis (DEG) and Weighted Gene Co-expression Network Analysis (WGCNA)
5. **Step 4 - Functional Enrichment**: GO/KEGG enrichment and GSEA analysis

## Project Structure

```
quantWF_2.0/
├── preprocess_func.R           # Core utility functions (PCA, MDS, density plots, outlier detection, PVCA, BIC)
├── DEG_troika.R               # Differential expression using limma/edgeR/DESeq2
├── wgcna_plot.R               # WGCNA visualization functions
├── wgcna_tuning.R             # WGCNA parameter tuning
├── compareCluster_template.RData  # Template for cluster comparison
│
├── step0_MS_cmd / step0_MS_R          # Step 0: MS data extraction (cmd/R versions)
├── step1_countTPM_cmd / step1_countTPM_R  # Step 1: RNA-seq with count+TPM
├── step1_MS_cmd / step1_MS_R          # Step 1: MS-based proteomics
├── step1_rna_cmd / step1_rna_R        # Step 1: RNA-seq (expression matrix input)
├── step2_batch_adj_cmd / step2_batch_adj_R    # Step 2: Batch adjustment
├── step2_cov4QTL_cmd / step2_cov4QTL_R        # Step 2: Covariates for QTL
├── step2_covariate_adj_cmd / step2_covariate_adj_R  # Step 2: General covariate adjustment
├── step3_DEG_classic_cmd / step3_DEG_classic_R    # Step 3: Classic DEG (t-test/wilcoxon)
├── step3_WGCNA                        # Step 3: WGCNA analysis
├── step4_GO&KEGGenrich.R              # Step 4: GO/KEGG enrichment
├── step4_GS_GSenrich.R                # Step 4: Gene set enrichment
│
├── script/                  # Alternative script versions
├── testdata/               # Test datasets
├── temp/                   # Temporary files
└── install_dependencies.R    # Package installation script
```

## Technology Stack

### Core R Packages

**Data Processing & Analysis:**
- `edgeR`, `DESeq2`, `limma` - Differential expression analysis
- `WGCNA` - Weighted gene co-expression network analysis
- `sva` - Surrogate variable analysis for hidden covariates
- `DDPNA` - Data-dependent proteomics normalization
- `preprocessCore` - Quantile normalization

**Statistics & Modeling:**
- `lme4`, `nlme` - Mixed effects models for PVCA
- `peer` - Probabilistic estimation of expression residuals (optional, Linux only)

**Visualization:**
- `ggplot2` (version 3.5.2), `ggpubr`, `ggrepel` - Publication-quality plots
- `corrplot` - Correlation matrix visualization
- `pheatmap` - Heatmaps

**Functional Enrichment:**
- `clusterProfiler` - GO/KEGG enrichment
- `org.Hs.eg.db`, `org.Mm.eg.db`, `org.Rn.eg.db` - Organism databases

**Utility:**
- `optparse` - Command-line argument parsing
- `data.table` - Fast data manipulation
- `BiocManager` - Bioconductor package management
- `rtracklayer` - GTF file parsing

### Installation

Run the `install_dependencies.R` script to install all dependencies:

```r
source("install_dependencies.R")
```

Or install manually:
```r
install.packages(c("optparse", "DDPNA", "data.table", "remotes", "ggplot2", 
                   "ggpubr", "ggrepel", "doParallel", "lme4", "nlme", 
                   "foreach", "devtools", "WGCNA", "corrplot"))
BiocManager::install(c("Biostrings", "sva", "rtracklayer", "DESeq2", 
                       "limma", "edgeR", "vsn", "impute", "Biobase",
                       "org.Hs.eg.db", "org.Mm.eg.db", "org.Rn.eg.db",
                       "clusterProfiler"))
```

## Usage Patterns

### Two Interface Modes

Each step provides two versions:

1. **`_cmd` files**: Command-line interface using `optparse` for scripting/automation
2. **`_R` files**: R-script interface for interactive analysis with variable injection

### Typical Workflow

#### For RNA-seq Data:

```r
# Step 1: Preprocessing
quantfile <- "counts.RData"  # Contains 'count' and 'tpm' matrices
covarfile <- "metadata.RData"  # Contains 'covar' data frame
output <- "my_analysis"
grpname <- "group"
norm.method <- "TMM"
step <- "234"  # 2=outlier removal, 3=gene filter, 4=normalization
source("step1_countTPM_R")

# Step 2: Covariate adjustment
quantfile <- "my_analysis_filter_TMM.RData"
keepvar <- "diagnose"  # Variables to preserve
method <- "sva"  # Hidden covariate method: sva/peer/pca
adjcov <- "obs&hid"  # Adjust both observed and hidden
evaluate_method <- "pvca&pca"  # Evaluation methods
source("step2_covariate_adj_R")

# Step 3: Differential expression
input <- "my_analysis_adj_final.RData"
groupA <- "control"
groupB <- "treatment"
log2FCcut <- 1
adjPcut <- 0.05
source("step3_DEG_classic_R")
```

#### For MS Proteomics Data:

```r
# Step 0: Extract from raw output
quantfile <- "protein_groups.txt"  # MaxQuant/Proteome Discoverer output
IDinfCol <- "Protein.IDs"  # Column for protein IDs
output <- "proteomics_analysis"
source("step0_MS_R")

# Step 1: Preprocessing
quantfile <- "proteomics_analysis_raw.RData"
step <- "1234"  # 1=norm, 2=outlier, 3=filter, 4=impute
sdout <- 3  # Outlier threshold
maxNA <- 0.8  # Max missing rate
source("step1_MS_R")
```

### Command Line Usage

```bash
# Step 0 - MS data extraction
Rscript step0_MS_cmd --quantfile="protein_groups.txt" \
                     --output="analysis" \
                     --IDinfCol="Protein.IDs"

# Step 1 - RNA-seq preprocessing
Rscript step1_rna_cmd --quantfile="counts.RData" \
                      --covarfile="metadata.RData" \
                      --output="rna_analysis" \
                      --grpname="group" \
                      --step="234" \
                      --norm.method="TMM"

# Step 2 - Covariate adjustment
Rscript step2_covariate_adj_cmd --quantfile="rna_analysis_filter_TMM.RData" \
                                --covarfile="metadata.RData" \
                                --output="adj_analysis" \
                                --keepvar="diagnose" \
                                --method="sva" \
                                --adjcov="obs&hid"
```

## Code Organization

### Main Modules

1. **`preprocess_func.R`** - Core utility library containing:
   - `pvca_evaluate()` / `pvca_plot()` - Principal variance component analysis
   - `pca_dat()` / `pcaplot()` / `pca_plot_var()` - PCA analysis and visualization
   - `runElbow()` / `runBE()` - Optimal PC number determination
   - `BIC_evaluate()` - Bayesian information criterion evaluation
   - `Outlier_Detect()` - Sample outlier detection using IAC/bicor
   - `densityplot()` - Distribution visualization with Shapiro test

2. **Step 1 Scripts** (`step1_*`) - Data type-specific preprocessing:
   - `step1_rna`: RNA-seq count data (CPM → filter → TMM/quantile norm)
   - `step1_countTPM`: RNA-seq with both count and TPM
   - `step1_MS`: Proteomics data (total intensity norm → KNN imputation)

3. **Step 2 Scripts** (`step2_*`) - Covariate handling:
   - `step2_covariate_adj`: General covariate adjustment
   - `step2_batch_adj`: Batch-specific adjustment
   - `step2_cov4QTL`: QTL-specific covariate preparation

4. **Step 3 Scripts** (`step3_*`) - Statistical analysis:
   - `step3_DEG_classic`: t-test/wilcoxon for normalized data
   - `DEG_troika`: limma/edgeR/DESeq2 for count data
   - `step3_WGCNA`: Co-expression network analysis

5. **Step 4 Scripts** (`step4_*`) - Functional analysis:
   - `func_enrich()`: GO/KEGG over-representation
   - `gsea_enrich()`: Gene set enrichment analysis
   - `gsea_1_step()`: Automated GSEA pipeline

### Data Flow

```
Raw Data → Step 0 (Extract) → Step 1 (Preprocess) → Step 2 (Adjust) → Step 3 (Analyze) → Step 4 (Enrich)
     ↓                              ↓                       ↓                  ↓                ↓
  .txt/.csv                  _raw.RData          _adj_obs.RData      _DEG.RData     enrichment
  (MaxQuant,                → _rmOut.RData       _adj_final.RData    _stat.RData    results
  counts, etc.)             → _filter.RData
                            → _norm.RData
```

## Development Conventions

### Coding Style

1. **Variable Naming**:
   - `opt$*` for option parameters
   - `quant` for expression matrix (features × samples)
   - `covar` for sample metadata
   - `IDinf` for feature annotation
   - `loginfo` for run logging

2. **File Naming Conventions**:
   - Output files follow: `{prefix}_{step}.{ext}`
   - RData files preserve all objects with `list=otherobj`
   - PDF plots saved to `plot/` directory

3. **Error Handling**:
   - Extensive `try()` wrapping for robustness
   - `inherits(x, "try-error")` for error detection
   - `stop()` with informative messages

4. **Version Control**:
   - Version and date headers in each script
   - Change log comments (e.g., `#250717 fixed limma log2FC`)

### Input/Output Standards

**Standard Objects in RData Files:**
- `quant` - Normalized expression matrix (genes × samples)
- `covar` - Sample covariates data frame (samples × variables)
- `IDinf` - Feature annotation data frame
- `loginfo` - Processing history list

**Covariate Matching Logic:**
Scripts automatically match `quant` columns to `covar` rows by:
1. Checking column names of `quant` against row names of `covar`
2. Checking column names of `quant` against covariate columns
3. Reordering to ensure alignment

## Testing

Test data is provided in `testdata/` directory:

- `count.txt`, `countTPM.RData` - RNA-seq test data
- `pro.txt`, `phos_data.txt` - Proteomics/phosphoproteomics data
- `*_covar.txt` - Sample metadata
- `测试脚本cmd版本.R`, `测试脚本R版本.R` - Test scripts

Run test analysis:

```r
setwd("testdata")
# Follow test script examples in testdata directory
```

## Key Parameters Reference

### Step 1 Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `step` | Steps to run (e.g., "234") | "234" |
| `sdout` | Outlier SD threshold | 3 |
| `min.value` | Minimum expression threshold | 1 (RNA), 0.1 (countTPM) |
| `norm.method` | Normalization method | "TMM" |
| `grpname` | Grouping variable | "group" |

### Step 2 Parameters

| Parameter | Description | Options |
|-----------|-------------|---------|
| `method` | Hidden covariate method | "sva", "peer", "pca" |
| `adjcov` | Adjustment mode | "obs&hid", "obs", "hid", "obs+hid" |
| `keepvar` | Variables to preserve | comma-separated |
| `svamethod` | SVA number estimation | "be", "leek" |
| `nPCAchoose` | PC selection method | "BE", "Elbow" |
| `evaluate_method` | Evaluation metrics | "pvca", "pca", "BIC", combinations |

### Step 3 DEG Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `log2FCcut` | Log2 fold change cutoff | 1 |
| `adjPcut` | Adjusted p-value cutoff | 0.05 |
| `do_ttest` | Also run t-test (large samples) | FALSE |
| `varequal` | Equal variance assumption | TRUE |

## Security Considerations

1. **No external network calls** - All analysis is local
2. **File path validation** - Scripts check file existence before loading
3. **No persistent credentials** - No authentication mechanisms
4. **Temporary files** - Written to `temp/` directory

## Troubleshooting

### Common Issues

1. **Missing packages**: Run `install_dependencies.R` to install dependencies
2. **Outlier detection fails**: Check `maxNA` parameter
3. **SVA fails**: Ensure `keepvar` is specified and has variation
4. **PEER unavailable**: Install from GitHub or use `method="pca"` instead

### Log Information

Each step appends to `loginfo` object containing:
- `record` - Human-readable processing summary
- `input` / `output` - File paths
- `packageVersion` - Package versions used
- `sessionInfo` - Full R session info
- `parameter` - All parameters used
- `running_time` - Timestamp

## Version History

- **v2.0** (July 2025) - Current version with unified interface
- Added PCA-based hidden covariate method
- Added BIC evaluation
- Improved INT transformation for QTL
- Enhanced outlier detection visualization

