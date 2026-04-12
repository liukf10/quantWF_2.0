# quantWF_2.0 使用说明文档

本文档提供 quantWF_2.0 框架的详细使用指南。

---

## 📑 目录

1. [入门指南](#-入门指南)
2. [数据准备](#-数据准备)
3. [Step 0: 数据提取](#-step-0-数据提取)
4. [Step 1: 数据预处理](#-step-1-数据预处理)
5. [Step 2: 协变量校正](#-step-2-协变量校正)
6. [Step 3: 统计分析](#-step-3-统计分析)
7. [Step 4: 功能富集](#-step-4-功能富集)
8. [输出文件说明](#-输出文件说明)
9. [故障排除](#-故障排除)

---

## 🚀 入门指南

### 环境配置

```r
# 设置工作目录
setwd("D:/func/quantWF_2.0/")

# 创建临时文件夹存放分析结果
if(!dir.exists("temp")) dir.create("temp")
setwd("temp")

# 加载依赖包（首次使用需安装）
source("../install_dependencies.R")
```

### 基本概念

quantWF_2.0 采用**变量注入**方式运行。在调用脚本前，需要预先定义相应的参数变量：

```r
# 定义参数
quantfile <- "your_data.RData"
output <- "analysis_name"
step <- "234"

# 运行脚本
source("step1_rna_R")
```

---

## 📊 数据准备

### 标准数据格式

#### 1. 定量数据 (quant)
- **格式**: 矩阵，行是基因/蛋白，列是样本
- **类型**: 
  - RNA-seq: 原始count值或CPM值
  - 蛋白质组: 强度值 (Intensity)
- **存储**: RData文件中的 `quant` 对象，或CSV文件

#### 2. 样本信息 (covar)
- **格式**: 数据框，行是样本，列是协变量
- **必需列**: 分组信息（如 group, diagnose）
- **存储**: RData文件中的 `covar` 对象，或CSV/TXT文件

#### 3. 注释信息 (IDinf) - 可选
- **格式**: 数据框，包含ID的注释信息
- **必需列**: `ID` 列与quant的行名对应

### 示例数据结构

```r
# 查看测试数据结构
load("testdata/RNA_data.RData")
head(quant[,1:5])       # 定量数据
head(covar)             # 样本信息
```

---

## 📥 Step 0: 数据提取

用于从质谱分析软件（MaxQuant、Proteome Discoverer等）的输出文件中提取定量数据。

### step0_MS_R

```r
# 参数设置
quantfile <- "testdata/phos_data.txt"   # 输入文件路径
output <- "phos"                         # 输出前缀
IDinfCol <- 1:7                          # ID信息所在列

# 运行
source("step0_MS_R")
```

**输出文件**:
- `phos_raw.RData`: 包含 quant, IDinf, covar(如有)

**参数说明**:

| 参数 | 说明 | 默认值 |
|:---|:---|:---|
| `quantfile` | 输入文件路径 | 必须提供 |
| `output` | 输出前缀 | "data" |
| `IDinfCol` | ID信息列号 | 所有非intensity列 |
| `NAvalue` | 缺失值标记 | NA |

---

## 🔧 Step 1: 数据预处理

根据数据类型选择对应的脚本：

### 1.1 RNA-seq 预处理 (step1_rna_R)

适用于RNA-seq原始count数据。

```r
# 参数设置
quantfile <- "testdata/RNA_data.RData"
covarfile <- "testdata/countTPM_covar.txt"  # 可选
output <- "rna"
funcdir <- "D:/func/quantWF_2.0/"             # preprocess_func.R路径
grpname <- "diagnose"                         # 分组列名，可选
sdout <- 3                                    # 异常样本阈值
min.value <- 1                                # 低表达阈值
norm.method <- "TMM"                          # 标准化方法
step <- "2345"                                # 执行步骤

# 运行
source("step1_rna_R")
```

**Step参数详解**:

| Step | 功能 | 说明 |
|:---:|:---|:---|
| 1 | CPM标准化 | 转换为CPM值（通常不需要单独执行） |
| 2 | 异常样本检测 | 基于曼哈顿距离，移除sdout倍标准差外的样本 |
| 3 | 低表达基因过滤 | 使用edgeR::filterByExpr |
| 4 | TMM/Quantile标准化 | 文库大小校正 |
| 5 | INT转换 | 逆正态转换，用于QTL分析 |

**输出文件**（假设output="rna"，outputdir="./output"）：

| 文件 | 路径 | 说明 |
|:---|:---|:---|
| `rna_rmOut.RData` | `outputdir/` | 移除异常样本后的数据 |
| `rna_rmOut_filter.RData` | `outputdir/` | 过滤低表达基因后的数据 |
| `rna_rmOut_filter_TMM.RData` | `outputdir/` | TMM标准化结果 |
| `rna_log_TMM.RData` | `outputdir/` | log2转换后结果（用于后续分析） |
| `rna_rmOut_filter_TMM_INT.RData` | `outputdir/` | INT转换结果（用于QTL） |
| `rna_report.md` | `outputdir/` | Markdown分析报告 |
| `rna_step1_pca.pdf` | `outputdir/plot/` | PCA图（多页PDF） |
| `rna_step1_pca_*.png` | `outputdir/assets/` | PCA图（PNG，按raw→cpm→filter→norm→INT顺序） |
| `rna_step1_distribution.pdf` | `outputdir/plot/` | 样本和基因分布图（PDF） |
| `rna_step1_box_*.png` | `outputdir/assets/` | 各阶段样本箱线图 |
| `rna_step1_density_*.png` | `outputdir/assets/` | 各阶段基因密度图 |

### 1.2 Count+TPM 预处理 (step1_countTPM_R)

适用于有count和TPM两种数据的场景（如GTEx流程）。

```r
quantfile <- "testdata/countTPM.RData"  # 包含count和TPM两个对象
output <- "test4QTL"
funcdir <- "D:/func/quantWF_2.0/"
min.value <- 0.1                        # TPM阈值
step <- "12345"
source("step1_countTPM_R")
```

**过滤逻辑**:
- count > 6 且 TPM > 0.1 定义为低表达
- 保留非低表达值占比 > 25% 的基因

### 1.3 质谱数据预处理 (step1_MS_R)

```r
quantfile <- "phos_raw.RData"
output <- "phos"
funcdir <- "D:/func/quantWF_2.0/"
sdout <- 3
maxNA <- 0.5                          # 最大允许缺失率
step <- "1234"
source("step1_MS_R")
```

**处理流程**:
1. 总强度归一化（每列和调为一致）
2. 异常样本检测（曼哈顿距离）
3. 高缺失率ID过滤（maxNA参数控制）
4. KNN缺失值填充

---

## 🎛️ Step 2: 协变量校正

### 2.1 一般协变量校正 (step2_covariate_adj_R)

用于差异表达分析前的协变量校正，包括已知协变量和隐藏协变量。

```r
quantfile <- "rna_log2_TMM.RData"
covarfile <- "testdata/countTPM_covar.txt"
output <- "svatest"
funcdir <- "D:/func/quantWF_2.0/"
keepvar <- "group,sex"              # 需要保留的协变量
adjcov <- "obs&hid"                 # 校正模式
method <- "sva"                     # 隐藏协变量方法
svamethod <- "be"                   # SVA估计方法
evaluate_method <- "pca&pvca&bic"   # 评估方法
source("step2_covariate_adj_R")
```

**关键参数详解**:

#### adjcov 校正模式

| 模式 | 说明 | 适用场景 |
|:---|:---|:---|
| `obs` | 只校正已知协变量 | 协变量影响明确 |
| `hid` | 只校正隐藏协变量 | 用于QTL分析 |
| `obs+hid` | 先校正已知，再校正隐藏 | 保护生物学变异 |
| `obs&hid` | 一起校正（推荐） | 差异表达分析 |

#### method 隐藏协变量方法

| 方法 | 特点 | 推荐场景 |
|:---|:---|:---|
| `sva` | 基于代理变量分析 | 有明确生物学分组 |
| `peer` | 概率估计 | 大样本（>100） |
| `pca` | 主成分分析 | 需要快速分析 |

#### 输出文件

| 文件 | 路径 | 说明 |
|:---|:---|:---|
| `*_adj_obs.RData` | `outputdir/` | 已知协变量校正结果（仅obs+hid模式） |
| `*_adj_final.RData` | `outputdir/` | 最终校正结果 |
| `*_report.md` | `outputdir/` | Markdown分析报告 |
| `*_cov_cor.pdf` | `outputdir/plot/` | 协变量相关性热图（PDF） |
| `*_cov_cor.png` | `outputdir/assets/` | 协变量相关性热图（PNG） |
| `*_step2_*_pca.pdf` | `outputdir/plot/` | 校正前后PCA图（PDF，多页） |
| `*_step2_*_pca_*.png` | `outputdir/assets/` | PCA图（PNG，按协变量分组） |
| `*_step2_*_pvca.pdf` | `outputdir/plot/` | PVCA分析图（PDF） |
| `*_step2_*_pvca.png` | `outputdir/assets/` | PVCA分析图（PNG） |
| `*_step2_*_bic.pdf` | `outputdir/plot/` | BIC分析图（PDF） |
| `*_step2_*_bic.png` | `outputdir/assets/` | BIC分析图（PNG） |
| `*_step2_box_*.png` | `outputdir/assets/` | 样本分布箱线图 |
| `*_step2_density_*.png` | `outputdir/assets/` | 基因密度图 |

### 2.2 批次校正 (step2_batch_adj_R)

用于校正明显的批次效应。

```r
quantfile <- "testdata/test_for_batch_adj.RData"
output <- "batchtest"
funcdir <- "D:/func/quantWF_2.0/"
keepvar <- "Group,Sex"
batch_name <- "Study"               # 批次列名
method <- "combat"                  # combat/combat_seq/limma
evaluate_method <- "pca&pvca&mds&bic"
source("step2_batch_adj_R")
```

### 2.3 QTL数据准备 (step2_cov4QTL_R)

为QTL分析准备协变量文件。

```r
quantfile <- "test4QTL_rmOut_filter_TMM_INT.RData"
covarfile <- "testdata/countTPM_covar.txt"
genotype_pcfile <- "testdata/countTPM_genoPC.txt"
IDfile <- "testdata/countTPM_IDinf.txt"
output <- "test4QTL"
funcdir <- "D:/func/quantWF_2.0/"
method <- "pca"
tools <- "both"                     # tensorQTL和QTLtools格式
evaluate_method <- "pca&pvca"
source("step2_cov4QTL_R")
```

**输出文件**:

| 文件 | 路径 | 说明 |
|:---|:---|:---|
| `*_pca_covariate4QTL.txt` | `outputdir/` | 协变量文件 |
| `*.tss.tensor.bed` | `outputdir/` | tensorQTL格式表达矩阵 |
| `*.tss.bed` | `outputdir/` | QTLtools TSS格式 |
| `*_report.md` | `outputdir/` | Markdown分析报告 |
| `PCA_variance_factor_*.pdf` | `outputdir/plot/` | PCA scree plot（PDF） |
| `PCA_variance_factor_*.png` | `outputdir/assets/` | PCA scree plot（PNG） |
| `*_covhid_cor.pdf` | `outputdir/plot/` | 协变量相关性热图（PDF） |
| `*_covhid_cor.png` | `outputdir/assets/` | 协变量相关性热图（PNG） |
| `*_pvca.pdf` | `outputdir/plot/` | PVCA分析图（PDF） |
| `*_pvca.png` | `outputdir/assets/` | PVCA分析图（PNG） |
| `*_bic.pdf` | `outputdir/plot/` | BIC分析图（PDF） |
| `*_bic.png` | `outputdir/assets/` | BIC分析图（PNG） |

---

## 📈 Step 3: 统计分析

### 3.1 经典差异表达分析 (step3_DEG_classic_R)

基于t检验或wilcoxon检验，适用于已标准化的数据（如经过Step 2校正的数据）。

```r
input <- "svatest_adj_final.RData"
output <- "差异分析"
grpname <- "Group"
groupA <- "control"                 # 对照组
groupB <- "treatment"               # 实验组
do_ttest <- FALSE                   # 是否做t检验
varequal <- TRUE                    # 方差齐性假设
log2FCcut <- 1                      # log2FC阈值
adjPcut <- 0.05                     # 校正p值阈值
source("step3_DEG_classic_R")
```

**样本量判断**:
- 单组样本 < 10: 强制使用t检验
- 单组样本 ≥ 10: 使用wilcoxon检验（大样本更稳健）

**输出文件**:
- `*_DEG_classic.RData`: 差异基因列表
- `*_diffstat_classic.RData`: 所有基因的统计结果

### 3.2 Count数据差异分析 (DEG_troika.R)

使用limma-voom、edgeR、DESeq2三种方法分析原始count数据。

```r
counts <- read.delim("testdata/count.txt", row.names = 1)
covar <- data.frame(
  sample = colnames(counts),
  group = gsub("[0-9]", "", colnames(counts))
)
output <- "DEG_results"
grpname <- "group"
groupA <- c("ctl", "ctl")
groupB <- c("exp", "expa")
log2FCcut <- 1
adjPcut <- 0.05
source("DEG_troika.R")
```

**输出文件**:
- `*_DEG_troika.RData`: 三种方法的差异基因
- `*_stat_troika.RData`: 三种方法的完整统计结果

### 3.3 WGCNA分析 (step3_WGCNA)

加权基因共表达网络分析。

```r
input <- "testdata/svatest_adj_final.RData"
output <- "WGNCATEST"
funcdir <- "D:/func/quantWF_2.0/"
tuning <- FALSE                     # 是否参数调优
networkType <- "unsigned"
TOMType <- "unsigned"
power <- NA                         # 自动选择
minModuleSize <- 50
mergeCutHeight <- 0.15
source("step3_WGCNA")
```

---

## 🧬 Step 4: 功能富集

### 4.1 GO/KEGG富集 (step4_GO&KEGGenrich.R)

```r
source("step4_GO&KEGGenrich.R")

# 单组基因富集
genes <- c("TP53", "BRCA1", "BRCA2", "EGFR", "PTEN")
result <- func_enrich(
  gene = genes,
  organism = "human",
  keyType = "SYMBOL",
  filename = "enrich_result"
)

# 查看结果
head(as.data.frame(result$BP))  # 生物学过程
head(as.data.frame(result$KEGG)) # KEGG通路
```

### 4.2 GSEA分析 (step4_GS_GSenrich.R)

```r
source("step4_GO&KEGGenrich.R")

# 基于差异分析结果做GSEA
load("差异分析_diffstat_classic.RData")
stat <- wilcox_result_all[[1]]  # 获取统计结果

result <- gsea_1_step(
  stat = stat,
  id = "SYMBOL",
  log2fc = "log2FC",
  pvalue = "p.value",
  organism = "human",
  keyType = "SYMBOL",
  filename = "gsea_result"
)
```

---

## 📁 输出文件说明

### RData文件命名规则

```
{前缀}_{步骤描述}_{方法}.RData
```

### 标准对象说明

| 对象名 | 类型 | 说明 |
|:---|:---|:---|
| `quant` | matrix | 标准化后的定量矩阵（基因×样本） |
| `covar` | data.frame | 样本协变量信息 |
| `IDinf` | data.frame | 基因/蛋白注释信息 |
| `loginfo` | list | 运行日志和参数记录 |

### 图形输出

图形输出分为两个目录：
- **`plot/`** 目录：保存高分辨率 PDF 文件（适合发表）
- **`assets/`** 目录：保存 PNG 文件（用于 Markdown 报告嵌入）

| 文件名模式 | 格式 | 内容 |
|:---|:---:|:---|
| `*_step1_pca_*.pdf/png` | PDF/PNG | Step 1各阶段PCA图（按 raw→cpm→rm→norm→INT 顺序排列） |
| `*_step1_distribution.pdf/png` | PDF/PNG | 样本箱线图和基因密度图 |
| `*_step1_box_*.png` | PNG | 各阶段样本箱线图（raw, cpm, filter, norm, INT） |
| `*_step1_density_*.png` | PNG | 各阶段基因密度图 |
| `*_outliersample.pdf` | PDF | 异常样本检测图 |
| `*_cov_cor.pdf/png` | PDF/PNG | 协变量相关性热图 |
| `*_step2_*_pca.pdf/png` | PDF/PNG | Step 2校正前后PCA对比 |
| `*_step2_*_pvca.pdf/png` | PDF/PNG | PVCA分析图 |
| `*_step2_*_bic.pdf/png` | PDF/PNG | BIC分析图 |
| `*_step2_box_*.png` | PNG | Step 2各阶段样本分布箱线图 |
| `*_step2_density_*.png` | PNG | Step 2各阶段基因密度图 |

---

## 🔧 故障排除

### 常见问题

#### Q1: 包安装失败
```r
# 尝试使用国内镜像
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
BiocManager::install("package_name", ask = FALSE)
```

#### Q2: 样本匹配失败
```r
# 检查样本名是否一致
colnames(quant)
rownames(covar)
# 确保covar有行名或有一列包含样本名
```

#### Q3: SVA计算失败
- 检查`keepvar`是否正确设置
- 确保保留的变量有足够的变异
- 检查样本量是否足够（建议>10）

#### Q4: PEER无法安装
PEER包只能在Linux系统安装。Windows/macOS用户可以：
- 使用`method="pca"`或`method="sva"`替代
- 使用WSL或Docker运行Linux环境

#### Q5: 内存不足
```r
# 减小BIC计算核心数
BIC_cores <- 10

# 或者跳过BIC评估
evaluate_method <- "pca&pvca"
```

### 调试技巧

```r
# 查看详细日志
verbose <- 1

# 检查中间结果
load("中间结果.RData")
ls()
head(quant)
```

---

## 📚 进一步阅读

- **测试脚本**: `testdata/测试脚本R中运行.R` 包含完整的可运行示例
- **测试脚本**: `testdata/测试脚本cmd运行.R` 包含命令行示例
- **函数文档**: `preprocess_func.R` 头部包含函数说明

---

<div align="center">

**如有问题请联系**: liukefu19@163.com

</div>
