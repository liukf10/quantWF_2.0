
#cmd操作需要
if(!requireNamespace("optparse")) install.packages("optparse")
#step1需要
if(!requireNamespace("DDPNA")) install.packages("DDPNA")
#step0 需要
if(!requireNamespace("data.table")) install.packages("data.table")
#基础包
if(!requireNamespace("BiocManager")) install.packages("BiocManager")
if(!requireNamespace("Biostrings")) BiocManager::install("Biostrings")
#甲基化芯片数据分析需要
#if(!requireNamespace("ChAMP")) BiocManager::install("ChAMP")
#作图需要
if(!requireNamespace("remotes")) install.packages("remotes")
if(!requireNamespace("ggplot2")) remotes::install_version("ggplot2", version = "3.5.2")
if(!requireNamespace("ggpubr")) install.packages("ggpubr")
if(!requireNamespace("ggrepel")) install.packages("ggrepel")
#BIC需要
if(!requireNamespace("doParallel")) install.packages("doParallel")
#step2 批次校正需要
if(!requireNamespace("sva")) BiocManager::install("sva")
if(!requireNamespace("lme4")) install.packages("lme4")
if(!requireNamespace("nlme")) install.packages("nlme")
if(!requireNamespace("foreach")) install.packages("foreach")
if(!requireNamespace("devtools"))  install.packages("devtools")
#if(!requireNamespace("peer"))  devtools::install_github("PMBio/peer") #不太容易安成功

#gtf文件读取需要
if(!requireNamespace("rtracklayer")) BiocManager::install("rtracklayer")
#step1rna和批次校正都需要
if(!requireNamespace("DESeq2")) BiocManager::install("DESeq2")
if(!requireNamespace("limma")) BiocManager::install("limma")
if(!requireNamespace("edgeR")) BiocManager::install("edgeR")
#WGCNA需要
if(!requireNamespace("WGCNA")) install.packages("WGCNA")
#TAMPOR needed
if(!requireNamespace("vsn")) BiocManager::install("vsn")

if(!requireNamespace("corrplot")) install.packages("corrplot")



