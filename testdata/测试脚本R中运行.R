path = "D:/func/quantWF_2.0/"
setwd(path) #根据该脚本位置设置路径
if(!dir.exists("temp")) dir.create("temp") #创建临时文件夹
setwd("temp")

loaded_packages <- grep("^package:", search(), value = TRUE)
loaded_packages <- gsub("^package:", "", loaded_packages)
loaded_packages <- setdiff(loaded_packages, c("base","stats","graphics","grDevices","utils","datasets","methods"))

##使用quantWF需要安装的包###########

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

##批次校正 根据是否有不同批次数据以及样本是否存在批次效应来选择 #########
##包含quant(gene(行）X样本（列）的定量数据
##包含covar 里面包含样本信息（性别,年龄，分组等协变量，还必须包含批次信息）
#支持的校正方法为sva中的combat和limma中的removeBatchEffect
#combat方法有combat（一般是log2之后偏向正态的数据）和针对测序数据combat_seq(需要全是正值，不要log2处理)
#结果中对象为 quant, covar(去除批次列), covar_batch(包含批次列)，以及其他quantfile里面包含的对象
#输出结果为'output'信息_batch_remove.RData 未设置output为data_batch_remove.RData
#数据评估
#校正结果评价方法:PCA,MDS查看样本分布(协变量对其影响); PVCA,BIC查看协变量解释变异程度
#PCA绘制PCA图(且保留PCA数据（可用于探索PCA和基因或者表型信息） output_batch_adj_pcaplot.pdf  和 output_batch_adj_pca.RData
#MDS绘制MDS图只画批次信息且保留PCA数据 output_batch_adj_MDSplot.pdf  和 output_batch_adj_mds.RData
#PCVA绘制PVCA图且保留PVCA数据 output_batch_adj_pvca.pdf  和 output_batch_adj_PVCA.RData
#BIC绘制BIC图且保留BIC数据 output_batch_adj_bic.pdf  和 output_batch_adj_BIC.RData
#绘制不同处理步骤的样本箱线图和基因密度图  output_batch_adj_sample&gene distribution.pdf
rm(list=ls()[!ls() %in% "path"])
#需要校正的定量数据，可以是csv文件(read.csv) 或者是RData 需要有quant对象
quantfile = paste0(path,"testdata/test_for_batch_adj.RData");  #必须提供
#covarfile   协变量文件, txt数据，可以作为covar对象存在quantfile RData
output = "batchtest"; #输出文件的前缀
funcdir = path; #preprocess_func.R路径 #必须提供
#需要保留的协变量（批次校正同时保留这些变量信息)   可以不设置则不保留任何协变量
keepvar = c("Group","Sex");  #covar中的某一列列名
batch_name = "Study";#需要校正的批次名称（covar中的某一列列名） #必须提供
#combat, combat_seq用于测序数据(未对数化，标准化处理的),limma
#method = "combat";  默认是combat 可以不设置
#只需要包含方法对应的字符即可 MDS,PCA,PVCA,BIC会各自存一个RData
evaluate_method = "pca&pvca&mds&bic" #默认是pca&pvca
#BIC_cores 默认30, BIC计算核数
source(paste0(path,"step2_batch_adj_R"))

##MS step0+1 质谱数据提取 + 质谱数据归一化######
#step0 质谱定量文件提取 
#提取IDinf 基于IDinfCol信息
#提取定量数据（基于列名包含intensity或abundance(首字母大小写均可识别)
#将IDinf的第一列的重复值添加_1,_2后缀，作为ID名注释到定量数据
#如果提供covarfile将会和quant样本顺序对齐并作为covar存储入RData
#结果中对象为IDinf,quant, covar(可能有)
#输出结果为output_raw.RData 未设置output为data_raw.RData
rm(list=ls()[!ls() %in% "path"])
#定量文本文件,支持fread读取   磷酸化数据
quantfile = paste0(path,"testdata/phos_data.txt");
#covarfile (可不提供) 协变量文件, txt数据(支持read.table读取)，可以作为covar对象存在quantfile RData
output = "phos" #输出文件的前缀
#IDinfCol 是把哪些列指定为ID的注释信息（基因代谢物或者分子量等信息）
IDinfCol = 1:7 #未指定该参数将把所有非intensity/abundance列设置为ID信息
#NAvalue #缺失值 默认为NA，如果设置为0，则会把定量中数值为0的值改为NA  
source(paste0(path,"step0_MS_R"))

#step1 质谱数据预处理
#第1步 纵向归一化（每列总和值调为一致）
#第2步 移除异常样本 采用曼哈顿距离 存为output_rmOut.RData 
#异常样本检测会绘制图 output_outliersample.pdf
#第3步 过滤高缺失率ID  存为output_rmOut_filter.RData 
#第4步 KNN缺失值填充   存为output_rmOut_impute.RData 
#结果对象为 quant, IDinf,covar(可能有), 以及其他quantfile里面包含的对象
#所存的quant均为log2结果，第四步的结果是先做了横向归一化后的log2值
#数据评估
#绘制不同处理步骤的样本箱线图（log2)和基因密度图（原值） output_sample&gene distribution.pdf
#有covar 则会绘制PCA图(log2的)且保留PCA数据（可用于探索PCA和基因或者表型信息） output_step1_pcaplot.pdf  和 output_step1_pca.RData

rm(list=ls()[!ls() %in% "path"])
#需要处理的定量数据，可以是csv文件(read.csv) 或者是RData 需要有quant和IDinf对象
quantfile = "phos_raw.RData" #必须提供
#covarfile (可不提供) 协变量文件, txt数据(支持read.table读取)，可以作为covar对象存在quantfile RData
output = "phos";  #输出文件的前缀
funcdir = path;  #preprocess_func.R路径  #必须提供
#sdout = 3 默认是方差大于3作为离群样本
maxNA = 0.5 #保留缺失率为多少的样本 默认是0.8 即位80%缺失样本保留
step = 1234 #执行步骤 可以只选其中任意1步或者2步 默认1234
source(paste0(path,"step1_MS_R"))


#蛋白组数据 step1 质谱数据预处理 和前面描述一样
rm(list=ls()[!ls() %in% "path"])
quantfile = paste0(path,"testdata/PRO_data.RData")
output = "pro";
funcdir = path; 
sdout = 2
step = 12
source(paste0(path,"step1_MS_R"))

rm(list=ls()[!ls() %in% "path"])
quantfile = paste0(path,"testdata/PRO_data.RData")
output = "pro";
funcdir = path; 
sdout = 2; maxNA = 0.3
step = 1234 
source(paste0(path,"step1_MS_R"))


##step1 转录组count数据预处理 #############
#第1步 纵向归一化（CPM校正法）
#第2步 移除异常样本 采用曼哈顿距离 存为output_rmOut.RData 
#异常样本检测会绘制图 output_outliersample.pdf
#第3步 过滤低表达基因  存为output_rmOut_filter.RData 
#第4步 TMM或quantile标准化   存为output_rmOut_TMM.RData or output_rmOut_quantile.RData
#同时存log2数据 存为output_log2_TMM.RData or output_log2_quantile.RData
#第5步 INT转换 正态转换 通常用于QTL分析 存为output_rmOut_filter_TMM_INT.RData
#结果对象为 quant, covar(可能有), 以及其他quantfile里面包含的对象
#数据评估
#绘制不同处理步骤的样本箱线图（log2)和基因密度图（原值） output_sample&gene distribution.pdf
#有covar 则会绘制PCA图(log2的)且保留PCA数据（可用于探索PCA和基因或者表型信息） output_step1_pcaplot.pdf  和 output_step1_pca.RData
rm(list=ls()[!ls() %in% "path"])
#定量数据，可以是csv文件(read.csv) 或者是RData 需要有quant
quantfile = paste0(path,"testdata/RNA_data.RData") #必须提供
#covarfile (可不提供) 协变量文件, txt数据(支持read.table读取)，可以作为covar对象存在quantfile RData
output = "rna" #输出文件的前缀
funcdir = path;  #preprocess_func.R路径
#sdout = 3 默认是方差大于3作为离群样本
#min.value = 1  默认是1 过滤低表达基因参数 edgeR::filterByExpr  基因表达值小于此参数为低表达值
#min.total  默认是3*min.value 过滤低表达基因参数 edgeR::filterByExpr 基因在所有样本中的表达之和小于此参数为低表达
#grpname 可选 分组信息 过滤低表达基因参数 edgeR::filterByExpr 最小组的样本数作为n
#grpname 提供时需要有covar,且该值与covar中的某一列名对应
#低表达值比例不超过70%样本数n的基因保留
#norm.method  TMM 或者 quantile 默认是TMM
step = 12345  #执行步骤 可以只选其中任意1步或者2步 默认234
source(paste0(path,"step1_rna_R"))

#包含grpname的测试
rm(list=ls()[!ls() %in% "path"])
quantfile = paste0(path,"testdata/RNA_data.RData")
output = "rna"
funcdir = path; 
grpname = "diagnose" #grpname参数是做低表达基因过滤的时候是否考虑分组信息 
step = 1234
source(paste0(path,"step1_rna_R"))

##step2 协变量校正 主要用于差异分析#########
##包含quant(gene(行）X样本（列）的定量数据
##包含covar 里面包含样本信息（性别,年龄，分组等协变量）
#covar中协变量不能有缺失值,不能是单一值
#校正方法为线性回归
#隐藏协变量提取方法为sva,peer或pca
#第一步，已知协变量校正（仅针对adjcov为obs和obs+hid情况)
#adjcov为obs 存储为output_adj_final.RData; adjcov为obs+hid 存储为 output_adj_obs.RData
#第二步，未知协变量发现和校正（仅针对adjcov为hid和obs+hid情况) 存储为output_adj_final.RData
#method为peer 存碎石图 peer_variance_factor_output.pdf
#method为pca 存碎石图 PCA_variance_factor_output.pdf
#结果中对象为 quant, covar, covar_hid(adjcov有hid才会存),以及其他quantfile里面包含的对象
#数据评估
#协变量之间的相关性图 output_{method}_cov_cor.pdf
#校正结果评价方法:PCA查看样本分布(协变量对其影响); PVCA,BIC查看协变量解释变异程度
#PCA绘制PCA图(且保留PCA数据（可用于探索PCA和基因或者表型信息） output-adj_{adjcov}_pcaplot.pdf  和 output-adj_{adjcov}_pca.RData
#PCVA绘制PVCA图且保留PVCA数据 output-adj_{adjcov}_pvca.pdf  和 output-adj_{adjcov}_PVCA.RData #使用已知协变量
#BIC绘制BIC图且保留BIC数据 output-adj_{adjcov}_bic.pdf  和 output-adj_{adjcov}_BIC.RData #使用已知协变量
#绘制不同处理步骤的样本箱线图和基因密度图  output-adj_{adjcov}_sample&gene distribution.pdf
rm(list=ls()[!ls() %in% "path"])
quantfile = paste0(path,"testdata/test_for_batch_adj.RData");
#covarfile  协变量文件, txt数据，可以作为covar对象存在quantfile RData
output = "svatest"; #输出文件的前缀
funcdir = path; #preprocess_func.R路径 #必须提供
keepvar = c("Group","Sex"); #需要保留的协变量, covar中的其他协变量将被校正 默认为diagnose
#如果不需要保留已知协变量 需要设置为NULL
#vartype #协变量特征, 默认是value, 分组数据不分排序，如果设置为factor,则会自动将非因子字符数据改为因子
#adjcov 协变量校正方法 obs 只校正已知协变量, hid 只校正未知协变量 obs+hid先校正已知再校正未知 obs&hid 一起校正
#method = "sva"; #未知协变量提取方法,默认为sva,可选sva,peer,pca方法
#svamethod #sva方法中协变量发现方法 默认be, 可选be或leek
#nFactor  #peer方法所选factor个数 默认会根据样本量选择数量（基于GTEx建议）
#nPCAchoose #PCA方法隐藏斜变量数量选择 默认BE 可选BE 或Elbow
#只需要包含方法对应的字符即可 MDS,PCA,PVCA,BIC会各自存一个RData
evaluate_method = "pca&pvca&bic" #默认是pca&pvca
#BIC_cores 默认30, BIC计算核数
source(paste0(path,"step2_covariate_adj_R"))


rm(list=ls()[!ls() %in% "path"])
quantfile = paste0(path,"testdata/test4cov4QTL.RData");
output = "pcatest";
funcdir = path;
keepvar = "age"
method = "pca"; #使用PCA方法
evaluate_method = "pca&pvca&bic"
source(paste0(path,"step2_covariate_adj_R"))

##step1+2 QTL数据准备 INT norm #########
#step1 使用count和TPM来过滤, follow GTEX流程
#第1步 去除全是0的基因
#第2步 移除异常样本 采用曼哈顿距离 存为output_rmOut.RData 
#异常样本检测会绘制图 output_outliersample.pdf
#第3步 过滤低表达基因  存为output_rmOut_filter.RData 
#第4步 TMM或quantile标准化   存为output_rmOut_TMM.RData or output_rmOut_quantile.RData
#同时存log2数据 存为output_log2_TMM.RData or output_log2_quantile.RData
#第5步 INT转换 正态转换 通常用于QTL分析 存为output_rmOut_filter_TMM_INT.RData
#绘制不同处理步骤的样本箱线图（log2)和基因密度图（原值） output_sample&gene distribution.pdf
#有covar 则会绘制PCA图(log2的)且保留PCA数据（可用于探索PCA和基因或者表型信息） output_step1_pcaplot.pdf  和 output_step1_pca.RData
rm(list=ls()[!ls() %in% "path"])
#定量数据， RData 需要有count和TPM 2个对象
quantfile = paste0(path,"testdata/countTPM.RData") #必须提供
#covarfile 协变量文件, txt数据(支持read.table读取)，可以作为covar对象存在quantfile RData
output = "test4QTL" #输出文件的前缀
funcdir = path;  #preprocess_func.R路径
#min.value  默认是0.1  过滤低表达基因参数,tpm > 0.1 为低表达值
#过滤低表达基因参数，count > 6 为低表达值(此参数默认 不可设置)
#ratio 默认是0.25   过滤低表达基因参数,非低表达值占比大于0.25
#grpname 可选 分组信息 过滤低表达基因参数 有至少一组样本满足 非低表达值占比大于0.25
#grpname 提供时需要有covar,且该值与covar中的某一列名对应
#norm.method  TMM或者quantile 默认是TMM
step = 12345 #执行步骤 可以只选其中任意1步或者2步 默认12345
source(paste0(path,"step1_countTPM_R"))

#step2 提取隐藏协变量 为QTL分析准备数据
##包含quant(gene(行）X样本（列）的定量数据
##包含covar 里面包含样本信息（性别,年龄，分组等协变量）
##包含genotype_pc 基因组主成分（只保留要校正的）
#包含IDinf
#隐藏协变量提取方法为pca和peer 默认是pca
#搜索隐藏协变量，将covar,genotype_pc和隐藏协变量合并
#提取IDinf和定量数据合并
#method为peer 存碎石图 peer_variance_factor_output.pdf
#method为pca 存碎石图 PCA_variance_factor_output.pdf
#存储协变量文件为 output_peer_covariate4QTL.txt或output_pca_covariate4QTL.txt
#tensorQTL 存储bed文件为 output.tss.tensor.bed
#qtltools 存储bed文件为output.tss.bed  和 output.genebody.bed
#数据评估
#协变量之间的相关性图 output_{method}_cov_cor.pdf
#校正结果评价方法:PCA,MDS查看样本分布(协变量对其影响); PVCA,BIC查看协变量解释变异程度
#PCVA绘制PVCA图且保留PVCA数据 output-adj_{adjcov}_pvca.pdf 和 output-adj_{adjcov}_PVCA.RData #使用所有协变量
#BIC绘制BIC图且保留BIC数据 output-adj_{adjcov}_bic.pdf 和 output-adj_{adjcov}_BIC.RData #使用所有协变量
rm(list=ls()[!ls() %in% "path"])
quantfile = "test4QTL_rmOut_filter_TMM_INT.RData"  #必须提供
covarfile = paste0(path,"testdata/countTPM_covar.txt") #协变量文件, txt数据(支持read.table读取)，可以作为covar对象存在quantfile RData
genotype_pcfile = paste0(path,"testdata/countTPM_genoPC.txt")#基因组主成分 txt数据(支持read.table读取)，genotype_pc存在quantfile RData
#IDfile可以是gtf文件，txt文件和RData 需要有ID,seqnames,start,end,gene_id,strand列 ID列和quant的行名对应
IDfile = paste0(path,"testdata/countTPM_IDinf.txt")#IDinf信息 txt数据(支持read.table读取)，genotype_pc存在quantfile RData
output = "test4QTL"; #输出文件的前缀
funcdir = path;  #preprocess_func.R路径
#method = "pca"; #未知协变量提取方法,默认为pca,可选peer,pca方法
#nFactor  #peer方法所选factor个数 默认会根据样本量选择数量（基于GTEx建议）
#nPCAchoose #PCA方法隐藏斜变量数量选择 默认BE 可选BE 或Elbow
#只需要包含方法对应的字符即可 MDS,PCA,PVCA,BIC会各自存一个RData
evaluate_method = "pca&pvca" #默认是pca&pvca
#BIC_cores 默认30, BIC计算核数
tools = "both" #下游QTL分析工具, tensorQTL,QTLtools,TSS.QTLtools,GB.QTLtools,all.QTLtools,both
#tensorQTL 只输出tensorQTL文件, QTLtools和TSS.QTLtools输出QTLtools标准文件，GB.QTLtools输出genebodyQTL文件, all.QTLtools输出TSS和genebody文件，all全输出
source(paste0(path,"step2_cov4QTL_R"))


rm(list=ls()[!ls() %in% "path"])
quantfile = paste0(path,"testdata/test4cov4QTL.RData");
output = "cov4QTLtest";
funcdir = path;
method = "pca";
BIC_cores = 2;
tools = "both"
source(paste0(path,"step2_cov4QTL_R"))

#指定genotype_pcfile和IDfile
rm(list=ls()[!ls() %in% "path"])
quantfile = paste0(path,"testdata/test4cov4QTL.RData");
genotype_pcfile = paste0(path,"testdata/test4cov4QTL_genotype_pc.txt")
IDfile = paste0(path,"testdata/test4cov4QTL_IDinf.txt")
output = "cov4QTLtest";
funcdir = path;
method = "pca";
BIC_cores = 2;
tools = "both"
source(paste0(path,"step2_cov4QTL_R"))


##step3 WGCNA ###########
rm(list=ls()[!ls() %in% "path"])
input = paste0(path,"testdata/svatest_adj_final.RData")
output = "WGNCATEST"  #输出文件的前缀
funcdir = path
#onlytune #默认是FALSE 只测试不同参数组合的结果，不返回最终结果
#tuning  #默认是FALSE  会测试不同power和其他参数组合下的结果，并选择一个合适参数运行最终结果
#WGCNA相关参数 见WGCNA::blockwiseModules参数介绍   https://edo98811.github.io/WGCNA_official_documentation/
#networkType  #默认是"unsigned" 网络类型 是否优选正相关
#TOMType   #默认是"unsigned"  TOM矩阵类型 是否优选正相关
#power #默认是NA  不设置时会使用pickSoftThreshold选取最佳power 
#deepSplit #默认是2  0-4 值越大 分割模块越多
#minModuleSize  #默认是50  最小模块的基因数
#minKMEtoStay  #默认是0.3  较低KME值会被移除出模块
#detectCutHeight #默认是0.995  修建树参数
#reassignThreshold #默认是1e-4
#mergeCutHeight #默认是0.15
#MaxMod0ratio #默认是0.5 #只有tuning 为TRUE才有用，会要求Mod0基因比例不超过这个值
source(paste0(funcdir,"step3_WGCNA"))

#DEG classic 使用t检验或者wilcox检验分析差异基因（需要完成step2之后的结果）########
#需要有quant,covar，使用FDR校正
#大样本做wilcox,只有单组样本小于10的时候会强制做ttest
#结果存为output_DEG_classic.RData(差异基因)和output_diffstat_classic.RData（所有基因差异统计结果）
rm(list=ls()[!ls() %in% "path"])
input = paste0(path,"testdata/svatest_adj_final.RData")
output = "差异分析"  #输出文件的前缀
grpname = "Group" #进行差异分析的名称（covar中的列名） 默认设置为group
#grpname列中只有2个组 可以不指定groupA,B 
#groupA  大于2组需指定  可以多个
#groupB  大于2组需指定  可以多个
#do_ttest 是否做t检验 默认为FALSE, 单组样本数小于10，则只做ttest 
#varequal 是否方差齐性 默认为TRUE 针对t.test
#log2FCcut 默认为1
#adjPcut 默认为0.05
source(paste0(path,"step3_DEG_classic_R"))

#DEG troika DEseq2,limma,edgeR的差异分析（需要用count数据）######
#需要一个count,一个covar, groupA 和 groupB
rm(list=ls()[!ls() %in% "path"])
counts <- read.delim(paste0(path,"testdata/count.txt"),row.names = 1)
covar <- gsub("[0-9]","",colnames(counts))
groupA = c("ctl","ctl"); groupB = c("exp","expa");
source("D:/func/DEG_troika.R")
