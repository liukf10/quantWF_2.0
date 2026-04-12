path = "D:/func/quantWF_2.0/"
setwd(path) #根据该脚本位置设置路径
if(!dir.exists("temp")) dir.create("temp") #创建临时文件夹
setwd("temp")
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

##########################
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

