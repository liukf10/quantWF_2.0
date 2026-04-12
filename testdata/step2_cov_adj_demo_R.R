path = "D:/func/quantWF_2.0/"
setwd(path) #根据该脚本位置设置路径
if(!dir.exists("temp")) dir.create("temp") #创建临时文件夹
setwd("temp")
##sva寻找隐藏协变量 #########
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
quantfile = paste0(path,"testdata/test4cov4QTL.RData");
#covarfile  协变量文件, txt数据，可以作为covar对象存在quantfile RData
output = "svatest"; #输出文件的前缀
funcdir = path; #preprocess_func.R路径 #必须提供
keepvar = c("age","sex_check"); #需要保留的协变量, covar中的其他协变量将被校正 默认为diagnose
#如果不需要保留已知协变量 需要设置为NULL
#vartype #协变量特征, 默认是value, 分组数据不分排序，如果设置为factor,则会自动将非因子字符数据改为因子
adjcov = "obs+hid" # obs,只校正已知协变量; hid,只校正未知协变量;obs+hid先校正已知再校正未知;obs&hid,一起校正
#method = "sva"; #未知协变量提取方法,默认为sva,可选sva,peer,pca方法
#svamethod #sva方法中协变量发现方法 默认be, 可选be或leek
#nFactor  #peer方法所选factor个数 默认会根据样本量选择数量（基于GTEx建议）
#nPCAchoose #PCA方法隐藏斜变量数量选择 默认BE 可选BE 或Elbow
#只需要包含方法对应的字符即可 MDS,PCA,PVCA,BIC会各自存一个RData
evaluate_method = "pca&pvca&bic" #默认是pca&pvca
#BIC_cores 默认30, BIC计算核数
source(paste0(path,"step2_covariate_adj_R"))

##PCA寻找隐藏协变量 #########
rm(list=ls()[!ls() %in% "path"])
quantfile = paste0(path,"testdata/test4cov4QTL.RData");
output = "pcatest";
funcdir = path;
keepvar = "age"
method = "pca"; #使用PCA方法
evaluate_method = "pca&pvca&bic"
source(paste0(path,"step2_covariate_adj_R"))

