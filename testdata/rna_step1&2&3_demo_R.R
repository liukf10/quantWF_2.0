path = "D:/func/quantWF_2.0/"
setwd(path) #根据该脚本位置设置路径
if(!dir.exists("temp")) dir.create("temp") #创建临时文件夹
setwd("temp")

##step1 转录组count数据预处理 #############
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
step = 1234  #执行步骤 可以只选其中任意1步或者2步 默认234
source(paste0(path,"step1_rna_R"))

##step2 协变量校正 #########
rm(list=ls()[!ls() %in% "path"])
quantfile = paste0(path,"temp/rna_log_TMM.RData");
#covarfile  协变量文件, txt数据，可以作为covar对象存在quantfile RData
output = "rna"; #输出文件的前缀
funcdir = path; #preprocess_func.R路径 #必须提供
keepvar = "diagnose"; #需要保留的协变量, covar中的其他协变量将被校正 默认为diagnose
evaluate_method = "pca&pvca&bic" #默认是pca&pvca
#BIC_cores 默认30, BIC计算核数
source(paste0(path,"step2_covariate_adj_R"))

##step3 差异分析 #########
rm(list=ls()[!ls() %in% "path"])
input = paste0(path,"temp/rna_adj_final.RData");
output = "rna"  #输出文件的前缀
funcdir = path;
grpname = "diagnose" #进行差异分析的名称（covar中的列名） 默认设置为group
#grpname列中只有2个组 可以不指定groupA,B 
#groupA  大于2组需指定  可以多个
#groupB  大于2组需指定  可以多个
#do_ttest 是否做t检验 默认为FALSE, 单组样本数小于10，则只做ttest 
#varequal 是否方差齐性 默认为TRUE 针对t.test
#log2FCcut 默认为1
#adjPcut 默认为0.05
source(paste0(path,"step3_DEG_classic_R"))

# 分步执行step1 ############################
###分步执行step1 按顺序 #############
rm(list=ls()[!ls() %in% "path"])
quantfile = paste0(path,"testdata/RNA_data.RData") #必须提供
output = "rna2" #输出文件的前缀
funcdir = path;  #preprocess_func.R路径
step = 12  #执行步骤 可以只选其中任意1步或者2步 默认234
source(paste0(path,"step1_rna_R"))


rm(list=ls()[!ls() %in% "path"])
quantfile = "rna2_rmOut.RData"
output = "rna2"
funcdir = path; 
grpname = "diagnose" #grpname参数是做低表达基因过滤的时候是否考虑分组信息 
step = 3
source(paste0(path,"step1_rna_R"))

rm(list=ls()[!ls() %in% "path"])
quantfile = "rna2_filter.RData"
output = "rna2"
funcdir = path; 
step = 4
source(paste0(path,"step1_rna_R"))

rm(list=ls()[!ls() %in% "path"])
quantfile = "rna2_TMM.RData"
output = "rna2"
funcdir = path; 
step = 5
source(paste0(path,"step1_rna_R"))
###分步执行step1 不按顺序#############
rm(list=ls()[!ls() %in% "path"])
quantfile = paste0(path,"testdata/RNA_data.RData") #必须提供
output = "rna3" #输出文件的前缀
funcdir = path;  #preprocess_func.R路径
step = 13  #执行步骤 可以只选其中任意1步或者2步 默认234
source(paste0(path,"step1_rna_R"))


rm(list=ls()[!ls() %in% "path"])
quantfile = "rna3_filter.RData"
output = "rna3"
funcdir = path; 
step = 2
source(paste0(path,"step1_rna_R"))

rm(list=ls()[!ls() %in% "path"])
quantfile = "rna2_rmOut.RData"
output = "rna2"
funcdir = path; 
step = 4
source(paste0(path,"step1_rna_R"))
