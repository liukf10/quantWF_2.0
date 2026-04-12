path = "D:/func/quantWF_2.0/"
setwd(path) #根据该脚本位置设置路径
if(!dir.exists("temp")) dir.create("temp") #创建临时文件夹
setwd("temp")
###step1 质谱数据预处理 一步执行#############
rm(list=ls()[!ls() %in% "path"])
quantfile = paste0(path,"testdata/PRO_data.RData")
output = "pro";
funcdir = path; 
sdout = 2
maxNA = 0.3
step = 1234
source(paste0(path,"step1_MS_R"))
###step2 质谱数据斜变量校正 需要执行完step1才能运行#######
rm(list=ls()[!ls() %in% "path"])
quantfile = paste0(path,"temp/pro_rmOut_filter_impute.RData");
#covarfile  协变量文件, txt数据，可以作为covar对象存在quantfile RData
output = "pro"; #输出文件的前缀
funcdir = path; #preprocess_func.R路径 #必须提供
keepvar = "group"; #需要保留的协变量, covar中的其他协变量将被校正 默认为diagnose
evaluate_method = "pca&pvca&bic" #默认是pca&pvca
#BIC_cores 默认30, BIC计算核数
source(paste0(path,"step2_covariate_adj_R"))
##分步骤执行 output不一样时候，会把前一步output内容迁移进来 ##############
###step1 质谱数据预处理 分步骤运行 ################
#蛋白组数据 step1 质谱数据预处理 和前面描述一样
rm(list=ls()[!ls() %in% "path"])
quantfile = paste0(path,"testdata/PRO_data.RData")
output = "pro2";
funcdir = path; 
sdout = 2
step = 12
source(paste0(path,"step1_MS_R"))

rm(list=ls()[!ls() %in% "path"])
quantfile = paste0(path,"temp/pro2_rmOut.RData")
output = "pro2_rmOut";
funcdir = path;
step = 3
source(paste0(path,"step1_MS_R"))

rm(list=ls()[!ls() %in% "path"])
quantfile = paste0(path,"temp/pro2_rmOut_filter.RData")
output = "pro2_rmOut_filter";
funcdir = path; 
maxNA = 0.3
step = 4
source(paste0(path,"step1_MS_R"))







##output一样的时候，名字只会是output+本身步骤的后缀########
###step1 质谱数据预处理 分步骤运行 ################
#蛋白组数据 step1 质谱数据预处理 和前面描述一样
rm(list=ls()[!ls() %in% "path"])
quantfile = paste0(path,"testdata/PRO_data.RData")
output = "pro3";
funcdir = path; 
sdout = 2
step = 12
source(paste0(path,"step1_MS_R"))

rm(list=ls()[!ls() %in% "path"])
quantfile = paste0(path,"temp/pro3_rmOut.RData")
output = "pro3";
funcdir = path;
step = 3
source(paste0(path,"step1_MS_R"))

rm(list=ls()[!ls() %in% "path"])
quantfile = paste0(path,"temp/pro3_filter.RData")
output = "pro3";
funcdir = path; 
maxNA = 0.3
step = 4
source(paste0(path,"step1_MS_R"))






##不按顺序执行 先执行filter 再执行步骤2,再执行步骤4########
###step1 质谱数据预处理 分步骤运行 ################
#蛋白组数据 step1 质谱数据预处理 和前面描述一样
rm(list=ls()[!ls() %in% "path"])
quantfile = paste0(path,"testdata/PRO_data.RData")
output = "pro4";
funcdir = path; 
step = 13
source(paste0(path,"step1_MS_R"))

rm(list=ls()[!ls() %in% "path"])
quantfile = paste0(path,"temp/pro4_filter.RData")
output = "pro4";
funcdir = path;
step = 2
source(paste0(path,"step1_MS_R"))

rm(list=ls()[!ls() %in% "path"])
quantfile = paste0(path,"temp/pro4_filter.RData")
output = "pro4";
funcdir = path; 
maxNA = 0.3
step = 4
source(paste0(path,"step1_MS_R"))

