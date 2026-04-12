path = "D:/func/quantWF_2.0/"
setwd(path) #根据该脚本位置设置路径
if(!dir.exists("temp")) dir.create("temp") #创建临时文件夹
setwd("temp")


rm(list=ls()[!ls() %in% "path"])
quantfile = paste0(path,"testdata/test_for_batch_adj.RData");  #必须提供
output = "batchtest"; #输出文件的前缀
funcdir = path; #preprocess_func.R路径 #必须提供
keepvar = c("Group","Sex");  #covar中的某一列列名
batch_name = "Study";#需要校正的批次名称（covar中的某一列列名） #必须提供
evaluate_method = "pca&pvca&mds&bic" #默认是pca&pvca
source(paste0(path,"step2_batch_adj_R"))
