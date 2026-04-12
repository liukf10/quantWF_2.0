path = "D:/func/quantWF_2.0/"
setwd(path) #根据该脚本位置设置路径
if(!dir.exists("temp")) dir.create("temp") #创建临时文件夹
setwd("temp")
#MS step0+1 质谱数据提取 + 质谱数据归一化######
###step0 质谱定量文件提取  ################
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
###step1 质谱数据预处理 ################
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


