#Version: 1.0
#created: Kefu liu(liukefu19@163.com)
#Maintainer: Kefu Liu
#Date: July 30, 2022
#Software: R


#wgcnadata, output,corType,maxBlockSize,sftthreshould = 0.85
library(DDPNA)
library(WGCNA)
objs <- ls()
if(!"wgcnadata" %in% objs) {stop("no wgcnadata.") }
if(!"output" %in% objs) {stop("no output") }
if(!"sftthreshould" %in% objs) {sftthreshould = 0.85;message(" the default sftthreshould: 0.85.") }
if(!"corType" %in% objs) {corType = "bicor"; message(" the default corType: bicor.")}
if(!"maxBlockSize" %in% objs) {maxBlockSize = 20000; message(" the default maxBlockSize: 20000.")}
if(!"minModNum" %in% objs) {minModNum = 10; message(" the default minModNum: 10.")}
if(!"maxModNum" %in% objs) {maxModNum = 50; message(" the default maxModNum: 50.")}
if(!"MaxMod0ratio " %in% objs) {MaxMod0ratio  = 0.5; message(" the default MaxMod0ratio: 0.5.")}

tuning = TRUE;
if(tuning) {
  if(!"sft" %in% objs) {
      sft = SoftThresholdScaleGraph(wgcnadata,filename = output,
                                    corFnc = bicor)
      sft$powerEstimate = min(sft$fitIndices$Power[sft$fitIndices$SFT.R.sq > sftthreshould])
  }
  if(is.null(sft$powerEstimate)) stop("No SFT reach the setting parameter.")
  if(is.infinite(sft$powerEstimate)) stop("No SFT reach the setting parameter.")
  cat("--the best WGCNA construct power is ", sft$powerEstimate,"--") 
  WGCNAadjust <- wgcnatest(wgcnadata, power = sft$powerEstimate,
                           corType = corType, maxBlockSize = maxBlockSize,
                           minModNum = minModNum, maxModNum = maxModNum,
                           MaxMod0ratio = MaxMod0ratio)
  write.csv(WGCNAadjust,file = paste0(output,"_WGCNAadjust.csv"))
}