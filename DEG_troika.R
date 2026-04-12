#Version: 2.1
#created: Kefu liu(liukefu19@163.com)
#Maintainer: Kefu Liu
#Date: March 25, 2026
#Software: R
#20250717 fixed limma and edgeR have ctl/ind log2FC
#20250918 add loginfo and message in last, and fix DESeq2 DEG extract error
#V2.1: 添加outputdir参数支持指定输出目录

#######
#log2FCcut = 1; adjPcut = 0.05;
#counts; covar;grpname = "group"

#说明############
#padj为NA
#Note on p-values set to NA: some values in the results table can be set to NA for one of the following reasons:
#1. 所有样品中为0, baseMean为0,log2FC,p,padj为NA
#2.有一个样品有异常值,p,padj设置成NA. Cooks's distance 来检测outlier.
#3. 如果基因没有通过DESeq2的自动做的independent filtering(有低的mean normalized count),只有padj会设置成NA.

# object check##################
allobj <- ls();
if(!"counts" %in% allobj) stop("No counts obj, see readme document.");
if(!"covar" %in% allobj) stop("No covar obj, see readme document.");
if(!"output" %in% allobj) {
  output = "data";
  message("no output, setting to 'data' as prefix")
}
if(!"outputdir" %in% allobj) {
  outputdir = ".";
  message("no outputdir, setting to current directory")
}
if(!"grpname" %in% allobj) {
  grpname = "group"; message("Setting grpname as default: group.")}
if(is.vector(covar)) {
  if(length(covar) != ncol(counts)) stop("covar length is not equal with counts column number.")
  message("covar is a vector, considered in the same order with counts.")
  covar <- data.frame(samplename = colnames(counts),grpname = covar);
  colnames(covar)[2] = grpname;
}
for (i in 1:(ncol(covar)+1)) {
  if(i == ncol(covar)+1) {
    if(all(colnames(counts) %in% rownames(covar)))
      break else i = 0;
  }
  if(all(colnames(counts) %in% covar[,i]))
    break
}
if(i == 0) stop("Not all sample match with covar.")
if(i == (ncol(covar)+1)) 
  pos <- match(colnames(counts),rownames(covar)) else {
    pos <- match(colnames(counts),covar[,i]);
    covar <- covar[pos,]; covar <- covar[-i]; 
    rownames(covar) <- colnames(counts);
  } #covar & counts are matched
if(any(!grpname %in% colnames(covar)))
  stop("grpname is not in covar, please check it.")
group_list <- covar[,grpname];
group = unique(group_list)


if(!"log2FCcut" %in% allobj) {
  log2FCcut = 1; message("Setting log2FCcut as default: 1.")}
if(!"adjPcut" %in% allobj) {
  adjPcut = 0.05; message("Setting adjPcut as default: 0.05.")}

if(length(group) != 2) {
  if(!"groupA" %in% allobj) stop("No groupA, see readme document.")
  if(!"groupB" %in% allobj) stop("No groupB, see readme document.")
  if(!all(groupA %in% group)) stop("Not all groupA belongs to group.")
  if(!all(groupB %in% group)) stop("Not all groupB belongs to group.")
  if(length(groupA) != length(groupB)) stop("groupA and groupB should have the same length.")
} else {
  if(!"groupA" %in% allobj | !"groupB" %in% allobj) {
    groupA = group[1]; groupB = group[2];
    message(paste0("setting groupA as ",group[1],"; groupB as ",group[2]))
  }
  if(!all(groupA %in% group)) stop("Not all groupA belongs to group.")
  if(!all(groupB %in% group)) stop("Not all groupB belongs to group.")
  if(length(groupA) != length(groupB)) stop("groupA and groupB should have the same length.")
}

###0.package check ########
existpck <- rownames(installed.packages());
if(!"limma" %in% existpck)
  stop("no limma package, please install the R pakcage before run me.")
if(!"edgeR" %in% existpck)
  stop("no edgeR package, please install the R pakcage before run me.")
if(!"DESeq2" %in% existpck)
  stop("no DESeq2 package, please install the R pakcage before run me.")

##limma #####
library(limma)
library(edgeR)
limma_result=list()
limma_result_all=list()

for(i in seq_along(groupB)){
  ctl = which(group_list == groupA[i])
  ind=which(group_list==groupB[i])
  design <- model.matrix(~0+factor(group_list[c(ctl,ind)]))
  colnames(design) <- levels(factor(c(groupA[i],groupB[i])))
  rownames(design) <- colnames(counts)[c(ctl,ind)]
  dge <- DGEList(counts=counts[,c(ctl,ind)])
  # calculate Norm Factor. default is TMM
  dge <- calcNormFactors(dge)
  # add 3 before logCPM, it can reduce low expression gene variation
  logCPM <- cpm(dge, log=TRUE, prior.count=3)
  v <- voom(dge,design, plot= TRUE, normalize="quantile")
  fit <- lmFit(v, design)
  com=unique(group_list[c(ind,ctl)])#250717 修改 使得上调为实验组上调
  comp=paste(com,collapse="-")
  cont.matrix <- makeContrasts(contrasts=c(comp),levels = design)
  fit2 <- contrasts.fit(fit,cont.matrix)
  fit2 <- eBayes(fit2)
  tmp <- topTable(fit2, coef=comp, n=Inf)
  DEG_limma_voom <- na.omit(tmp)
  #head(DEG_limma_voom)
  result=cbind(rownames(DEG_limma_voom),DEG_limma_voom)
  colnames(result)[1]="SYMBOL"
  limma_result_all[[i]]=result;
  ind2=which(abs(DEG_limma_voom$logFC)> log2FCcut & DEG_limma_voom$adj.P.Val< adjPcut)
  message(paste("found",length(ind2), "DEG in" ,paste(com,collapse=" vs "), "by limma."))
  loginfo <- paste0("log2FC cutoff: ",log2FCcut, "; adjust Pvalue cutoff: ",adjPcut,
                   "\n found ",length(ind2), " DEG in" ,paste(com,collapse=" vs "), " by limma.\n")
  
  if(length(ind2) == 0) 
    limma_result[[i]]=NA else 
      limma_result[[i]]=result[ind2,]
  names(limma_result_all)[i] = paste0(groupA[i],"vs",groupB[i])
  names(limma_result)[i] = paste0(groupA[i],"vs",groupB[i])
  
}

##edgeR #######
#### 第一步，构建edgeR的DGEList对象，并过滤
library(edgeR)
edgeR_result=list()
edgeR_result_all=list()

for(i in seq_along(groupB)){
  ctl = which(group_list == groupA[i])
  ind=which(group_list==groupB[i])
  design <- model.matrix(~0+factor(group_list[c(ctl,ind)]))
  colnames(design) <- levels(factor(c(groupA[i],groupB[i])))
  rownames(design) <- colnames(counts)[c(ctl,ind)]
  DEG <- DGEList(counts=counts[,c(ctl,ind)])
  keep <- filterByExpr(DEG,c(ctl,ind))
  table(keep)
  DEG <- DEG[keep ,keep.lib.sizes=FALSE]
  #DEG$samples$lib.size <- colSums(DEG$counts)
  DEG <- calcNormFactors(DEG)
  #### 第二步，差异表达分析
  DEG <- estimateGLMCommonDisp(DEG,design)
  DEG <- estimateGLMTrendedDisp(DEG, design)
  DEG <- estimateGLMTagwiseDisp(DEG, design)
  # 拟合线性模型
  fit <- glmFit(DEG, design)
  # 进行差异分析
  # 1,-1意味着前比后
  com=unique(group_list[c(ind,ctl)])
  comp=paste(com,collapse="-")
  cont.matrix <- makeContrasts(contrasts=c(comp),levels = design)
  lrt <- glmLRT(fit, contrast=cont.matrix) 
  #### 第三步，提取过滤差异分析结果
  edgeR_DEG <- topTags(lrt, n=nrow(DEG))
  edgeR_DEG <- as.data.frame(edgeR_DEG)
  result=cbind(rownames(edgeR_DEG),edgeR_DEG)
  colnames(result)[1]="SYMBOL"
  edgeR_result_all[[i]]=result
  ind2=which(abs(edgeR_DEG$logFC)> log2FCcut &edgeR_DEG$FDR < adjPcut);
  message(paste("found",length(ind2), "DEG in" ,paste(com,collapse=" vs "),"by edgeR."))
  loginfo <- paste(loginfo,
                   "found",length(ind2), "DEG in" ,paste(com,collapse=" vs "), "by edgeR.\n")
  
  if(length(ind2) == 0) 
    edgeR_result[[i]]=NA else
      edgeR_result[[i]]=result[ind2,]
  names(edgeR_result_all)[i] = paste0(groupA[i],"vs",groupB[i])
  names(edgeR_result)[i] = paste0(groupA[i],"vs",groupB[i])
  
}



##DESeq2 #########
suppressMessages(library(DESeq2))
#### 第一步，构建DESeq2的DESeq对象
counts_round=round(counts)
DESeq2_result=list()
DESeq2_result_all=list()
for(i in seq_along(groupB)){
  ctl = which(group_list == groupA[i])
  ind=which(group_list==groupB[i])
  com=unique(group_list[c(ctl,ind)])
  comp=paste(com,collapse="-")
  condition=factor(c(group_list[c(ctl,ind)]))
  colData <- data.frame(row.names=colnames(counts[,c(ctl,ind)]),group_list=condition)
  
  dds <- DESeqDataSetFromMatrix(countData = counts_round[,c(ctl,ind)],colData = colData,
                                design = ~group_list)
  #### 第二步，进行差异表达分析
  
  dds2 <- DESeq(dds)
  # 保存差异表达分析结果
  #### 第二步，提取差异分析结果
  res <- results(dds2)
  res <- res[order(res$padj),]
  DEG <- as.data.frame(res)
  result=cbind(rownames(DEG),DEG)
  colnames(result)[1]="SYMBOL";
  DESeq2_result_all[[i]]=result
  ind2=which(abs(DEG$log2FoldChange)>log2FCcut & DEG$padj < adjPcut & !is.na(DEG$padj)); #250918 old is ind
  message(paste("found",length(ind2), "DEG in" ,paste(com,collapse=" vs "),"by DESeq2."))
  loginfo <- paste(loginfo,
                   "found",length(ind2), "DEG in" ,paste(com,collapse=" vs "), "by DESeq2.\n")
  
  if(length(ind2) == 0) 
    DESeq2_result[[i]]=NA else 
      DESeq2_result[[i]]=result[ind2,]
  names(DESeq2_result_all)[i] = paste0(groupA[i],"vs",groupB[i])
  names(DESeq2_result)[i] = paste0(groupA[i],"vs",groupB[i])
}


##save ########
save(limma_result,edgeR_result,DESeq2_result,
     file = file.path(outputdir, paste0(output,"_DEG_troika.RData")))
save(limma_result_all,edgeR_result_all,DESeq2_result_all,
     file = file.path(outputdir, paste0(output,"_stat_troika.RData")))

loginfo <- paste(loginfo,
                 "   saved in",file.path(outputdir, paste0(output,"_stat_troika.RData")),"\n")
cat("************************************\n")
cat(loginfo)
cat("************************************\n")