#Version: 2.0
#created: Kefu liu(liukefu19@163.com)
#Maintainer: Kefu Liu
#Date: July 30, 2025
#Software: R

existpck <- rownames(installed.packages());
if(!"corrplot" %in% existpck)
  stop("no corrplot package, please install the R pakcage before run me.")
if(!"ggplot2" %in% existpck)
  stop("no ggplot2 package, please install the R pakcage before run me.")
if(!"ggpubr" %in% existpck)
  stop("no ggpubr package, please install the R pakcage before run me.")
if(!"ggrepel" %in% existpck)
  stop("no ggrepel package, please install the R pakcage before run me.")

suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(ggrepel))

##pvca###########
pvca_confounder_effect <- function(data,meta,myExpressionSet = NULL,batch.factors = NA,
                                   pct_threshold = 0.6) {
  if(is.null(myExpressionSet)) {
    metadata <- data.frame(labelDescription=colnames(meta),
                           row.names=colnames(meta))
    suppressPackageStartupMessages(require(Biobase))
    phenoData <- new("AnnotatedDataFrame",data=meta,varMetadata=metadata)
    myExpressionSet <- ExpressionSet(assayData=as.matrix(data),phenoData=phenoData,annotation="uniprot")#annotation随便写一个字符串就好，因为nsFilter函数要求Expression类必须有annotation
    myExpressionSet#
  }
  if(is.na(batch.factors)) batch.factors = colnames(meta);
  if(!all(batch.factors %in% colnames(meta))) 
    stop("wrong batch.factors setting")
  pvcaObj <- pvcaBatchAssess (myExpressionSet, batch.factors = batch.factors, 
                              threshold = pct_threshold)
  pvcaObj
}
pvca_plot <- function(pvcaObj, title = "", filename = NA,
                      width = 12,height = 6){
  if(!is.na(filename))  
    png(paste0(filename,".png"),width = 320*width,height = 320*height,res = 300)
  bp <- barplot(pvcaObj$dat, xlab = "Effects",
                ylab = "Weighted average proportion variance",
                ylim= c(0,1.1),col = c("blue"), las=2,  
                main=paste0("PVCA estimation (",title,")"))
  axis(1, at = bp, labels = pvcaObj$label, xlab = "Effects", cex.axis = 0.8)
  values = pvcaObj$dat
  new_values = round(values , 3)
  text(bp,pvcaObj$dat,labels = new_values, pos=3, cex = 0.8)
  if(!is.na(filename))  
    graphics.off()
  if(is.na(filename)) pvcaObj
}

pvca_evaluate <- function(pvcadata,exp_design,pct_threshold = 0.6, output = NA) {
  if(!requireNamespace("lme4", quietly = TRUE)) 
    stop("no lme4 package, please install the R pakcage before run me.")
  suppressPackageStartupMessages(library(lme4))
  ## Load data 
  dataRowN <- nrow(pvcadata)
  dataColN <- ncol(pvcadata)
  expDesignRowN <- nrow(exp_design)
  expDesignColN <- ncol(exp_design)
  myColNames <- names(exp_design)
  
  # Center the data (center rows)
  pvcadata_Centered <- matrix(data = 0, nrow = dataRowN, ncol = dataColN)
  pvcadata_Centered_transposed = apply(pvcadata, 1, scale, center = TRUE, scale = FALSE)
  pvcadata_Centered = t(pvcadata_Centered_transposed)
  # Compute correlation matrix: gene variation
  theDataCor <- cor(pvcadata_Centered)
  # Obtain eigenvalues 
  eigenData <- eigen(theDataCor)
  eigenValues = eigenData$values
  ev_n <- length(eigenValues)
  eigenVectorsMatrix = eigenData$vectors
  eigenValuesSum = sum(eigenValues)
  percents_PCs = eigenValues /eigenValuesSum
  #Merge experimental file and eigenvectors for n components 
  pc_n <- max(3,min(which(cumsum(percents_PCs) > pct_threshold)))
  {#my_counter_2 = 0
  #my_sum_2 = 1
  #for (i in ev_n:1){
  #  my_sum_2  = my_sum_2 - percents_PCs[i]
  #  if ((my_sum_2) <= pct_threshold ){
  #    my_counter_2 = my_counter_2 + 1
  #  }
  #}
  #if (my_counter_2 < 3){
  #  pc_n  = 3
  #}else {
  #  pc_n = my_counter_2 
  #}
  }
  # pc_n is the number of principal components to model
  pc_data_matrix <- as.numeric(eigenVectorsMatrix[,1:pc_n])
  {#pc_data_matrix <- matrix(data = 0, nrow = (expDesignRowN*pc_n), ncol = 1)
  #mycounter = 0
  #for (i in 1:pc_n){
  #  for (j in 1:expDesignRowN){
  #    mycounter <- mycounter + 1
  #    pc_data_matrix[mycounter,1] = eigenVectorsMatrix[j,i]
  #  }
  #}
  }
  AAA <- data.frame(exp_design[rep(1:expDesignRowN,pc_n),])
  colnames(AAA) <- colnames(exp_design)
  Data <- cbind(AAA,pc_data_matrix)
  #Edit these variables according to your factors
  variables <- c(colnames(exp_design))
  for (i in 1:length(variables)) {
    #Data$variables[i] <- as.factor(Data$variables[i]) 
    if(is.character(Data[,variables[i]]))
      Data[,variables[i]] <- as.factor(Data[,variables[i]])#250802
  }
  # Mixed linear model 
  op <- options(warn = (-1))
  # effects_n = expDesignColN + choose(expDesignColN, 2) + 1
  effects_n = expDesignColN + 1
  randomEffectsMatrix <- matrix(data = 0, nrow = pc_n, ncol = effects_n)
  model.func <- c()
  index <- 1
  for (i in 1:length(variables)) {
    mod = paste("(1|", variables[i], ")", sep = "")
    model.func[index] = mod
    index = index + 1
  }
  function.mods <- paste(model.func, collapse = " + ")
  for (i in 1:pc_n) {
    y = (((i - 1) * expDesignRowN) + 1)
    funct <- paste("pc_data_matrix", function.mods, sep = " ~ ")
    Rm1ML <- lmer(funct, Data[y:(((i - 1) * expDesignRowN) + 
                                   expDesignRowN), ], REML = TRUE, control=lmerControl(check.nobs.vs.nlev = "ignore",check.nobs.vs.rankZ = "ignore",check.nobs.vs.nRE="ignore"),verbose = FALSE, 
                  na.action = na.omit)
    randomEffects <- Rm1ML
    randomEffectsMatrix[i, ] <- c(unlist(VarCorr(Rm1ML)), 
                                  resid = sigma(Rm1ML)^2)
  }
  effectsNames <- c(names(getME(Rm1ML, "cnms")), "resid")
  colnames(randomEffectsMatrix) <- effectsNames
  rownames(randomEffectsMatrix) <- colnames(eigenVectorsMatrix)[1:pc_n]
  ## Standardize Variance
  randomEffectsMatrixStdze <- (randomEffectsMatrix)/apply(randomEffectsMatrix,1,sum)
  # Compute Weighted Proportions 
  randomEffectsMatrixWtProp <- (randomEffectsMatrixStdze)*percents_PCs[1:pc_n]
  # Compute Weighted Ave Proportions
  randomEffectsSums <- apply(randomEffectsMatrixWtProp,2,sum)
  randomEffectsMatrixWtAveProp <- randomEffectsSums/sum(randomEffectsSums)
  {#randomEffectsMatrixStdze <- matrix(data = 0, nrow = pc_n, ncol = effects_n)
  #for (i in 1:pc_n){
  #  mySum = sum(randomEffectsMatrix[i,])
  #  for (j in 1:effects_n){
  #    randomEffectsMatrixStdze[i,j] = randomEffectsMatrix[i,j]/mySum	
  #  }
  #}
  
  #randomEffectsMatrixWtProp <- matrix(data = 0, nrow = pc_n, ncol = effects_n)
  #for (i in 1:pc_n){
  #  weight = eigenValues[i]/eigenValuesSum
  #  for (j in 1:effects_n){
  #    randomEffectsMatrixWtProp[i,j] = randomEffectsMatrixStdze[i,j]*weight
  #  }
  #}
  
  #randomEffectsSums <- matrix(data = 0, nrow = 1, ncol = effects_n)
  #randomEffectsSums <-colSums(randomEffectsMatrixWtProp)
  #totalSum = sum(randomEffectsSums)
  #randomEffectsMatrixWtAveProp <- matrix(data = 0, nrow = 1, ncol = effects_n)
  #for (j in 1:effects_n){
  #  randomEffectsMatrixWtAveProp[j] = randomEffectsSums[j]/totalSum 	
  #}
  }
  #plot
  pvcaObj <- data.frame(var = effectsNames,PVCA_Prop =as.vector(randomEffectsMatrixWtAveProp))
  pvcaObj$var <- factor(pvcaObj$var, levels = effectsNames)
  pic <- ggplot(pvcaObj,aes(x = var, y = PVCA_Prop)) +
    geom_col() + 
    geom_text(aes(label = round(PVCA_Prop,3)), vjust = -0.5,
              hjust = 0.5, color = "black") +
    ylim(c(0, 1)) +
    labs(x = "Covariates",y = "Proportion of PVCA") +
    theme_bw()
  if(!is.na(output)) pic <- pic + labs(title = output)
  list(pvca_res = list(randomEffectsMatrix, percents_PCs = percents_PCs[1:pc_n]),pic = pic)
}


###pca#######
#v2.0
pca_dat <- function(exprM, scale = TRUE,keep_dim = 20,...){
  sum <- apply(exprM,1,sum,na.rm=TRUE);
  #row have only 1 value will remove 
  Not1 <- apply(exprM,1,function(x) length(unique(x))!=1) #240515 fix error when gene have the same value
  exprM <- exprM[!is.na(sum)&Not1,];
  data <- t(as.matrix(exprM));
  # do PCA 
  data.pca <- try(prcomp(data,scale = scale,...),silent = TRUE)
  if(inherits(data.pca,"try-error")) {
    message("have NA value, setting to 0 to run pca.")
    data[is.na(data)] <- 0
    #251126 fix error when gene have the same value
    Not1 <- apply(data,2,function(x) length(unique(x))!=1) 
    data <- data[,Not1];
    data.pca <- try(prcomp(data,scale = scale,...),silent = TRUE)
  }
  if (!inherits(data.pca,"try-error")) {
    if(ncol(data.pca$rotation) < keep_dim) keep_dim = ncol(data.pca$rotation)
    data.pca$rotation <- data.pca$rotation[,1:keep_dim]
    data.pca$x <- data.pca$x[,1:keep_dim]
  } else {data.pca = NA; message("error in run PCA analysis.")}
  data.pca
}
#exprM can be a prcomp data or expression matrix
pcaplot<- function(exprM, group, name, scale = TRUE,label=TRUE, addzeroline = FALSE,
                   point_size=4,point_shape=19,...){
  if(missing(group))  group <- gsub("[0-9]","",colnames(exprM));
  if(missing(name)) name = "PCA plot"
  if(class(exprM)[1] == "prcomp") {data.pca <- exprM;
  exprM <- t(data.pca$x)
  } else if(class(exprM)[1] %in% "data.frame" | class(exprM)[1] %in% "matrix") {
    if(ncol(exprM) != length(group)) stop("group length is not equal with exprM sample num.")
    sum <- apply(exprM,1,sum);
    Not1 <- apply(exprM,1,function(x) length(unique(x))!=1) #240515 fix error when gene have the same value
    exprM <- exprM[!is.na(sum)&Not1,];
    data <- t(as.matrix(exprM));
    # do PCA 
    data.pca <- prcomp(data,scale = scale)
  }
  # fetch the proportion of PC1 and PC2
  pc <- as.data.frame(data.pca$x);
  importance <- summary(data.pca)$importance;
  p1 <- as.numeric(sprintf("%.3f",importance[2,1]))*100;
  p2 <- as.numeric(sprintf("%.3f",importance[2,2]))*100;
  pc$group = group;  
  if(length(label) == ncol(exprM)) {
    pc$names = label; label = TRUE
  } else pc$names = rownames(pc);
  if(!is.logical(label) && length(label) != ncol(exprM)) {
    label = FALSE; warning("label length is not equal with exprM sample num, foced setting label to FALSE")
  } 
  # draw PCA plot figure
  xlab=paste0("PC1(",p1,"%)");
  ylab=paste0("PC2(",p2,"%)");
  title = paste0("PCA plot (",name,")");
  p <- ggplot(pc,aes(PC1,PC2)) + 
    theme_bw() +
    theme(panel.grid.major =element_blank(), 
          panel.grid.minor = element_blank())+
    geom_point(aes(color=group),size=point_size,shape=point_shape,...) + 
    labs(x=xlab,y=ylab,title=title) + 
    theme(plot.title=element_text(hjust=0.5),legend.title = element_blank())
  if(isTRUE(label)) 
    p <- p + geom_text(aes(label=names,vjust=2,hjust=0.5),size=3)
  if(addzeroline) {
    p <- p + geom_hline(yintercept=0,linetype=4,color="grey") + 
      geom_vline(xintercept=0,linetype=4,color="grey")}
  p
}
pcasave <- function(pca, filename, width=8, height=6){
  png(paste0(filename,".png"),width=width*80*4, height=height*80*4,res=300)
  print(pca)
  graphics.off()
}

#v2.0
pca_plot_var <- function(pcadata, meta, title = "",
                         point_size=4,point_shape=19, na.value = "grey50",
                         gradient_color = rainbow(5),...) {
  if(class(pcadata)!="prcomp") stop("pcadata must be a prcomp object.")
  if(all(rownames(pcadata$x) != rownames(meta)))
    stop("pcadata and meta must have the same sample and sequence.")
  for(i in 1:ncol(meta)) {
    if(length(unique(meta[,i])) > 1) {
      if(is.numeric(meta[,i])) {
        p1 <- pcaplot(exprM = pcadata,  group = meta[,i],
                      name = title,label = FALSE,
                      point_size=point_size,point_shape=point_shape,...) +
          theme(legend.title = element_text()) +labs(color = colnames(meta)[i]) +
          scale_color_gradientn(colors = gradient_color,na.value = na.value)
      } else {
        p1 <- pcaplot(exprM = pcadata, group= meta[,i],
                      name = title,label = FALSE,
                      point_size=point_size,point_shape=point_shape,...) +
          theme(legend.title = element_text()) +labs(color = colnames(meta)[i])
      }
      print(p1)}
  }
}

#Copy from PCAForQTL
runElbow<-function(X=NULL,
                   prcompResult=NULL){
  
  #Obtain prcompResult.
  if(is.null(prcompResult)){
    if(is.null(X)){
      stop("Please input X or prcompResult.")
    }else{
      cat("Running PCA...\n")
      prcompResult<-prcomp(X,center=TRUE,scale.=TRUE) #This may take a moment.
    }
  }else{
    if(class(prcompResult)!="prcomp"){
      stop("prcompResult must be a prcomp object returned by the function prcomp().")
    }
  }
  
  importanceTable<-summary(prcompResult)$importance
  x<-1:ncol(importanceTable) #PC indices.
  y<-importanceTable[2,] #PVEs.
  
  #Given x  and y, calculate the distance between each point and the diagonal line (the line connecting the first and last points).
  x1<-x[1] #First point.
  y1<-y[1] #First point.
  
  x2<-x[length(x)] #Last point.
  y2<-y[length(y)] #Last point.
  
  x0<-x
  y0<-y
  
  distancesDenominator<-sqrt((x2-x1)^2+(y2-y1)^2)
  distancesNumerator<-abs((x2-x1)*(y1-y0)-(x1-x0)*(y2-y1))
  distances<-distancesNumerator/distancesDenominator
  # plot(distances)
  
  numOfPCsChosen<-which.max(distances) #12.
  names(numOfPCsChosen)<-NULL
  return(numOfPCsChosen)
}
runBE<-function(X,B=20,alpha=0.05,
                mc.cores=min(B,parallel::detectCores()-1),
                verbose=FALSE){
  
  if(alpha<0 || alpha>1){
    stop("alpha must be between 0 and 1.")
  }
  
  n<-nrow(X) #Number of observations.
  p<-ncol(X) #Number of features.
  d<-min(n,p) #This is the total number of PCs.
  
  if(verbose) cat("Running PCA on permuted data...\n")
  
  # testStatsPerm<-matrix(data=NA,nrow=d,ncol=B) #PC by permutation.
  # for(b in 1:B){
  #   # b<-3
  #   if(verbose) cat("b=",b," out of ",B," permutations...\n",sep="")
  #
  #   #Permute each column of X. That is, permute the observations in each feature.
  #   XPermuted<-matrix(data=NA,nrow=n,ncol=p)
  #   for(j in 1:p){
  #     # j<-7
  #     XPermuted[,j]<-sample(x=X[,j],size=n,replace=FALSE)
  #   }
  #
  #   prcompResultPerm<-prcomp(x=XPermuted,center=TRUE,scale.=TRUE) #Key step.
  #   importanceTablePerm<-summary(prcompResultPerm)$importance
  #   testStatsPerm[,b]<-importanceTablePerm[2,] #The second row is PVE.
  # }
  
  results<-parallel::mclapply(1:B,FUN=function(b){
    # b<-3
    if(verbose) cat("b=",b," out of ",B," permutations...\n",sep="")
    
    #Permute each column of X. That is, permute the observations in each feature.
    XPermuted<-matrix(data=NA,nrow=n,ncol=p)
    for(j in 1:p){
      # j<-7
      XPermuted[,j]<-sample(x=X[,j],size=n,replace=FALSE)
    }
    
    prcompResultPerm<-prcomp(x=XPermuted,center=TRUE,scale.=TRUE) #Key step.
    importanceTablePerm<-summary(prcompResultPerm)$importance
    PVEsPerm<-importanceTablePerm[2,]
    return(PVEsPerm)
  },mc.cores=mc.cores) #results is a list of vectors.
  temp<-unlist(results)
  testStatsPerm<-matrix(data=temp,nrow=d,byrow=FALSE) #PC by permutation.
  
  if(verbose) cat("Running PCA on the unpermuted data...\n")
  prcompResult<-prcomp(x=X,center=TRUE,scale.=TRUE) #Key step.
  importanceTable<-summary(prcompResult)$importance
  PVEs<-importanceTable[2,]
  # Compare PVEs to testStatsPerm.
  # temp<-(testStatsPerm>=PVEs) #temp is calculated as desired.
  pValues<-(rowSums(testStatsPerm>=PVEs)+1)/(B+1) #The p-value for the jth PC is calculated as, roughly speaking, the proportion of permutations where the PVE of the jth PC is greater than or equal to PVE_j.
  
  for(j in 2:d){ #Enforce monotone increase of the p-values.
    if(pValues[j]<pValues[j-1]){
      pValues[j]<-pValues[j-1]
    }
  }
  
  numOfPCsChosen<-sum(pValues<=alpha)
  toReturn<-list(pValues=pValues,alpha=alpha,numOfPCsChosen=numOfPCsChosen)
  return(toReturn)
}
makeScreePlot<-function(prcompResult,labels,values,
                        titleText=NULL,subtitleText=NULL,
                        maxNumOfPCsToPlot=max(100,values),colors=NULL){
  
  if(class(prcompResult)!="prcomp") stop("prcompResult must be a prcomp object returned by the function prcomp().")
  if(!is.character(labels)) stop("labels must be a string or vector of strings.")
  if(!is.numeric(values)) stop("values must be an integer or vector of integers.")
  if(length(labels)!=length(values)) stop("labels and values must have the same length.")
  
  importanceTable<-summary(prcompResult)$importance
  dataPlot<-as.data.frame(t(importanceTable))
  dataPlot$PCIndex<-1:nrow(dataPlot)
  colnames(dataPlot)<-c("sd","PVE","cumulativePVE","PCIndex")
  # dataPlot<-dataPlot%>%filter(PCIndex<=maxNumOfPCsToPlot)
  dataPlot<-dataPlot[which(dataPlot$PCIndex<=maxNumOfPCsToPlot),] #Avoid using dplyr.
  
  dataPlotVline<-data.frame(matrix(nrow=length(labels),ncol=2))
  colnames(dataPlotVline)<-c("Category","Value")
  dataPlotVline$Category<-labels
  dataPlotVline$Category<-factor(dataPlotVline$Category,levels=dataPlotVline$Category)
  dataPlotVline$Value<-values
  if(is.null(colors)){
    colors<-RColorBrewer::brewer.pal(n=8,name="Accent")[1:length(labels)]
  }
  
  p<-ggplot(data=dataPlot,aes(x=PCIndex,y=PVE))+
    geom_point()+
    geom_vline(data=dataPlotVline, #Key.
               aes(xintercept=Value,color=Category),
               linewidth=0.75)+ #Default size is 0.5.
    scale_x_continuous(breaks=seq(0,nrow(dataPlot),by=10))+
    scale_y_continuous(labels=scales::number_format(accuracy=0.01))+
    scale_color_manual(values=colors)+
    labs(x="PC index")+
    theme(
      text=element_text(size=15),
      
      plot.title=element_text(hjust=0.5), #Plot title.
      plot.subtitle=element_text(hjust=0.5), #Plot subtitle.
      
      panel.background=element_blank(), #Panel background.
      panel.grid.major.y=element_blank(), #Horizontal grid lines.
      panel.grid.major.x=element_blank(), #Vertical grid lines.
      
      axis.line=element_line(), #Axis lines (default size is 0.5),
      
      legend.title=element_blank(), #Legend title.
      legend.key=element_blank(), #Legend key background.
      legend.position="bottom", #Legend position.
    )
  
  if(!is.null(titleText)) p<-p+labs(title=titleText)
  if(!is.null(subtitleText)) p<-p+labs(subtitle=subtitleText)
  
  # p
  return(p)
}
filterKnownCovariates<-function(knownCovariates,inferredCovariates,unadjustedR2_cutoff=0.9,
                                verbose=FALSE){
  
  if(nrow(knownCovariates)!=nrow(inferredCovariates)){ #Check whether the sample sizes are equal in case the rownames are empty.
    stop("knownCovariates and inferredCovariates must have the same number of rows (i.e., observations).")
  }
  if(!identical(rownames(knownCovariates),rownames(inferredCovariates))){ #identical(NULL,NULL) returns TRUE.
    stop("The rownames of knownCovariates and inferredCovariates must match.")
  }
  if(unadjustedR2_cutoff<0 || unadjustedR2_cutoff>1){
    stop("unadjustedR2_cutoff must be between 0 and 1.")
  }
  
  R2s<-groupPredict(dataResponse=knownCovariates,dataPredictors=inferredCovariates,R2Type="unadjusted") #Use unadjusted because we don't want to penalize for model complexity here.
  indicesOfKnownCovariatesToKeep<-which(R2s<unadjustedR2_cutoff) #This may be integer(0).
  if(verbose){
    cat(ncol(knownCovariates)-length(indicesOfKnownCovariatesToKeep)," out of the ",ncol(knownCovariates)," known covariates has/have been filtered out.\n")
  }
  
  toReturn<-knownCovariates[,indicesOfKnownCovariatesToKeep] #This syntax works even when all of the known covariates are filtered out.
  return(toReturn)
}
groupPredict<-function(dataResponse,dataPredictors,R2Type=c("unadjusted","adjusted")){
  # set.seed(1)
  # n<-100
  # dataResponse<-matrix(data=rnorm(n=n*3),nrow=n) #100*3.
  # dataPredictors<-matrix(data=rnorm(n=n*10),nrow=n) #100*10.
  # R2Type<-"unadjusted"
  
  R2s<-rep(NA,ncol(dataResponse))
  for(j in 1:ncol(dataResponse)){
    # j<-1
    mod<-lm(dataResponse[,j]~as.matrix(dataPredictors))
    # summary(mod)
    if(R2Type=="unadjusted"){
      R2s[j]<-summary(mod)$r.squared
    }else{
      R2s[j]<-summary(mod)$adj.r.squared
    }
  }
  
  return(R2s)
}

###BIC##########
BIC_evaluate <- function(quant,covar,num_cores = 30, output = NULL) {
  if (!requireNamespace("doParallel", quietly = TRUE)) {
    stop("No doParallel package, cannot run BIC. Please install: install.packages('doParallel')", 
         call. = FALSE)
  }
  if (!requireNamespace("foreach", quietly = TRUE)) {
    stop("No foreach package, cannot run BIC. Please install: install.packages('foreach')", 
         call. = FALSE)
  }
  if(is.null(covar)){
    stop("covar is NULL, cannot run BIC_evaluate.")
  }
  if(ncol(quant) != nrow(covar)) stop("quant and covar have not the same number.")
  for (i in 1:ncol(covar)) {
    if(!class(covar[,i]) %in% c("numeric","integer"))
      covar[,i] <- as.numeric(as.factor(covar[,i]))
  }
  suppressPackageStartupMessages(library(doParallel))
  suppressPackageStartupMessages(library(foreach))
  d <- cbind(t(quant), covar)
  nGene <- nrow(quant)
  nSample <- nrow(covar)
  # Set up parallel computing
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  #pcs <- paste(paste0("d$PC", 1:35), collapse = " + ")
  results <- foreach(i = 1:nGene, .packages = c("stats"), .combine = c) %dopar% {
    fmla <- as.formula(paste(colnames(d)[i],"~ ", paste(colnames(covar), collapse= "+")))
    resBIC <- try(step(lm(fmla, data = d), k = log(nSample), direction = "both"),silent = TRUE)
    if(!inherits(resBIC,"try-error"))
      return(list(names(attr(resBIC$terms, "dataClasses"))[-1])) else return(NA)
    #formula <- as.formula(paste("d[, i] ~ d$PC1 + d$PC2 + d$Group + d$Sex"))
    #resBIC <- step(lm(formula, data = d), k = log(nSample), direction = "both")
    #log_file <- paste("BIC/log.BIC.", i, sep = "")
    #capture.output(resBIC, file = log_file)
    # Return BIC variables
  }
  if(all(is.na(results))) stop("cannot run BIC.")
  names(results) <- rownames(quant)
  stopCluster(cl)
  # Combine results
  covBIC <- unlist(results)
  # Prepare for plotting
  #save.image(file="bicresult.RData")
  bicplotdata <- data.frame((table(covBIC)[colnames(covar)] / nGene))
  #level is arrange by origin sequence
  bicplotdata$covBIC <- factor(colnames(covar), levels = colnames(covar))
  bicplotdata$Freq[is.na(bicplotdata$Freq)] <- 0
  #bicplotdata$covBIC <- factor(bicplotdata$covBIC, levels = bicplotdata$covBIC)
  pic <- ggplot(bicplotdata,aes(x = covBIC, y = Freq)) +
    geom_col() + 
    geom_text(aes(label = round(Freq,3)), vjust = -0.5,
              hjust = 0.5, color = "black") +
    ylim(c(0, 1)) +
    labs(x = "Covariates",y = "Proportion of genes decreased BIC") +
    theme_bw()
  if(!is.null(output)) {
    pic <- pic + labs(title = output)
  }
  list(BICres = results, pic = pic)
}
##MDS######
MDS_dat <- function(quant) {
  sum <- apply(quant,1,sum,na.rm=TRUE);
  Not1 <- apply(quant,1,function(x) length(unique(x))!=1) 
  quant <- quant[!is.na(sum)&Not1,];
  # do MDS 
  data.mds <- try(cmdscale(dist(t(quant)), eig = T),silent = TRUE)
  if(inherits(data.mds,"try-error")) {
    message("have NA value, setting to 0 to run MDS.")
    quant[is.na(quant)] <- 0
    data.mds <- try(cmdscale(dist(t(quant)), eig = T),silent = TRUE)
  }
  if (inherits(data.mds,"try-error")){
    data.mds = NA; message("error in run MDS analysis.")}
  data.mds
}
MDS_plot <- function(data.mds, covar=NULL, group=NA, title = "",
                     filename = NA, width=8, height=6){
  drawcolor = FALSE;
  if(!inherits(data.mds,"list")) stop("data.mds must be a list object.")
  if(!is.null(covar) & !is.na(group)) {
    if(!group %in% colnames(covar)) 
      stop(paste0("no ", group ," column in covar."))
    if(all(rownames(data.mds$points) != rownames(covar)))
      stop("data.mds and covar must have the same sample and sequence.")
    drawcolor = TRUE;
  } else message("missing covar or group, no group info will draw on fig.")
  mds <- data.mds
  pc1 = mds$eig[1]^2 / sum(mds$eig^2)   
  pc2 = mds$eig[2]^2 / sum(mds$eig^2)
  if(!is.na(filename)) png(paste0(filename,".png"),width=width*80*4, height=height*80*4,res=300)
  if(drawcolor) {
    covar[,group] = as.factor(covar[,group])
    plot(mds$points, col=as.numeric(as.factor(covar[,group])), pch=20, 
         main=paste0(title," Multidimensional Scaling Plot"), 
         xlab = paste("eig1 (", signif(100*pc1,3), "%)", sep=""),
         ylab = paste("eig2 (", signif(100*pc2,3),"%)",sep=""))
    legend("topright", (levels(covar[,group])), col=unique(as.numeric(covar[,group])), pch=16, cex=0.5)
  } else 
    plot(mds$points, pch=20, 
         main=paste0(title," Multidimensional Scaling Plot"),
         xlab = paste("eig1 (", signif(100*pc1,3), "%)", sep=""),
         ylab = paste("eig2 (", signif(100*pc2,3),"%)",sep=""))
  if(!is.na(filename)) graphics.off()
}


##density plot########
# sample distribution 
density_plot <- function(quant, covar=NULL, group=NA, title = "",
                         filename = NA, width=8, height=6) {
  drawcolor = FALSE;
  if(!is.null(covar) & !is.na(group)) {
    if(!group %in% colnames(covar)) 
      stop(paste0("no ", group ," column in covarfile."))
    if(nrow(covar) != ncol(quant)) stop("quant and covar is not match.")
    drawcolor = TRUE;
    group <- (as.factor((covar[,group])))
  } else message("missing covar or group, no group info will draw on fig.")
  
  if(!is.na(filename)) png(paste0(filename,".png"),width=width*80*4, height=height*80*4,res=300)
  ymax = 0
  for(i in 1:ncol(quant)) {
    ymax2 = max(density(quant[,i], na.rm=T)$y);
    if(ymax2 > ymax) k = i;
    ymax = max(ymax,ymax2)
  }
  i = k;
  if(drawcolor) 
    plot(density(quant[,i], na.rm=T), col = as.numeric(group)[i], main=paste0(title,": Hist of expression"), xlab = "expression") else
      plot(density(quant[,i], na.rm=T), main=paste0(title,": Hist of expression"), xlab = "exp")
  for(i in (1:ncol(quant))[-k]) {
    if(drawcolor) lines(density((quant[,i]), na.rm=T), col = as.numeric(group)[i]) else
      lines(density((quant[,i]), na.rm=T))
  } 
  if(drawcolor) legend("topright", levels(group), text.col =unique(as.numeric(group),cex=1))
  if(!is.na(filename)) graphics.off()
}

#data evaluate pick 4 gene to draw the plot.
densityplot <- function(quant,output, savefile = TRUE,
                        width = 12, height = 10,verbose = 0){
  if(!"opt" %in% ls()) {
    opt <- list();
    opt$grpname = NA;
    if("group" %in% ls()) {rm(group);message("group obj is force to removed")}
  }
  if(!dir.exists("plot")) dir.create("plot")
  #251126 fix some error
  noNAnum <- rowSums(!is.na(quant))
  #remove all NA row
  allNA <- noNAnum == 0
  quant <- quant[!allNA,]
  if(nrow(quant) == 0) stop("something wrong. All NA value.")
  #remove Not 1 row
  Not1  = apply(quant,1,function(x) length(unique(na.omit(x))) != 1)
  quant <- quant[Not1,]
  if(nrow(quant) == 0) stop("something wrong. value is unqique each row.")
  #ncol > 30 but each row have many NA (noNA num > 30) row < 5) will adjust to 29 column.
  noNAnum <- rowSums(!is.na(quant))
  if(ncol(quant) > 30) {
    if(sum(noNAnum >= 30) < 5) {
      quant <- quant[noNAnum < 30,]
      quant <- as.data.frame(
        t(apply(quant,1, function(x) {
        noNAvalue = head(na.omit(x),29);
        padded_value <- c(noNAvalue, rep(NA, 29 - length(noNAvalue)));
        return(padded_value)})))
    }
     }
  meanD = apply(quant,1,mean,na.rm=TRUE)
  
  if(ncol(quant) < 30) {
    if(verbose > 0) cat("sample less than 30, will merge by mean value in density plot.\n")
    pickn = ceiling(30/ncol(quant));
    sortmu <- sort(meanD)
    diff = (sortmu[-1]-sortmu[-nrow(quant)])
    ratio = abs(1000*diff/(sortmu[-1])) < 1
    continuepos = rle(ratio)
    ##ratio<1 and length > pickn
    continuepos2 = which(continuepos$lengths>pickn & continuepos$values == TRUE)
    if(length(continuepos2) < 4 ) stop("No enough data can draw.")
    posmerge = match(names(continuepos2)[1], names(sortmu))
    xname = paste0(names(sortmu)[posmerge - 1:pickn],collapse = "/")
    g = data.frame(ID1 = as.vector(as.matrix(quant[names(sortmu)[posmerge - 1:pickn],])))
    
    posmerge = match(names(continuepos2)[length(continuepos2)], names(sortmu))
    xname = c(xname,paste0(names(sortmu)[posmerge - 1:pickn],collapse = "/"))
    g = data.frame(g,ID2 = as.vector(as.matrix(quant[names(sortmu)[posmerge - 1:pickn],])))
    
    posmerge = match(names(continuepos2)[trunc(length(continuepos2)/3)], names(sortmu))
    xname = c(xname,paste0(names(sortmu)[posmerge - 1:pickn],collapse = "/"))
    g = data.frame(g,ID3 = as.vector(as.matrix(quant[names(sortmu)[posmerge - 1:pickn],])))
    
    posmerge = match(names(continuepos2)[trunc(2*length(continuepos2)/3)], names(sortmu))
    xname = c(xname,paste0(names(sortmu)[posmerge - 1:pickn],collapse = "/"))
    g = data.frame(g,ID4 = as.vector(as.matrix(quant[names(sortmu)[posmerge - 1:pickn],])))
    
    p1 = ggplot(g,aes(x = ID1))+ theme_bw()+
      xlab(paste0("maxmean: ",xname[1],"(p=",signif(shapiro.test(g[,1])$p.value,digits=2),")")) 
    p2 = ggplot(g,aes(x = ID2))+ theme_bw() +
      xlab(paste0("minmean: ",xname[2],"(p=",signif(shapiro.test(g[,2])$p.value,digits=2),")")) 
    p3 = ggplot(g,aes(x = ID3))+theme_bw()+
      xlab(paste0("33%mean: ",xname[3],"(p=",signif(shapiro.test(g[,3])$p.value,digits=2),")")) 
    p4 = ggplot(g,aes(x = ID4))+ theme_bw()+
      xlab(paste0("66%mean: ",xname[4],"(p=",signif(shapiro.test(g[,4])$p.value,digits=2),")")) 
    p1 = p1 + geom_density(fill = "gray",color = "black")
    p2 = p2 + geom_density(fill = "gray",color = "black")
    p3 = p3 + geom_density(fill = "gray",color = "black")
    p4 = p4 + geom_density(fill = "gray",color = "black")
    
  } else {
    quant <- quant[noNAnum >= 30,]
    meanD <- meanD[noNAnum >= 30]
    varD = apply(quant,1,var,na.rm = TRUE)
    #pick maxvar minvar maxmean minmean 
    pos = order(meanD,decreasing = TRUE)[1]
    pos = c(pos,order(meanD)[1])
    pos = c(pos,order(varD,decreasing = TRUE)[1])
    if(any(pos[3] == pos[-3]))
      pos[3] = order(varD,decreasing = TRUE)[2]
    pos = c(pos,order(varD)[1])
    if(any(pos[4] == pos[-4]))
      pos[4] = order(varD)[2]
    g = data.frame(t(quant[pos,]))
    colnames(g) <- c("ID1","ID2","ID3","ID4")
    xname = rownames(quant)[pos]
    if(!is.na(opt$grpname)) {
      group = covar[,opt$grpname]
      g = data.frame(g,group)}
    
    p1 = ggplot(g,aes(x = ID1))+ theme_bw()+
      xlab(paste0("maxmean: ",xname[1],"(p=",signif(shapiro.test(g[,1])$p.value,digits=2),")")) 
    p2 = ggplot(g,aes(x = ID2))+ theme_bw() +
      xlab(paste0("minmean: ",xname[2],"(p=",signif(shapiro.test(g[,2])$p.value,digits=2),")")) 
    p3 = ggplot(g,aes(x = ID3))+theme_bw()+
      xlab(paste0("maxvar: ",xname[3],"(p=",signif(shapiro.test(g[,3])$p.value,digits=2),")")) 
    p4 = ggplot(g,aes(x = ID4))+theme_bw()+
      xlab(paste0("minvar: ",xname[4],"(p=",signif(shapiro.test(g[,4])$p.value,digits=2),")")) 
    if("group" %in% ls()) {
      p1 = p1 + geom_density(aes(fill = group),color = "black",alpha = 0.4)
      p2 = p2 + geom_density(aes(fill = group),color = "black",alpha = 0.4)
      p3 = p3 + geom_density(aes(fill = group),color = "black",alpha = 0.4)
      p4 = p4 + geom_density(aes(fill = group),color = "black",alpha = 0.4)
    } else{
      p1 = p1 + geom_density(fill = "gray",color = "black",na.rm = TRUE)
      p2 = p2 + geom_density(fill = "gray",color = "black",na.rm = TRUE)
      p3 = p3 + geom_density(fill = "gray",color = "black",na.rm = TRUE)
      p4 = p4 + geom_density(fill = "gray",color = "black",na.rm = TRUE)
    }
  }
  toptitle = paste0("density_plot_",output)
  pall = ggarrange(p1, p2, p3, p4, nrow = 2, ncol = 2, labels = c('A', 'B', 'C', 'D'), 
                   font.label = list(color = 'black'))
  if(savefile) pdf(file = paste0("plot/density_plot_",output,".pdf"),width = width, height = height)
  print(pall)    
  if(savefile) dev.off()
}

###Outlier_detect##########

Outlier_Detect <- function(data, iteration = NA, intensity = "intensity",maxNA=0.5,
                           distmethod = "manhattan", A.IAC = FALSE, sdout = 2,
                           plot = TRUE, filename = NULL,
                           text.labels = NA, abline.col = "red", abline.lwd = 1){
  distmethod <- match.arg(distmethod,
                          c("manhattan","euclidean", "canberra","correlation","bicor"));
  if(!is.list(data)) stop("data should be a list and contain inf and intensity.")
  if( !"inf" %in% names(data) ) stop("data have no inf.")
  if(length(intensity) != 1) stop("intenstiy should be a single vector.")
  if( !intensity %in% names(data) ) stop(paste0("data have no ",intensity,"."))
  inf <- data.frame(data[["inf"]]); 
  quant <- data[[intensity]]
  if (!is.na(miss.value)) quant[quant == miss.value] = NA;
  #remove rows which have 0,1,2 value
  NArow <- rowSums(!is.na(quant))
  quant <- quant[NArow > 2,]
  inf <- data.frame(inf[NArow > 2, ]); 
  #if sample have only 4 values and missing rate > 0.9 will remove and warning.
  NAcol <- colSums(!is.na(quant))
  NAcol_ratio <- 1-NAcol/nrow(quant)
  #it means some sample may have some problem
  if(max(NAcol_ratio) > 0.9 & min(NAcol) < 4) {
    pos <- which(NAcol_ratio > 0.9 & NAcol < 4)
    warning(paste0("These samples have 0.9 missing rate and less than 4 values was removed first: \n",paste(colnames(quant)[pos],collapse = ", "),".\n"))
    quant <- quant[,-pos]
  }
  outliersample = 1; outlier = NULL; iter = 0;
  if (length(text.labels)!= ncol(quant) & !all(is.na(text.labels))) {
    text.labels = NA;
    warning("the number of text.labels is differ with sample number and use sample name instead of.")
  }
  if (all(is.na(text.labels))) text.labels <- colnames(quant);
  pic <- NULL # outlier sd plot will save in one file 250728
  #outliersample remove
  while (outliersample > 0){
    iter <- iter + 1;
    #remove all NA cols
    if(any(NAcol < 4)) stop("Some sample have almost all NA value.")    
    if (distmethod %in% c("correlation","bicor"))
    {
      # intra-assay correlation
      if(distmethod == "correlation") IAC <- cor(quant, use = "pairwise.complete.obs");
      #210728 add bicor function
      if(distmethod == "bicor") {
        if (!requireNamespace("WGCNA", quietly = TRUE)) {
          stop ("WGCNA in Bioconductor needed when using bicor method. Please install it.",
                call. = FALSE)}
        IAC <- WGCNA::bicor(quant, use='pairwise.complete.obs')}
      if (A.IAC) IAC <- ((1+IAC)/2)^2;
      # 1-IAC as distance
      dist <- as.dist(1 - IAC);
    } else {
      dist <- dist(t(quant), method = distmethod);
      IAC <- as.matrix(dist);
    }
    IAC[is.na(IAC)] <- 0
    # Another way to visualize outliers is to calculate the mean IAC for each array and examine this distribution
    meanIAC <- apply(IAC, 2, mean);#hist(meanIAC,breaks=10)
    sdCorr <- sd(meanIAC);#sd
    numbersd <- abs(meanIAC - mean(meanIAC)) / sdCorr;
    #if meanIAC is normal distribution,+-2 means 95%
    out_pos <- as.numeric(which(numbersd > sdout));
    out_name <- text.labels[out_pos];
    #outliers 250728
    if (plot) {
      #250728
      sddata <- data.frame(x = 1:length(numbersd),numbersd, samples = text.labels)
      sddata$lab <- ""; sddata$lab[out_pos] <- sddata$samples[out_pos]
      sddata$col <- "black"; sddata$col[out_pos] <- "red"
      sddata$pch <- 1; sddata$pch[out_pos] <- 2
      pic <- c(pic,list(ggplot(sddata,aes(y = numbersd,x=x)) +
                          geom_point(color = sddata$col,size =sddata$pch) +
                          geom_text_repel(mapping = aes(label = .data[["lab"]]))+
                          labs(title = paste0("Sample varation in outlier detection ", iter))+
                          geom_line(y = sdout,col = abline.col,linewidth = abline.lwd) +
                          theme_bw()))
      names(pic)[iter] <- paste0("Detect_",iter)
    }
    outliersample <- length(out_pos);
    outlier <- c(outlier, out_name);
    if (outliersample > 0) {
      quant <- quant[ ,-out_pos];
      text.labels <- text.labels[-out_pos];
      NArow <- rowSums(!is.na(quant))
      quant <- quant[NArow > 1,]
      inf <- data.frame(inf[NArow >1, ])
      NAcol <- colSums(!is.na(quant))
      if(any(NAcol < 10)) stop("Some sample have all NA value.") 
    }
    if (!is.na(iteration) & iter >= iteration) outliersample <- 0; 
    message("Outlier detection ",iter," times")
  }
  colnames(inf) <- colnames(data$inf)
  #remove rows with large than maxNA ratio.
  NArow <- rowSums(is.na(quant))
  NArow_ratio <- NArow/ncol(quant)
  inf$NA_num <- NArow
  pos <- which(NArow_ratio < maxNA)
  inf <- inf[pos,]
  quant <- quant[pos,]
  if (!is.na(miss.value)) quant[is.na(quant)] <- miss.value; 
  if(plot) {
    if (is.character(filename) & length(filename) == 1){
      if(!dir.exists("plot")) dir.create("plot"); #210615fix exist plot dir
      pdf(paste0("plot/",filename, " outliersample.pdf"))
    }
    for( i in names(pic))
      print(pic[[i]])
    if (is.character(filename) & length(filename) == 1) dev.off()
  }
  data <- list(inf = inf, intensity = quant, outlier =outlier)
}

