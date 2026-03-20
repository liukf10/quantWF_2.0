#Version: 2.0
#created: Kefu liu(liukefu19@163.com)
#Maintainer: Kefu Liu
#Date: Sep 21, 2025
#Software: R


#covar can be NULL
#quantfile must be RData or csv file. 
### package check ########
existpck <- rownames(installed.packages());
if(!"optparse" %in% existpck)
  stop("no optparse package, please install the R pakcage before run me.")
if(!"lme4" %in% existpck)
  stop("no lme4 package, please install the R pakcage before run me.")

### command line parser######
library("optparse")
option_list = list(
  make_option("--quantfile", action="store", default=NA,
              type='character',help="quant file path"),
  make_option("--covarfile", action="store", default=NA,
              type='character',help="meta file path"),
  make_option("--genotype_pcfile", action="store", default=NA,
              type='character',help="genotype_pc file path"),
  make_option("--IDfile", action="store", default=NA,
              type='character',help="IDinf file path"),
  make_option("--output", action="store", default=NA,
              type='character',help="output prefix name"),
  make_option("--funcdir", action="store", default=NA,
              type='character',help="preprocess function file path"),
  make_option("--method", action="store", default=NA,
              type='character',help="estimate hidden factor method: pca or peer, default is pca."),
  make_option("--svamethod", action="store", default="be",
              type='character',help="sva method: be or leek, default is be."),
  make_option("--nFactor", action="store", default= NA,
              type='numeric',help="nFactor in peer, default is depend on sample size."),
  make_option("--nPCAchoose", action="store", default= NA,
              type='character',help="PCAforQTLs PC number choose, BE or Elbow method, default is BE."),
  make_option("--evaluate_method", action="store", default=NA,
              type='character',help="evaluate the result: pca,pvca and bic, default is pvca&pca."),
  make_option("--BIC_cores", action="store", default= NA,
              type='numeric',help="parallel computing cores, default is 30."),
  make_option("--tools", action="store", default= NA,
              type='character',help="tools for QTL analysis: tensorQTL or QTLtools. QTLtools have Genebody and TSS, you can use GB/TSS/all to refer, Setting tools as default: tensorQTL.")
)

opt = parse_args(OptionParser(option_list=option_list))
rmobj <- c(names(opt),"opt","option_list","allobj","allobj2","loginfo","existpck","funcpath")
allobj <- ls(); allobj <- allobj[!allobj %in% c("rmobj",rmobj)]

### 0.1 parameters evaluate ######
if(!dir.exists("plot")) dir.create("plot")
if(is.na(opt$output)) {
  opt$output = "data";
  message("no output, setting to 'data' as prefix")
}
if(is.na(opt$funcdir)) {
  warning('funcdir is setting to NA, means preprocess_func.R in the default dir');
  if(file.exists("preprocess_func.R")) 
    source("preprocess_func.R") else 
      stop("No preprocess_func.R file in fundir path.")
} else {
  if(tail(unlist(strsplit(opt$funcdir,"")),1) == "/"){
    funcpath = paste0(opt$funcdir,"preprocess_func.R");
  } else {
    opt$funcdir <-  paste0(opt$funcdir,"/")
    funcpath = paste0(opt$funcdir,"preprocess_func.R");}
  if(file.exists(funcpath)) 
    source(funcpath) else
      stop("No preprocess_func.R file in fundir path.")
}
if(is.na(opt$method)) {
  opt$method = "pca"; message("Setting method as default: pca.")}
opt$method <- match.arg(opt$method,c("peer","pca"));
if(opt$method == "peer"){
  if(!"peer" %in% existpck)
    stop("no peer package, peer only install in linux system.")
}
if(is.na(opt$nPCAchoose)) {
  opt$nPCAchoose = "BE"; message("Setting evaluate_method as default: BE.")}
opt$nPCAchoose <- match.arg(opt$svamethod,c("BE","Elbow"));
if(is.na(opt$evaluate_method)) {
  opt$evaluate_method = "pvca&pca"; message("Setting evaluate_method as default: pvca&pca.")}
if(is.na(opt$BIC_cores)) {
  opt$BIC_cores = 30; message("Setting BIC_cores as default: 30.")}
if(is.na(opt$tools)) {
  opt$tools = "tensorQTL"; message("tools for QTL analysis: tensorQTL or QTLtools.  QTLtools have Genebody and TSS, you can use GB/TSS/all to refer, Setting tools as default: tensorQTL.")}
opt$tools <- match.arg(opt$tools,c("tensorQTL","QTLtools","TSS.QTLtools","GB.QTLtools","all.QTLtools","both"));

# opt$nfactor will check in peer section
### 0.2 read quant & covar&genotype_pc&IDinf info ####
#no quantfile should have quant obj
if(is.na(opt$quantfile)) {
  if("quant" %in% allobj) 
    message("using quant obj to do the process.") else
      stop ("no quantfile or quant object, see help info");
}
#check quantfile and load. RData or csv
if(!is.na(opt$quantfile)) {
  if(grepl("\\.RData$",opt$quantfile)) {
    if(file.exists(opt$quantfile)) 
      load(opt$quantfile) else 
        stop(paste0("No ",opt$quantfile," file."))
  } else if (file.exists(opt$quantfile))
    quant <- read.csv(opt$quantfile, row.names = 1) else 
      stop(paste0("No ",opt$quantfile," file."))
  if(!"quant" %in% ls()) 
    stop ("no quant object, see help info");
}
if(!is.numeric(as.matrix(quant))) 
  stop("quant are not all numeric.")
#check covarfile and load. (if offered) RData or csv
if(!is.na(opt$covarfile)) { #230201 updated
  if(grepl("\\.RData$",opt$covarfile)) {
    if(file.exists(opt$covarfile)) 
      load(opt$covarfile) else 
        stop(paste0("No ",opt$covarfile," file."))
  } else if (file.exists(opt$covarfile))
    covar <- read.table(opt$covarfile) else 
      stop(paste0("No ",opt$covarfile," file."))
}
#check genotype_pcfile and load. (if offered) RData or csv
if(!is.na(opt$genotype_pcfile)) { 
  if(grepl("\\.RData$",opt$genotype_pcfile)) {
    if(file.exists(opt$genotype_pcfile)) 
      load(opt$genotype_pcfile) else 
        stop(paste0("No ",opt$genotype_pcfile," file."))
  } else if (file.exists(opt$genotype_pcfile))
    genotype_pc <- read.table(opt$genotype_pcfile) else 
      stop(paste0("No ",opt$genotype_pcfile," file."))
}
if(!is.na(opt$IDfile)) { 
  if (file.exists(opt$IDfile)) {
    #gtf file only extract gene
    if(grepl("\\.RData$",opt$IDfile)) {
      load(opt$IDfile)
    } else if(grepl("gtf",opt$IDfile)) {
      if(!"rtracklayer" %in% existpck)
        stop("no rtracklayer package, gtf file read need rtracklayer.")
      suppressPackageStartupMessages(library(rtracklayer))
      IDinf <- import(opt$IDfile)
      IDinf <- as.data.frame(IDinf)
      IDinf <- subset(IDinf,type %in% "gene")
      IDinf$ID <- gsub("\\..*","",IDinf$gene_id)
      IDinf$ID[duplicated(IDinf$ID)] <- paste0(IDinf$ID[duplicated(IDinf$ID)],"_1")
      for( i in ncol(IDinf):1) {
        if(length(unique(IDinf[,i])) == 1) IDinf <- IDinf[,-i]
      }
    } else IDinf <- read.table(opt$IDfile) 
  } else stop(paste0("No ",opt$IDfile," file."))
}

otherobj <- ls(); 
otherobj <- otherobj[!otherobj %in% c("rmobj",rmobj)]
otherobj <- otherobj[!sapply(otherobj, function(x) is.function(get(x)))]
if(!"loginfo" %in% ls()) loginfo = list();

### 0.3 obj check ###########
if(!"covar" %in% ls()) stop ("no covar object, please check the data.");
if(!"genotype_pc" %in% ls()) stop ("no genotype_pc object, please check the data.");
if(!"IDinf" %in% ls()) stop ("no IDinf object, please check the data.");

#covar quant match
if("covar" %in% ls()) {
  if(!is.null(covar)) {
    #when covar is a vector will not match and transfer to data.frame 
    if(is.vector(covar)) {
      if(length(covar) == ncol(quant)) {
        covar <- data.frame(covar);
        colnames(covar) = "group" 
        rownames(covar) = colnames(quant);
        warning("covar is vector, please make sure the order is consistent with quant colnames.")
      } else stop("covar is vector and the number is not the same with quant column number.")
    } else if(length(dim(covar)) == 2) 
    { #quant and covar match one-by-one
      for (i in 1:(ncol(covar)+1)) {
        if(i == ncol(covar)+1) {
          if(all(colnames(quant) %in% rownames(covar)))
            break else i = 0;
        }
        if(all(colnames(quant) %in% covar[,i]))
          break
      }
      if(i == 0) stop("Not all sample match with covar.")
      if(i == (ncol(covar)+1)) {
        pos <- match(colnames(quant),rownames(covar))
        #230310 
        if(ncol(covar) == 1) {
          keepcolvar <- colnames(covar)
          covar <- data.frame(covar[pos,])
          colnames(covar) <- keepcolvar;
          rownames(covar) <- colnames(quant);
        } else covar <- covar[pos,];
      } else {
        pos <- match(colnames(quant),covar[,i]);
        covar <- covar[pos,]; covar <- covar[-i]; 
        rownames(covar) <- colnames(quant);
      } #covar & quant are matched 
    }
  }
}
#check and remove single value col in covar
if(!is.null(covar)) {
  remove_covar_col = NULL
  for( i in 1:ncol(covar)) {
    if(length(unique(covar[,i])) == 1) 
      remove_covar_col = c(remove_covar_col, i) 
  }
  if(!is.null(remove_covar_col)) {
    remove_covar_name = colnames(covar)[remove_covar_col];
    covar <- covar[,-remove_covar_col];
    message(paste0(paste0(remove_covar_name,collapse = ",")," have removed because only one element. \n"))
  }
}
#check if NA in covar
if(!is.null(covar)) {
  for( i in 1:ncol(covar))if(any(is.na(covar[,i]))) stop("NA not allowed in covar.") 
}
#genotype_pc covar match and merge 
colnames(genotype_pc) <- paste0("genopc", 1:ncol(genotype_pc))
pos <- match(rownames(covar),rownames(genotype_pc))
if(any(is.na(pos))) stop("Not all sample in covar match with genotype_pc file.")
genotype_pc <- genotype_pc[pos,]
covar<- cbind(covar,genotype_pc)
#IDinf check
#match IDinf and quant
if(!"ID" %in% colnames(IDinf)) 
  stop("IDinf should have ID col.") else {
    pos <- match(rownames(quant),IDinf$ID)
    if(any(is.na(pos))) stop("Not all ID in IDinf match with quant.")
    IDinf <- IDinf[pos,]
  }
if(!"./plot" %in% list.dirs()) {dir.create("plot")}
### 1. hidden covariate search #########
cat('-----data preprocess 2 hidden covariate search----- \n')

if(opt$method == "peer") {
  if(is.na(opt$nFactor))  {
    if(ncol(quant) < 150) 
      nFactor = 15 else if (ncol(quant) < 250)
        nFactor = 30 else if (ncol(quant) < 350)
          nFactor = 45 else nFactor = 60
          message("Setting nFactor depend on sample size (GTEx suggestion): ", nFactor)
  } else nFactor = opt$nFactor
  if(nFactor < 1) stop("wrong nFactor.")
  suppressPackageStartupMessages(library(peer))
  expr = t(as.matrix(quant)) 
  model = PEER()  # create the model object
  invisible(PEER_setPhenoMean(model,expr))# set the observed data
  invisible(PEER_setNk(model,nFactor)) # gradient number of factors
  invisible(PEER_getNk(model))
  invisible(PEER_setAdd_mean(model, TRUE))  # include an additional factor (covariate) to account for the mean expression
  #covar1 = NULL 1.obs+hid and no keep; 2.covar is NULL
  if(!is.null(covar)) {
    invisible(PEER_setCovariates(model, as.matrix(covar)))  # adding covariates has no effect on the model?
  }
  #PEER_setNmax_iterations(model, 100)  # If the model is not converged after the default 1,000 iterations, and the variance of the residuals keeps decreasing, choose a higher value of iterations, e.g., 10,000.
  time = system.time(PEER_update(model))  # perform the inference
  
  factors = PEER_getX(model)  # inferred confounders samples x PEER factors
  factors=factors[,-1];
  rownames(factors) = rownames(expr);
  colnames(factors) = paste0("peer",1:ncol(factors))
  weights = PEER_getW(model)  # their weights
  precision = PEER_getAlpha(model)     # precision (inverse variance) of the weights
  residuals = PEER_getResiduals(model) # the residual dataset
  #plot(precision)
  #PEER_plotModel(model)
  #Variance_factor_plot
  {
    Variance = apply(factors,2,var);
    Variance<-sort(Variance, decreasing=TRUE); 
    Variance<-100*Variance/sum(Variance)
    pdf(paste0("plot/peer_variance_factor_",opt$output,".pdf"),width=5,height=5)
    plot(Variance,type="o",pch=19,col="red",xlab="Factors")
    dev.off()
  }
  if(is.null(covar))
    covar_hid <- data.frame(factors) else
      covar_hid <- data.frame(covar,factors)
}
if(opt$method == "pca"){
  expr = t(as.matrix(quant)) 
  prcompResult<-prcomp(expr,center=TRUE,scale.=TRUE)
  #This should take less than a minute.
  PCs<-prcompResult$x 
  dim(PCs)
  importanceTable<-summary(prcompResult)$importance
  PVEs<-importanceTable[2,]
  #sum(PVEs) #Theoretically, this should be 1.
  #pdf("PVEs.pdf")
  #plot(PVEs,xlab="PC index",ylab="PVE")
  resultRunElbow<-runElbow(prcompResult=prcompResult)
  #print(resultRunElbow)
  K_elbow<-resultRunElbow
  resultRunBE<-runBE(expr,B=25,alpha=0.05,mc.cores=1)
  #print(resultRunBE$numOfPCsChosen)
  K_BE<-resultRunBE$numOfPCsChosen 
  pdf(paste0("plot/PCA_variance_factor_",opt$output,".pdf"),width=5,height=5)
  makeScreePlot(prcompResult,labels=c("Elbow","BE"),values=c(K_elbow,K_BE),
                titleText=opt$output)
  dev.off()
  if(opt$nPCAchoose=="BE") PCsTop<-PCs[,1:K_BE] #default is choose K_BE
  if(opt$nPCAchoose=="Elbow") PCsTop<-PCs[,1:K_elbow]
  if(!is.null(covar)){
    if(ncol(covar) != 1) {
      covar2 <- covar;
      for (i in 1:ncol(covar2)) {
        if(!class(covar2[,i]) %in% c("numeric","integer"))
          covar2[,i] <- as.numeric(as.factor(covar2[,i]))
      }
    }
    if(ncol(covar) == 1) {
      covar2 <- covar;
      if(!inherits(covar[,1],c("numeric","integer")))
        covar2[,1] <- as.numeric(as.factor(covar2[,1]))
    }
    knowncovFiltered<-filterKnownCovariates(covar2,PCsTop,unadjustedR2_cutoff=0.9)
    PCsTop<-scale(PCsTop)
    covar_hid <-data.frame(covar[,colnames(covar) %in% colnames(knowncovFiltered)],PCsTop)
  }
} 

write.table(t(covar_hid),file=paste0(opt$output,"_",opt$method,"_covariate4QTL.txt"),col.names=NA,quote=F)
record <- paste0(opt$output,"_",opt$method,"_covariate4QTL.txt was saved.")
message(record)
record <- paste0(record,"\n")
### 2 merge info ############

IDinf$seqnames <- gsub("chr","",IDinf$seqnames)
pos <- IDinf$seqnames %in% 1:22
quant <- quant[pos,]
IDinf <- IDinf[pos,]
IDinf$TSS <- IDinf$start
IDinf$TSS[IDinf$strand == "-"] <- IDinf$end[IDinf$strand == "-"]
IDinf$TSSstart = IDinf$TSS -1;
IDinf$TSSend = IDinf$TSS;
IDinf <- IDinf[,c("seqnames","TSSstart","TSSend","gene_id","strand","start","end")]
IDinf$seqnames <- as.numeric(IDinf$seqnames)
colnames(IDinf)[1] <- "#chr"

###For tensorQTL 
if(grepl("tensorQTL|both",opt$tools)) {
  tssbed <- cbind(IDinf[,1:4],quant)
  colnames(tssbed)[2:3] <- c("start","end")
  tssbed <- tssbed[order(tssbed$start),]
  tssbed <- tssbed[order(tssbed$`#chr`),]
  
  # Write TSS-based BED file
  write.table(
    tssbed,
    file = paste0(opt$output,".tss.tensor.bed"),
    quote = FALSE, sep = "\t",row.names = FALSE, col.names = TRUE)
  
  record <- paste0(record, opt$output,".tss.tensor.bed was saved for tensorQTL\n")
  message("Done! The file 'tss.tensor.bed' for tensorQTL have been generated successfully.")
}
###  For QTltools 
if(grepl("QTLtools|both",opt$tools)) {
  if(grepl("tss|all|both",opt$tools)|opt$tools == "QTLtools") {
    message("Creating TSS-based BED file...")
    tssbed <- cbind(IDinf[,c(1,2,3,4,4,5)],quant)
    colnames(tssbed)[2:5] <- c("start","end","PID","GID") 
    tssbed <- tssbed[order(tssbed$start),]
    tssbed <- tssbed[order(tssbed$`#chr`),]
    
    # Write TSS-based BED file
    write.table(
      tssbed,
      file = paste0(opt$output,".tss.bed"),
      quote = FALSE,sep = "\t",row.names = FALSE,col.names = TRUE)
    record <- paste0(record, opt$output,".tss.bed was saved for QTLtools\n")
    message("Done! The file 'expression.tss.bed' for QTLtools have been generated successfully.")
  }
  if(grepl("GB|all|both",opt$tools)) {
    message("Creating Genebody-based BED file...")
    genebody_bed <- cbind(IDinf[,c(1,6,7,4,4,5)],quant)
    colnames(genebody_bed)[4:5] <- c("PID","GID") 
    #Write Genebody-based BED 
    write.table(
      genebody_bed,
      file = paste0(opt$output,".genebody.bed"),
      quote = FALSE,sep = "\t", row.names = FALSE,col.names = TRUE)
    record <- paste0(record, opt$output,".genebody.bed was saved for QTLtools\n")
    message("Done! The file 'expression.genebody.bed' for QTLtools have been generated successfully.")
  } 
}

###  covariants correlation #######

if(!is.null(covar_hid)){
  if(ncol(covar_hid) != 1) {
    covar_hid2 <- covar_hid;
    for (i in 1:ncol(covar_hid2)) {
      if(!class(covar_hid2[,i]) %in% c("numeric","integer"))
        covar_hid2[,i] <- as.numeric(as.factor(covar_hid2[,i]))
    }
  }
}
if(!is.null(covar_hid)) {
  cat("- Draw covar correlation plot - \n")
  #plot 
  suppressPackageStartupMessages(require(corrplot))
  col2 <- colorRampPalette(c("#053061","#2166AC","#4393C3","#92C5DE",
                             "#D1E5F0", "#FFFFFF","#FDDBC7","#F4A582",
                             "#D6604D","#B2182B","#67001F"))
  pdf(file = paste0("plot/",opt$output,"_",opt$method,"_covhid_cor.pdf"),
      width = 10, height = 10)
  if(!is.null(covar_hid))
    if(ncol(covar_hid) != 1)
      corrplot(cor(covar_hid2),p.mat = cor.mtest(covar_hid2)$p, col = col2(200),
               insig = "p-value", order = "original",method = "ellipse",
               type = "lower", tl.pos = "d")
  dev.off()
}
###  pvca  ########
if(grepl("pvca|PVCA",opt$evaluate_method)) {
  if("lme4" %in% existpck) { #pvca
    if(!is.null(covar)) {
      cat("- Draw pvca plot - \n")
      if(ncol(covar) > 0) {
        pdf(paste0("plot/",opt$output,"_",opt$method,"_pvca.pdf"),width = 12,height = 6)
        suppressMessages(pvcaresult <- pvca_evaluate(quant,covar_hid2, pct_threshold = 0.6,output = "cov4QTL"))
        print(pvcaresult$pic)
        dev.off()
        message("pvca plot was finished.") 
      }
    }
  } else  cat(" - no lme4 pcakge, not draw PVCA plot.- \n")
}

### BIC #######
if(grepl("bic|BIC",opt$evaluate_method)) {
  cat("- Draw BIC plot - \n")
  pdf(paste0("plot/",opt$output,"_",opt$method,"_bic.pdf"),width = 12,height = 6)
  BICresult <- BIC_evaluate(quant,covar_hid2,num_cores = opt$BIC_cores,output = paste0("BIC estimation"))
  print(BICresult$pic)
  dev.off()
  message("BIC plot was finished.") 
}
cat(record)
cat("----step2 cov4QTL runing finished ------\n")