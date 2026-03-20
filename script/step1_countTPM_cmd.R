#Version: 2.0
#created: Kefu liu(liukefu19@163.com)
#Maintainer: Kefu Liu
#Date: July 30, 2025
#Software: R
#using count and TPM

min.count = 6
distmethod = "manhattan";
### 0.package check ########
existpck <- rownames(installed.packages());
if(!"optparse" %in% existpck)
  stop("no optparse package, please install the R pakcage before run me.")
if(!"tibble" %in% existpck)
  stop("no tibble package, please install the R pakcage before run me.")
if(!"edgeR" %in% existpck)
  stop("no edgeR package, please install the R pakcage before run me.")
if(!"preprocessCore" %in% existpck)
  stop("no preprocessCore package, please install the R pakcage before run me.")
options(ggalt.quiet = TRUE)
options(ggpubr.quiet = TRUE)
options(edgeR.quiet = TRUE)

### 0. command line parser ######
library("optparse")
option_list = list(
  make_option("--quantfile", action="store", default=NA,
              type='character',help="quant file path"),
  make_option("--covarfile", action="store", default=NA,
              type='character',help="meta file path"),
  make_option("--output", action="store", default=NA,
              type='character',help="output prefix name"),
  make_option("--funcdir", action="store", default =NA,
              type='character',help="preprocess function file path"),
  make_option("--sdout", action="store", default=NA,
              type='numeric',help="The number of outlier cutoff value"),
  make_option("--min.value", action="store", default=NA,
              type='numeric',help="min.value used in genefilter."),
  make_option("--NAratio", action="store", default=NA,
              type='numeric',help="NAratio used in genefilter."),
  make_option("--grpname", action="store", default=NA,
              type='character',help="group name info used in genefilter."),
  make_option("--norm.method", action="store", default=NA,
              type='character',help="norm.method: TMM or quantile"),
  make_option("--step", action="store", default=NA,
              type='character',help="1.cpm norm; 2.removeOutlier;3.genefilter;4.normalization")
)
opt = parse_args(OptionParser(option_list=option_list))
rmobj <- c(names(opt),"opt","option_list","allobj","allobj2","loginfo","existpck","funcpath")
allobj <- ls(); allobj <- allobj[!allobj %in% c("rmobj",rmobj)]

### 0.1 parameters evaluate ######
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
if(is.na(opt$sdout)) {
  opt$sdout = 3;
  message("setting sdout as default: 3")
}
if(is.na(opt$min.value)) {
  opt$min.value = 0.1;
  message("setting min.value as default: 0.1")
}
if(is.na(opt$ratio)) {
  opt$ratio = 0.25;
  message("setting ratio as default: 0.25")
}
if(is.na(opt$norm.method)) {
  opt$norm.method = "TMM"; message("Setting norm.method as default: TMM.")
}
opt$norm.method <- match.arg(opt$norm.method,c("TMM","quantile"))
if(is.na(opt$step)) {
  opt$step = "234"; 
  message("Setting step as default: 234.")
}
if(opt$step == 1) stop("step is setting to 1, should be have other steps.")

### 0.2 read quant & covar info ####
#no quantfile should have quant obj
if(is.na(opt$quantfile)) {
  if("count" %in% ls()) 
    message("using count obj to do the process.") else
      stop ("no quantfile or count object, see help info");
  if("tpm" %in% ls()) 
    message("using tpm obj to do the process.") else
      stop ("no quantfile or tpm object, see help info");
}
#check quantfile and load. RData or csv
if(!is.na(opt$quantfile)) {
  if(grepl("\\.RData$",opt$quantfile)) {
    if(file.exists(opt$quantfile)) 
      load(opt$quantfile) else 
        stop(paste0("No ",opt$quantfile," file."))
  } else stop(paste0(opt$quantfile," must be .RData file."))
  if(!"count" %in% ls()) 
    stop ("no count object, see help info");
  if(!"tpm" %in% ls()) 
    stop ("no tpm object, see help info");
}
if(!is.numeric(as.matrix(count))) 
  stop("count are not all numeric.")
if(!is.numeric(as.matrix(tpm))) 
  stop("tpm are not all numeric.")
#tpm and count should have same rownames and colnames
if(any(colnames(tpm) != colnames(count)))
  stop("tpm and count are not have matched sample.")
if(any(rownames(tpm) != rownames(count)))
  stop("tpm and count are not have matched sample.")
#check covarfile and load. (if offered) RData or csv
if(!is.na(opt$covarfile)) { #230201 updated
  if(grepl("\\.RData$",opt$covarfile)) {
    if(file.exists(opt$covarfile)) 
      load(opt$covarfile) else 
        stop(paste0("No ",opt$covarfile," file."))
  } else if (file.exists(opt$covarfile))
    covar <- read.csv(opt$covarfile) else 
      stop(paste0("No ",opt$covarfile," file."))
}
otherobj <- ls(); 
otherobj <- otherobj[!otherobj %in% c("rmobj",rmobj)]
otherobj <- otherobj[!sapply(otherobj, function(x) is.function(get(x)))]
if(!"loginfo" %in% ls()) loginfo = list();

### 0.3 obj check (covar can not exist) ################
quant <- tpm
#check when no covar, grpname should be NA 
if(!"covar" %in% ls()) {
  if(is.na(opt$grpname)) 
    message("No covar and no grpname, do gene filter will not use covar info.") else
      stop("No covar but have grpname.")
}
#if have covar, match the count&tpm
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
#have covar and check grpname
if(!is.na(opt$grpname)) {
  if(any(!opt$grpname %in% colnames(covar))) {
    stop("grpname is not in covar.")
  }
}
#make names to readable
colnames(tpm) <- colnames(count) <- make.names(colnames(tpm))
rownames(tpm) <- rownames(count) <- make.names(rownames(tpm))
rownames(covar) <- colnames(tpm)
quant <- tpm

miss.value = 0;
savedname = NULL
step1r_final_rec = NULL
if(!"./plot" %in% list.dirs()) {dir.create("plot")}
### 1. remove all zero genes #########
cat('-----RNAseq data preprocess 1----- \n')
#remove genes which have no value
s1 <- ncol(quant);  g1 <- nrow(quant);
if("covar" %in% ls()) covar_raw <- covar
#match IDinf and quant
if("IDinf" %in% ls()) {
  if(!"ID" %in% colnames(IDinf)) 
    message("input have IDinf but no ID col, not process IDinf.") else {
      pos <- match(rownames(quant),IDinf$ID)
      if(any(is.na(pos))) message("input have IDinf, but not all ID in IDinf match with quant, so IDinf have NA number.")
      IDinf <- IDinf[pos,]
      IDinf$ID <- rownames(quant)
    }
}
noNArow <- rowSums(quant != 0)
if(any(noNArow == 0)) {
  message("remove ", sum(noNArow == 0)," rows, which have less than 1 values before outlier removal.")
  quant <- quant[noNArow > 0,]
  count <- quant
  tpm <- tpm[noNArow > 0,]
}
#for plot
int_raw <- quant

### 2. removeOutlier ########

if(grepl("2",opt$step)) {
  savedname = paste0(savedname,"_rmOut");
  quant[quant ==0] <- NA
  data <- list(inf = data.frame(ID = rownames(quant)), intensity = data.frame(quant,check.names = FALSE))
  data_rmOut <- try(Outlier_Detect(data,maxNA = 1,
                                   sdout = opt$sdout,filename = opt$output,
                                   distmethod = distmethod),silent = TRUE)
  if(class(data_rmOut)=="try-error") stop("Outlier detection have some bugs.")
  quant <- data_rmOut$intensity
  quant[is.na(quant)] <- 0
  count <- count[,colnames(count) %in% colnames(quant)]
  int_rm <- quant #250730 for plot
  s2 <- ncol(quant);
  plotfile <- list.files("plot/",pattern = paste(opt$output,"outliersample"))
  loginfo1 = list(record = paste0("using Data_impute to remove outlier, sdout:",opt$sdout,", distmethod:",distmethod,
                                  "\n do ",length(plotfile)-1," times to remove the outlier.\n ",
                                  "remove ",s1-s2," samples and remain ", s2, " sample X ",g1," genes.\n"),
                  input = if(is.na(opt$covarfile)) opt$quantfile else paste(opt$quantfile,"and",opt$covarfile),
                  output = paste0(opt$output,savedname,".RData"),
                  packageVersion =paste("DDPNA",packageVersion("DDPNA")),
                  sessionInfo = sessionInfo(),
                  srcipt = "step1_rna v2.0", parameter = opt,
                  running_time = Sys.time())
  loginfo = c(loginfo, step1.1_removeOutlier = list(loginfo1))

  #if outlier detection process change sample name, it will break in step 1.
  if("covar" %in% ls()) {
    if(all(colnames(quant) %in% rownames(covar))) {
      covar = data.frame(covar[rownames(covar) %in% colnames(quant),])
      if(ncol(covar) == 1) rownames(covar) <- colnames(quant);
      if(ncol(covar) == 1) colnames(covar) = colnames(covar_raw)
    } else {
      warning("covar is not the same with quant after outlier removal. Please check the ",opt$output ,"_rmOut.RData")
      opt$step = 2
      message("Running is break in step 1.")
    }
    save(quant, covar,loginfo,list = otherobj,file =paste0(opt$output,savedname,".RData")) #230201 updated
  } else 
    save(quant, loginfo,list = otherobj, file =paste0(opt$output,savedname,".RData")) #230201 updated
  
  cat('---1. sample outlier was removed--- \n')
  cat(paste0("remove ",s1-s2," samples and remain ", s2, " sample X ",g1," genes. \n"))
  step1r_final_rec = c(step1r_final_rec, "*** sample outlier was removed ***\n", 
                       paste0("remove ",s1-s2," samples and remain ", s2, " sample X ",g1," genes. \n"))
}

### 3. gene filter count > 6 & tpm > min.value########
if(grepl("3",opt$step)){
  if(is.na(opt$grpname)) {
    #filter without grpname
    keep_tpm <- (rowSums(quant > opt$min.value)) > opt$ratio*ncol(quant)
    keep_counts <- rowSums(count >= min.count) > opt$ratio * ncol(count)
    keep <- keep_tpm & keep_counts
    cat("( value >",opt$min.value,") >",opt$ratio,"ratio in all samples will keep.\n")
  } else {
    #filter have grpname
    group <- covar[,opt$grpname];
    if(length(unique(group)) == 1) stop("group have only 1 value.")
    for(i in unique(group)) {
      keep_tpm <- (rowSums(quant[,group == i] > opt$min.value)) > opt$ratio*ncol(quant[,group == i])
      keep_counts <- rowSums(count[,group == i] >= 6) > opt$ratio * ncol(count[,group == i])
      if(i == unique(group)[1]) 
        keep <- keep_tpm & keep_counts else 
          keep <- keep | (keep_tpm & keep_counts)      
    }
    cat("( value >",opt$min.value,") >",opt$ratio,"ratio in anyone group will keep.\n")
  }
  s2 <- ncol(quant); 
  savedname = paste0(savedname,"_filter");
  quant <- quant[keep,]
  int_genefilter <- quant #for plot
  g2 <- nrow(quant)
  
  loginfo2 = list(record = paste0("Gene filter, min.count:",opt$min.value,", groupname:",opt$grpname,"\n",
                                  "remove ",g1-g2," genes and remain ", s2, " sample X ",g2," genes. \n"),
                  input = if(is.na(opt$covarfile)) opt$quantfile else paste(opt$quantfile,"and",opt$covarfile),
                  output = paste0(opt$output,savedname,".RData"),
                  sessionInfo = sessionInfo(),
                  srcipt = "step1_rna v2.0", parameter = opt,
                  running_time = Sys.time())
  loginfo = c(loginfo, step1.2_genefilter = list(loginfo2))
  #remove gene in IDinf
  if("IDinf" %in% ls()) {
    if("ID" %in% colnames(IDinf)) {
      pos <- match(rownames(quant),IDinf$ID)
      IDinf <- IDinf[pos,]
      IDinf$ID <- rownames(quant)
    }
  }
  if("covar" %in% ls())
    save(quant, covar,loginfo,list = otherobj, file =paste0(opt$output,savedname,".RData")) else #230201 updated
      save(quant,loginfo,list = otherobj, file =paste0(opt$output,savedname,".RData"))  #230201 updated
  
  cat('---2. gene filter was finished--- \n')
  cat(paste0("remove ",g1-g2," genes and remain ", s2, " sample X ",g2," genes. \n"))
  step1r_final_rec = c(step1r_final_rec, 
                       "*** gene filter was finished by filterByExpr ***\n", 
                       paste0("remove ",g1-g2," genes and remain ", s2, " sample X ",g2," genes. \n")) 
  
}

### 4. normalize############
if(grepl("4",opt$step)) {
  savedname = paste0(savedname,"_",opt$norm.method);
  if(opt$norm.method == "TMM") {
    library(edgeR)
    library(tibble)
    y <- DGEList(counts=quant)
    ##Perform TMM normalization and transfer to CPM (Counts Per Million)
    y <- calcNormFactors(y, method="TMM")
    y <- y$counts
    #count_norm=as.data.frame(cpm(y)) # fix a bug using cpm to make TMM unuse.250908
    #quant → tmm  → cpm = quant → cpm
    #quant → tmm   != quant → cpm → tmm (because quant remove some genes and totalcount is not the same.)
    quant <- y
  }
  if(opt$norm.method == "quantile"){
    library(preprocessCore)
    y <-normalize.quantiles(as.matrix(quant))
    y <- data.frame(y)
    colnames(y) <- colnames(quant)
    rownames(y) <- rownames(quant);
    quant <- y
  }
  int_norm <- quant
  loginfo3 = list(record = paste0("using",opt$norm.method,"to normalize, log2(x+1)\n",
                                  "density plot was draw and save in plot/densitty_plot_",opt$output,".pdf\n"),
                  input = if(is.na(opt$covarfile)) opt$quantfile else paste(opt$quantfile,"and",opt$covarfile),
                  output = paste0(opt$output,savedname,".RData"),
                  packageVersion =if(opt$norm.method =="TMM") paste("edgeR",packageVersion("edgeR")) else paste("preprocessCore",packageVersion("preprocessCore")),
                  sessionInfo = sessionInfo(),
                  srcipt = "step1_rna v2.0", parameter = opt,
                  running_time = Sys.time())
  loginfo = c(loginfo, step1.3_normalize = list(loginfo3))
  #save TMM result
  if("covar" %in% ls()) #230201 updated
    save(quant, covar,loginfo,list = otherobj, file =paste0(opt$output,savedname,".RData")) else
      save(quant,loginfo,list = otherobj, file =paste0(opt$output,savedname,".RData"))
  #log2 result
  minvalue = 1;
  quant <- log2(quant*signif(1/opt$min.value,digits = 2)+signif(minvalue,digits = 2))
  int_log2 <- quant
  loginfo$step1.3_normalize$output = paste0(opt$output,"_log_",opt$norm.method,".RData")
  if("covar" %in% ls()) #230201 updated
    save(quant, covar,minvalue,loginfo,list = otherobj, file = paste0(opt$output,"_log_",opt$norm.method,".RData")) else
      save(quant, minvalue,loginfo,list = otherobj, file = paste0(opt$output,"_log_",opt$norm.method,".RData"))
  cat('---3. normalization was finished. --- \n')
  cat(minvalue, "was used in log transform.\n")
  cat("Final file was ", opt$output,"_log_",opt$norm.method,".RData\n",sep = "")
  step1r_final_rec = c(step1r_final_rec, 
                       "*** normalization was finished ***\n", 
                       paste0(minvalue, "was used in log transform.\n"),
                       paste0("Final file was ", opt$output,"_log_",opt$norm.method,".RData\n")) 
}

### 5.inverse normalization transform: only for QTL analysis ############
if(grepl("5",opt$step)) {
  quant <- t(apply(quant, 1, function(x) {
    r <- rank(x, na.last = "keep", ties.method = "average")
    n <- sum(!is.na(x))
    return(qnorm(r / (n + 1)))
  }))
  loginfo4 = list(record = paste0("using inverse normalization transform to make the value have normal distribution.\n"),
                  input = if(is.na(opt$covarfile)) opt$quantfile else paste(opt$quantfile,"and",opt$covarfile),
                  output = paste0(opt$output,savedname,"_INT.RData"),
                  sessionInfo = sessionInfo(),
                  srcipt = "step1_rna v2.0", parameter = opt,
                  running_time = Sys.time())
  loginfo = c(loginfo, step1.4_normalize = list(loginfo4))
  int_INT <- quant
  #save TMM result
  if("covar" %in% ls()) #230201 updated
    save(quant, covar,loginfo,list = otherobj, file =paste0(opt$output,savedname,"_INT.RData")) else
      save(quant,loginfo,list = otherobj, file =paste0(opt$output,savedname,"_INT.RData"))
  cat('---4. INT was finished. It is only for QTL analysis. --- \n')
}

### PCA #####
#  genefilter process will not draw pca
if("covar" %in% ls()) {
  if(opt$step != "3") {
    plotpdfname = "_step1"
    pca_list <- list(raw =int_raw,
                     rm = if(grepl("2",opt$step)) int_rm,
                     norm = if(grepl("4",opt$step)) int_norm,
                     INT = if(grepl("5",opt$step)) int_INT)
    if(!is.null(covar)) {
      cat("- Draw PCA plot - \n")
      for(i in names(pca_list)) {
        if(!is.null(pca_list[[i]])) pca_list[[i]] <- pca_dat(pca_list[[i]])
      }
      message("PCA is calcluated.")
      if(grepl("1|2",opt$step))
        save(covar,covar_raw,pca_list,file = paste0(opt$output,plotpdfname,"_pca.RData")) else
          save(covar,pca_list,file = paste0(opt$output,plotpdfname,"_pca.RData"))
      message("pca result was saved.")
      pdf(paste0("plot/",opt$output,plotpdfname,"_pcaplot.pdf"), width = 10, height = 8)
      for(i in names(pca_list)) {
        pcadata <- pca_list[[i]]
        if(!is.null(pcadata)) {
          if(inherits(pcadata,"prcomp")) {
            if(i %in% c("raw")) pca_plot_var(pcadata, covar_raw, title = i)
            if(i %in% c("rm","norm","INT")) pca_plot_var(pcadata, covar, title = i)
          } else 
            warning(paste0(i, "step data run pca have error."))
        }
      }
      dev.off()
      message("pca plot was finished.") 
    } else cat(" - covar is NULL, will not draw PCA plot.- \n")
  } else cat("- only do gene filter, will not draw PCA plot. - \n")
} else cat(" - No covar, will not draw PCA plot.- \n")

### boxplot #######
cat("- Draw boxplot/density plot - \n")
# random pick 50 samples to draw boxplot, rmOutlier will not darw
int <- quant
if(ncol(int) > 50) {message("Samples large than 50, random pick 50 sample to draw intensity distribution.")
  pos <- sample(1:ncol(int),50)
  int_raw2 <- int_raw[,pos]
  if(grepl("3",opt$step)) int_genefilter2 <- int_genefilter[,pos]
  if(grepl("4",opt$step)) int_norm2 <- int_norm[,pos]
  if(grepl("5",opt$step)) int_INT2 <- int_INT[,pos]
  width = ceiling(50/5);
} else {
  width = 10
  int_raw2 <- int_raw
  if(grepl("3",opt$step)) int_genefilter2 <- int_genefilter
  if(grepl("4",opt$step)) int_norm2 <- int_norm
  if(grepl("5",opt$step)) int_INT2 <- int_INT
}
pdf(paste0("plot/",opt$output,"_sample&gene distribution.pdf"),width = width,height = 6)
boxplot(log2(int_raw2+1), main = "raw intensity distribution in sample",ylab = "log2 value");
if(grepl("3",opt$step)) boxplot(log2(int_genefilter2+1), main = "gene filter data distribution in sample",ylab = "log2 value");
if(grepl("4",opt$step)) boxplot(log2(int_norm2+1), main = "normalization data distribution in sample",ylab = "log2 value");
if(grepl("5",opt$step)) boxplot(int_INT2, main = "INT normalization data distribution in sample",ylab = "value");

try(densityplot(int_raw,"raw data",savefile = FALSE,verbose = 1))
if(grepl("3",opt$step)) try(densityplot(int_genefilter,"gene filter data",savefile = FALSE))
if(grepl("4",opt$step)) try(densityplot(int_norm,"normalization data",savefile = FALSE))
dev.off()

cat("\n")
cat(step1r_final_rec)
cat("--- step1_countTPM all finished ---\n")

