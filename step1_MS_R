#Version: 2.0
#created: Kefu liu(liukefu19@163.com)
#Maintainer: Kefu Liu
#Date: July 26, 2025
#Software: R

distmethod = "manhattan";
### package check ########
existpck <- rownames(installed.packages());
if(!"DDPNA" %in% existpck)
  stop("no DDPNA package, please install the R pakcage before run me.")
if(!"Biobase" %in% existpck)
  stop("no Biobase package, please install the R pakcage before run me.")
if(!"impute" %in% existpck)
  stop("no impute package, please install the R pakcage before run me.")



### object check ##################
allobj <- ls()
if(!"quantfile" %in% allobj )
  quantfile = NA;
if(!"covarfile" %in% allobj )
  covarfile = NA;
if(!"output" %in% allobj)
  output = NA;
if(!"funcdir" %in% allobj)
  funcdir = NA;
if(!"sdout" %in% allobj)
  sdout = NA;
if(!"maxNA" %in% allobj)
  maxNA = NA;
if(!"step" %in% allobj)
  step = NA;
if(!"value" %in% allobj)
  value = NA;

cat('-----MS-based data preprocess 1----- \n')
opt = list(quantfile = quantfile, covarfile = covarfile,
           output = output, funcdir = funcdir, 
           sdout = sdout, maxNA = maxNA, step = step,
           value = value);
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
if(is.na(opt$sdout)) {
  opt$sdout = 3; message("Setting sdout to default: 3")
}
if(length(opt$maxNA) > 1) stop("maxNA should be one numeric.")
if(is.na(opt$maxNA)) {
  opt$maxNA = 0.8; message("Setting maxNA to default: 0.8")
}
if(!is.numeric(opt$maxNA)) stop("maxNA must be numeric.")
if(opt$maxNA > 1 & opt$maxNA < 0) {stop("maxNA must set beteween 0 to 1.") }
if(is.na(opt$step)) {
  opt$step = "1234"; message("Setting step to default: 1234")
  message("1.norm by all expression; 2.remove outlier; 3.filter high missing rows; 4.missing value imputation")
}
if(opt$step == 1) stop("step is setting to 1, should be have other steps.")
if(is.na(opt$value)) {
  opt$value = "value"; message("Setting value to default: value")
  message("value: using log2 value to run PCA and boxplot; ratio: using ratio to run PCA and boxplot.")
}
opt$value <- match.arg(opt$value, c("value","ratio"))
suppressPackageStartupMessages(require(DDPNA))

### 0.2 read quant & covar info ####
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
otherobj <- ls();
otherobj <- otherobj[!otherobj %in% c("rmobj",rmobj)]
otherobj <- otherobj[!sapply(otherobj, function(x) is.function(get(x)))]
if(!"loginfo" %in% ls()) loginfo = list();

### 0.3 obj check (covar can not exist) ################
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

step1p_final_rec = NULL
miss.value = NA;
if(!"IDinf" %in% ls()) stop("no IDinf, it seem not MS data.")
if(nrow(IDinf) != nrow(quant)) stop("IDinf and quant row number are not the same.")
if(!"./plot" %in% list.dirs()) {dir.create("plot")}
### 1. data norm ######
cat('-----MS data preprocess 1----- \n')
if("data_rmOut_imp" %in% otherobj) quant <- data_rmOut_imp$intensity
s1 <- ncol(quant);  g1 <- nrow(quant);
cat(paste0("raw data have ", s1, " sample X ", g1," rows \n"))
#remove genes which have no value
noNArow <- rowSums(!is.na(quant))
if(any(noNArow != 0)) {
  message("remove ", sum(noNArow == 0)," rows, which have no values before outlier removal.")
  quant <- quant[noNArow > 0,]
  IDinf <- IDinf[noNArow > 0,]
  if("data_rmOut" %in% otherobj) {
    data_rmOut$inf <- data_rmOut$inf[noNArow > 0,];
    data_rmOut$intensity <- data_rmOut$intensity[noNArow > 0,]; }
}
if(opt$value == "value") quant_raw <- log2(quant) else quant_raw = quant #recored the step res
IDinf_raw <- IDinf
#do norm   also rawdata
if(grepl("1",opt$step)) {
  cat("--- norm all ---\n")
  norm.factor <- colSums(quant, na.rm = TRUE);
  norm.factor <- norm.factor/mean(norm.factor);
  quant <- data.frame(t(t(quant)/norm.factor), check.names = FALSE)
  colnames(quant) <- colnames(quant_raw) #240313 modified to avoid change colname
  cat("---norm all was finished.--- \n")
  step1p_final_rec = c(step1p_final_rec," *** norm all  ***\n")
  if(opt$value == "value") quant_norm <- log2(quant) else quant_norm = quant #recored the step res
}
#remove rows which have 0,1,2 value
NArow <- rowSums(!is.na(quant))
if(any(NArow < 3)) {
  message("remove ", sum(NArow < 3)," rows, which have less than 3 values before outlier removal.")
  quant_raw <- quant_raw[NArow > 2,]
  IDinf_raw <- IDinf_raw[NArow > 2,]
  quant <- quant[NArow > 2,]
  IDinf <- IDinf[NArow > 2,]
  if("data_rmOut" %in% ls()) {
    data_rmOut$inf <- data_rmOut$inf[NArow > 2,];
    data_rmOut$intensity <- data_rmOut$intensity[NArow > 2,]; }
  if(grepl("1",opt$step)) if(opt$value == "value") quant_norm <- log2(quant) else quant_norm = quant #recored the step res
}

### 2. remove outlier ######
# covar will follow with outlier detect
if("covar" %in% ls()) covar_raw <- covar
if(grepl("2",opt$step)) {
  #remove outlier
  data_rmOut = try(Outlier_Detect(list(inf = IDinf,intensity = quant), maxNA = 1,
                                  sdout = opt$sdout,filename = opt$output,
                                  distmethod = distmethod, A.IAC = FALSE),silent = TRUE)
  if(class(data_rmOut)=="try-error") {
    message(data_rmOut)
    stop("Outlier detection have some bugs.")}
  #check colNAratio
  NAcol <- colSums(!is.na(data_rmOut$intensity))
  NAcolratio <- 1-(NAcol/nrow(data_rmOut$intensity))
  if(max(NAcolratio) > 0.8) {
    message("After remove outlier, some samples have 80% missing rate.")
  }
  #after remove some samples, may some rows have all missing value
  #remove rows which have 0,1,2 value
  NArow <- rowSums(!is.na(data_rmOut$intensity))
  if(any(NArow < 3)) {
    message("remove ", sum(NArow < 3)," rows, which have less than 3 values after outlier removal.")
    data_rmOut$intensity <- data_rmOut$intensity[NArow > 2,]
    data_rmOut$inf <- data_rmOut$inf[NArow > 2, ]
  }
  #quant save in two types
  if(opt$value == "value") quant <- log2(data_rmOut$intensity) else quant = data_rmOut$intensity
  IDinf <- data_rmOut$inf
  s2 = ncol(quant); g2 <- nrow(quant)
  if(grepl("1",opt$step)) record = "do normalization by total intensity.\n" else record = NULL
  loginfo2 = list(record = paste0(record, "Remove outlier, sdout:",opt$sdout,", distmethod:",distmethod,
                                  "\n remove ",s1-s2," samples and remain ", s2, " sample X ",g2," rows.\n"),
                  input = if(is.na(opt$covarfile)) opt$quantfile else paste(opt$quantfile,"and",opt$covarfile),
                  output = paste0(opt$output,"_rmOut.RData"),
                  packageVersion =paste("DDPNA",packageVersion("DDPNA")),
                  sessionInfo = sessionInfo(),
                  srcipt = "step1_MS V2.0", parameter = opt,
                  running_time = Sys.time())
  loginfo = c(loginfo, step1_removeOutlier = list(loginfo2))
  if("covar" %in% ls()) {
    if(all(colnames(quant) %in% rownames(covar))) {
        covar = data.frame(covar[rownames(covar) %in% colnames(quant),])
        if(ncol(covar) == 1) rownames(covar) <- colnames(quant);
        if(ncol(covar) == 1) colnames(covar) = colnames(covar_raw)
      } else {
        warning("covar is not the same with quant after outlier removal. Please check the ",opt$output ,"_rmOut.RData")
        opt$step = 2
        message("Running is break in step 2.")
      }
    covar_rm <- covar
    save(covar,quant,IDinf,data_rmOut,loginfo,list = otherobj, file = paste0(opt$output,"_rmOut.RData"))
  } else save(quant,IDinf,data_rmOut,loginfo,list = otherobj, file = paste0(opt$output,"_rmOut.RData"))
  quant_rm <- quant
  IDinf_rm <- IDinf
  cat('---sample outlier was removed--- \n')
  cat(paste0("remove ",s1-s2," samples and remain ", s2, " sample X ",g2," rows \n"))
  step1p_final_rec = c(step1p_final_rec,
                         "*** sample outlier was removed  ***\n",
                         paste0("remove ",s1-s2," samples and remain ", s2, " sample X ",g2," rows \n"))
} else if(!"data_rmOut" %in% ls()) data_rmOut <- list(inf = IDinf,intensity = quant)
### 3. gene filter ######
if(grepl("3",opt$step)) {
  if(!"data_rmOut" %in% ls())
    stop ("no data_rmOut object, cannot do gene filter");
  if(!"s2" %in% ls()) {
    s2 = ncol(data_rmOut$intensity); g2 <- nrow(data_rmOut$intensity);
  }
  #remove rows based on maxNA
  NArow <- rowSums(is.na(quant))
  NArowratio <- NArow/ncol(quant)
  quant <- quant[NArowratio < opt$maxNA,]
  data_rmOut$intensity <- data_rmOut$intensity[NArowratio < opt$maxNA,]
  IDinf <- data.frame(IDinf[NArowratio < opt$maxNA,])
  colnames(IDinf) <- colnames(data_rmOut$inf)
  data_rmOut$inf <- IDinf
  s3 = ncol(quant); g3 <- nrow(quant)
  loginfo3 = list(record = paste0("Remove low expression rows, maxNA:",opt$maxNA,
                                  "remove ",g2-g3," rows and remain ", s2, " sample X ",g3," rows \n",
                                  "sample abundance plot, density plot and pca plot is saved in plot dir to check the data distribtuion"),
                  input = if(is.na(opt$covarfile)) opt$quantfile else paste(opt$quantfile,"and",opt$covarfile),
                  output = paste0(opt$output,"_rmOut_filter.RData"),
                  packageVersion =paste("DDPNA",packageVersion("DDPNA")),
                  sessionInfo = sessionInfo(),
                  srcipt = "step1_MS V2.0", parameter = opt,
                  running_time = Sys.time())
  loginfo = c(loginfo, step1_gene_filter = list(loginfo3))
  #save quant info
  if("covar" %in% ls()) {
    covar_rm <- covar
    save(covar,quant,IDinf,data_rmOut,loginfo,list = otherobj, file = paste0(opt$output,"_rmOut_filter.RData"))
  } else save(quant,IDinf,data_rmOut,loginfo,list = otherobj, file = paste0(opt$output,"_rmOut_filter.RData"))
  quant_filter <- quant
  IDinf_filter <- IDinf

  cat('---gene filter was finished.--- \n');
  cat(paste0("remove ",g2-g3," rows and remain ", s3, " sample X ",g3," rows \n"))
  step1p_final_rec = c(step1p_final_rec,
                       "*** high missing rows were removed  ***\n",
                       paste0("remove ",g2-g3," rows and remain ", s3, " sample X ",g3," rows \n"))
} else if(!"data_rmOut" %in% ls()) data_rmOut <- list(inf = IDinf,intensity = quant)

### 4. imputation ######
if(grepl("4",opt$step)) {
  if(!"data_rmOut" %in% ls())
    stop ("no data_rmOut object, cannot do imputation.");
  data_rmOut_imp = Data_impute(data_rmOut, splNExt=FALSE, intensity = "intensity",
                               miss.value = miss.value, maxNAratio = opt$maxNA,
                               removeOutlier =FALSE,impute = TRUE);
  loginfo4 = list(record = paste0("\n using Data_impute to do knn imputation.\n ",
                                  "knn imputation will change by different random seed.\n",
                                  "sample abundance plot, density plot and pca plot is saved in plot dir to check the data distribtuion"),
                  input = if(is.na(opt$covarfile)) opt$quantfile else paste(opt$quantfile,"and",opt$covarfile),
                  output = paste0(opt$output,"_rmOut_impute.RData"),
                  packageVersion =paste("DDPNA",packageVersion("DDPNA")),
                  sessionInfo = sessionInfo(),
                  srcipt = "step1_MS V2.0", parameter = opt,
                  running_time = Sys.time())
  loginfo = c(loginfo, step1_imputation = list(loginfo4))
  #save quant info
  if(opt$value == "value") quant <- data_rmOut_imp$log2_value else quant = data_rmOut_imp$intensity
  if(opt$value == "value") quant_imp <- log2(data_rmOut_imp$intensity) else quant_imp = data_rmOut_imp$intensity
  IDinf <- data_rmOut_imp$inf
  otherobj <- otherobj[otherobj != "data_rmOut"]
  if("covar" %in% ls()) {
    if(all(colnames(quant) %in% rownames(covar))) {
      covar = data.frame(covar[rownames(covar) %in% colnames(quant),])
      if(ncol(covar) == 1) rownames(covar) <- colnames(quant);
      if(ncol(covar) == 1) colnames(covar) = colnames(covar_raw)
    } else {
      warning("covar is not the same with quant after outlier removal. Please check the ",opt$output ,"_rmOut_imp.RData")
      opt$step = 4
      message("Running is break in step 4.")
    }
    save(covar,quant,IDinf,data_rmOut_imp,loginfo,list = otherobj, file = paste0(opt$output,"_rmOut_impute.RData"))
  } else save(quant,IDinf,data_rmOut_imp,loginfo,list = otherobj, file = paste0(opt$output,"_rmOut_impute.RData"))
  cat('---missing value was imputed--- \n');
  step1p_final_rec = c(step1p_final_rec,
                         "*** missing value was imputed  ***\n")
}

### pca #####
if("covar" %in% ls()) {
  plotpdfname = "_step1"
  pca_list <- list(raw =quant_raw,
                   norm = if(grepl("1",opt$step)) quant_norm,
                   rm = if(grepl("2",opt$step)) quant_rm,
                   filter = if(grepl("3",opt$step))  quant_filter,
                   imp = if(grepl("4",opt$step)) quant_imp)
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
          if(i %in% c("norm","raw")) pca_plot_var(pcadata, covar_raw, title = i)
          if(i %in% c("rm","filter")) pca_plot_var(pcadata, covar_rm, title = i)
          if(i == "imp") pca_plot_var(pcadata, covar, title = i)
        } else
            warning(paste0(i, "step data run pca have error."))
      }
    }
    dev.off()
    message("pca plot was finished.")
  } else cat(" - covar is NULL, will not draw PCA plot.- \n")
}


### boxplot ######
cat("- Draw boxplot/density plot - \n")
# random pick 50 samples
if(ncol(quant) > 50) {message("Samples large than 50, random pick 50 sample to draw intensity distribution.")
  pos <- sample(1:ncol(quant),50)
  quant_raw2 <- quant_raw[,pos]
  if(grepl("1",opt$step)) quant_norm2 <- quant_norm[,pos]
  if(grepl("3",opt$step)) quant_filter2 <- quant_filter[,pos]
  if(grepl("4",opt$step)) quant_imp2 <- quant_imp[,pos]
  width = ceiling(50/5);
} else { 
  width = 10
  quant_raw2 <- quant_raw
  if(grepl("1",opt$step)) quant_norm2 <- quant_norm
  if(grepl("3",opt$step)) quant_filter2 <- quant_filter
  if(grepl("4",opt$step)) quant_imp2 <- quant_imp
}
pdf(paste0("plot/",opt$output,"_sample&gene distribution.pdf"),width = width,height = 6)
boxplot((quant_raw2), main = "raw intensity distribution in sample",ylab = "log2 value");
if(grepl("1",opt$step)) boxplot((quant_norm2),main = "norm intensity distribution in sample",ylab = "log2 value");
if(grepl("3",opt$step)) boxplot((quant_filter2),main = "gene filter intensity distribution in sample",ylab = "log2 value");
if(grepl("4",opt$step)) boxplot((quant_imp2), main = "impute data distribution in sample",ylab = "log2 value");
if(opt$value == "value") {
  try(densityplot(2^quant_raw,"raw data",savefile = FALSE,verbose = 1))
  if(grepl("1",opt$step)) try(densityplot(quant_norm,"norm data",savefile = FALSE))
  if(grepl("3",opt$step)) try(densityplot(2^quant_filter,"gene filter data",savefile = FALSE))
  if(grepl("4",opt$step)) try(densityplot(2^quant_imp,"impute data",savefile = FALSE))
} else {
  try(densityplot(quant_raw,"raw data",savefile = FALSE,verbose = 1))
  if(grepl("1",opt$step)) try(densityplot(quant_norm,"norm data",savefile = FALSE))
  if(grepl("3",opt$step)) try(densityplot(quant_filter,"gene filter data",savefile = FALSE))
  if(grepl("4",opt$step)) try(densityplot(quant_imp,"impute data",savefile = FALSE))
}
dev.off()

cat("\n")
cat(step1p_final_rec)
cat("--- step1_MS all finished ---\n")

