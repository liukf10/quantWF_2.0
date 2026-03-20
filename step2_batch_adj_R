#Version: 2.0
#created: Kefu liu(liukefu19@163.com)
#Maintainer: Kefu Liu
#Date: July 30, 2025
#Software: R

#covarfile must be csv file 
#quantfile must be RData or csv file.

### package check ########
existpck <- rownames(installed.packages());



### object check ##################
allobj <- ls();
if(!"quantfile" %in% allobj)
  quantfile = NA;
if(!"covarfile" %in% allobj) 
  covarfile = NA;
if(!"output" %in% allobj) 
  output = NA;
if(!"funcdir" %in% allobj)  funcdir = NA;
if(!"keepvar" %in% allobj)  keepvar = NA
if(!"batch_name" %in% allobj) batch_name = NA;
if(!"ref_batch" %in% allobj) ref_batch = NA;
if(!"method" %in% allobj)  method = NA; 
if(!"evaluate_method" %in% allobj) {evaluate_method = NA;}
if(!"BIC_cores" %in% allobj) {BIC_cores = NA;}
opt = list(quantfile= quantfile,
           covarfile = covarfile, 
           output = output, 
           funcdir = funcdir,
           keepvar = keepvar,
           batch_name = batch_name,
           ref_batch=ref_batch,
           method = method,
           evaluate_method =evaluate_method,
           BIC_cores = BIC_cores)
rmobj <- c(names(opt),"opt","option_list","allobj","allobj2","loginfo","existpck","funcpath")
allobj <- ls(); allobj <- allobj[!allobj %in% c("rmobj",rmobj)]

### 0.1 parameters evaluate (batch) ######
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
if(is.na(opt$batch_name)) {
  opt$batch_name = "batch"; message("Setting batch_name as default: batch.")}
if(is.na(opt$ref_batch)) {
  opt$ref_batch = NULL; message("Setting ref_batch as default: NULL.")}
if(is.na(opt$method)) {
  opt$method = "combat"; message("Setting method as default: combat.")}
opt$method <- match.arg(opt$method,c("combat","limma","combat_seq"));
if(grepl("combat",opt$method)) {
  if(!"sva" %in% existpck)
    stop("no sva package, please install the R pakcage before run me.")
  suppressPackageStartupMessages(library(sva))
}
if(opt$method == "limma") {
  if(!"limma" %in% existpck)
    stop("no limma package, please install the R pakcage before run me.")
  suppressPackageStartupMessages(library(limma))
}
if(is.na(opt$evaluate_method)) {
  opt$evaluate_method = "pvca&pca"; message("Setting evaluate_method as default: pvca&pca.")}
if(is.na(opt$BIC_cores)) {
  opt$BIC_cores = 30; message("Setting BIC_cores as default: 30.")}
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
### 0.3 obj check (covar should have batchname) ###########
if(length(opt$batch_name) > 1) stop("wrong batch_name: only one batch allowed.") #230202 updated
if(!"covar" %in% ls()) stop ("no covar object, please check the data.");
# covar is vector, put batchname into it.
if("covar" %in% ls()) {
  if(!is.null(covar)) {
    #when covar is a vector will not match and transfer to data.frame 
    if(is.vector(covar)) {
      if(length(covar) == ncol(quant)) {
        covar <- data.frame(covar);
        if(is.null(opt$batch_name)) stop("No batch_name.")
        colnames(covar) = opt$batch_name 
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
if(!opt$batch_name %in% colnames(covar)) 
  stop(paste0("wrong batch_name: ",opt$batch_name," not exist in covar.")) #230202 updated
#keepvar checking (default is setting to NULL )
if(length(opt$keepvar) == 1 & is.na(opt$keepvar)[1])
{opt$keepvar = NULL;
message("keepvar is setting to NULL.")
} else {
  opt$keepvar <- unlist(strsplit(opt$keepvar,","));
  if(any(!opt$keepvar %in% colnames(covar)))
    stop("wrong keepvar name: keepvar is not in covar.")
  }
if(!is.null(opt$ref_batch)) {
  message("only use ref_batch in Combat.")
  if(length(opt$ref_batch) != 1) stop("ref_batch should be one vector.")
  if(!opt$ref_batch %in% covar[,opt$batch_name])
    stop("ref_batch should be one of vector in batch_name column.")
}
cat('-----data preprocess 2 batch effect adjust----- \n')
if(!"./plot" %in% list.dirs()) {dir.create("plot")}
### 1. batch adjust #########
if(!is.null(opt$keepvar)) {
  fmla <- as.formula(paste("~ ", paste(opt$keepvar, collapse= "+")));
  mod = model.matrix(fmla, data=covar)
} else mod = NULL
batch = covar[,opt$batch_name]
covar_batch <- covar
#remove batch column
covar <- covar[!colnames(covar) %in% opt$batch_name]

if(opt$method == "combat") {
  suppressMessages(quant_mega <- ComBat(as.matrix(quant), batch=batch, mod=mod, prior.plots = F, ref.batch = opt$ref_batch))
  quant_mega = as.data.frame(quant_mega)
} else if(opt$method == "combat_seq") {
  invisible(capture.output(quant_mega <- ComBat_seq(as.matrix(quant), batch=batch, covar_mod=mod)))
  quant_mega = as.data.frame(quant_mega)
} else {
  if(is.null(mod)) mod <- matrix(1,ncol(quant))
  quant_mega <- removeBatchEffect(as.matrix(quant), batch = batch, design = mod)
  quant_mega = as.data.frame(quant_mega)
}

### 2. save #####
record = if(grepl("combat",opt$method)) 
  "using ComBat to remove batch effect.\n " else 
    "using removeBatchEffect to remove batch effect.\n "
loginfo1 = list(record = paste0(record, 
                                "draw pca, density, and MDS plot before and after batch effect removal in plot dir.\n "),
                input = if(is.na(opt$covarfile)) opt$quantfile else paste(opt$quantfile,"and",opt$covarfile),
                output = paste0(opt$output,"_batch_remove.RData"),
                packageVersion =if(grepl("combat",opt$method)) paste("sva",packageVersion("sva")) else paste("limma",packageVersion("limma")),
                sessionInfo = sessionInfo(),
                srcipt = "step2_batch_adj v2.0", parameter = opt,
                running_time = Sys.time())
loginfo = c(loginfo, step2_batch_adj = list(loginfo1))
quant_raw <- quant
quant <- quant_mega
save(quant,covar,loginfo,covar_batch, list = otherobj, file = paste0(opt$output,"_batch_remove.RData"))

###  pvca  ########
if(grepl("pvca|PVCA",opt$evaluate_method)) {
plotpdfname = "_batch_adj"
pvcalist <- list(batch_adj_before =quant_raw,
               batch_adj_after = quant_mega)
if("lme4" %in% existpck) { #pvca
  if(!is.null(covar_batch)) {
    cat("- Draw pvca plot - \n")
    if(ncol(covar_batch) > 0) {
      pdf(paste0("plot/",opt$output,plotpdfname,"_pvca.pdf"),width = 12,height = 6)
      for(ii in names(pvcalist)) {
        suppressMessages(pvcaresult <- pvca_evaluate(pvcalist[[ii]],covar_batch, pct_threshold = 0.6,output = ii))
        pvcalist[[ii]] <- pvcaresult$pvca_res
        print(pvcaresult$pic)
      }
      dev.off()
      save(pvcalist,file = paste0(opt$output,plotpdfname,"_PVCA.RData") )
      message("pvca plot was finished.") 
    }
  }
} else  cat(" - no lme4 pcakge, not draw PVCA plot.- \n")
}
### MDS ######
if(grepl("MDS|mds",opt$evaluate_method)) {
plotpdfname = "_batch_adj"
MDS_list <- list(batch_adj_before =quant_raw,
                  batch_adj_after = quant_mega)
if(!is.null(covar_batch)) {
  cat("- Draw MDS plot - \n")
  for(i in names(MDS_list)) {
    MDS_list[[i]] <- MDS_dat(MDS_list[[i]])
  }
  save(covar_batch,MDS_list,file = paste0(opt$output,plotpdfname,"_MDS.RData"))
  message("MDS result was saved.")
  pdf(paste0("plot/",opt$output,plotpdfname,"_MDSplot.pdf"), width = 10, height = 8)
  for(i in names(MDS_list)) {
    mdsdata <- MDS_list[[i]]
    if(!is.null(mdsdata)) {
      if(inherits(mdsdata,"list"))
        MDS_plot(mdsdata, covar_batch, group=opt$batch_name,title = i) else 
          warning(paste0(i, "step data run MDS have error."))
    }
  }
  dev.off()
  message("MDS plot was finished.") 
} else cat(" - covar is NULL, will not draw MDS plot.- \n")
}
##  pca #########
if(grepl("pca|PCA",opt$evaluate_method)) {
plotpdfname = "_batch_adj"
pca_list <- list(batch_adj_before =quant_raw,
                 batch_adj_after = quant_mega)
if(!is.null(covar)) {
  cat("- Draw PCA plot - \n")
  for(i in names(pca_list)) {
    pca_list[[i]] <- pca_dat(pca_list[[i]])
  }
  message("PCA is calcluated.")
  save(covar_batch,pca_list,file = paste0(opt$output,plotpdfname,"_pca.RData"))
  message("pca result was saved.")
  pdf(paste0("plot/",opt$output,plotpdfname,"_pcaplot.pdf"), width = 10, height = 8)
  for(i in names(pca_list)) {
    pcadata <- pca_list[[i]]
    if(!is.null(pcadata)) {
      if(inherits(pcadata,"prcomp"))
        pca_plot_var(pcadata, covar_batch, title = i) else 
          warning(paste0(i, "step data run pca have error."))
    }
  }
  dev.off()
  message("pca plot was finished.") 
} else cat(" - covar is NULL, will not draw PCA plot.- \n")
}
### boxplot&densityplot ############
cat("- Draw boxplot/density plot - \n")
# random pick 50 samples to draw boxplot, rmOutlier will not darw
if(ncol(quant) > 50) {message("Samples large than 50, random pick 50 sample to draw intensity distribution.")
  pos <- sample(1:ncol(quant),50)
  quant_1 <- quant_raw[,pos]
  quant_2 <- quant_mega[,pos]
  width = ceiling(50/5);
} else {
  width = 10
  quant_1 <- quant_raw
  quant_2 <- quant_mega
}
pdf(paste0("plot/",opt$output,"_batch_adj_sample&gene distribution.pdf"),width = width,height = 6)
if(opt$method == "combat_seq") 
  density_plot(log2(quant_raw+1), covar_batch, group=opt$batch_name, title = paste0(opt$output,"_batch_rm_before")) else 
    density_plot(quant_raw, covar_batch, group=opt$batch_name, title = paste0(opt$output,"_batch_rm_before"))
if(opt$method == "combat_seq")
  density_plot(log2(quant_mega+1), covar_batch, group=opt$batch_name, title = paste0(opt$output,"_batch_rm_after")) else
    density_plot(quant_mega, covar_batch, group=opt$batch_name, title = paste0(opt$output,"_batch_rm_after"))
boxplot(quant_1, main = "batch adj before data distribution in sample",ylab = "value");
boxplot(quant_2, main = "batch adj after data distribution in sample",ylab = "value");
try(densityplot(quant_raw,"batch adj before",savefile = FALSE))
try(densityplot(quant_mega,"batch adj after",savefile = FALSE))
dev.off()

### BIC #######
if(grepl("bic|BIC",opt$evaluate_method)) {
plotpdfname = "_batch_adj"
BIClist <- list(batch_adj_before =quant_raw,
                 batch_adj_after = quant_mega)
cat("- Draw BIC plot - \n")
pdf(paste0("plot/",opt$output,plotpdfname,"_bic.pdf"),width = 12,height = 6)
for(i in names(BIClist)) {
  BICresult <- BIC_evaluate(BIClist[[i]],covar_batch,num_cores = opt$BIC_cores,output = paste0("BIC estimation (",i,")"))
  BIClist[[i]] <- BICresult$BICres
  print(BICresult$pic)
  }
dev.off()
save(BIClist,file = paste0(opt$output,plotpdfname,"_BIC.RData") )
message("BIC plot was finished.") 
}


cat("--- all finished ---\n")
