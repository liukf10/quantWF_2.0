existpck <- rownames(installed.packages());

#放在opt之后,记录脚本运行前，env包含的变量
rmobj <- c(names(opt),"opt","option_list","allobj","allobj2","loginfo","existpck","funcpath")
allobj <- ls(); allobj <- allobj[!allobj %in% c("rmobj",rmobj)]

##  0.1  ####################
# output 固定设置
if(is.na(opt$output)) {
  opt$output = "data";
  message("no output, setting to 'data' as prefix")
}
#funcdir 固定设置 读取preprocess_func.R
if(is.na(opt$funcdir)) {
  warning('funcdir is setting to NA, check preprocess_func.R in the default dir');
  if(file.exists("preprocess_func.R")) 
    source("preprocess_func.R") else 
      stop("No preprocess_func.R file in default path.")
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

###检查读取quant, 检查读取covar 记录otherobj和loginfo###
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
    covar <- read.csv(opt$covarfile) else 
      stop(paste0("No ",opt$covarfile," file."))
}
otherobj <- ls(); 
otherobj <- otherobj[!otherobj %in% c("rmobj",rmobj)]
otherobj <- otherobj[!sapply(otherobj, function(x) is.function(get(x)))]
if(!"loginfo" %in% ls()) loginfo = list();

####covar 匹配到 quant  vector时设置为group列  batch_adj不用这个代码########
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

#只在step2需要
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
#keepvar   可以是vector 或者逗号分隔
if(length(opt$keepvar) == 1 & is.na(opt$keepvar)[1])
{opt$keepvar = NULL;
message("keepvar is setting to NULL.")
} else {
  opt$keepvar <- unlist(strsplit(opt$keepvar,","));
  if(any(!opt$keepvar %in% colnames(covar)))
    stop("wrong keepvar name: keepvar is not in covar.")
}





#存文件前##########
#check IDinf
if("IDinf" %in% ls()) {
  if(!"ID" %in% colnames(IDinf)) 
    message("input have IDinf but no ID col, not process IDinf.") else {
      pos <- match(rownames(quant),IDinf$ID)
      if(any(is.na(pos))) message("input have IDinf, but not all ID in IDinf match with quant, so IDinf have NA number.")
      IDinf <- IDinf[pos,]
      IDinf$ID <- rownames(quant)
    }
}
### pvca 必须有covar #############
covar <- covar_batch; plotpdfname = "_batch_adj"
pvca4plot <- list(batch_adj_before =quant_raw,
                  batch_adj_after = quant_mega)
if("lme4" %in% existpck) { #pvca
  if(!is.null(covar)) {
    cat("- Draw pvca plot - \n")
    pct_threshold <- 0.6
    if(ncol(covar) > 0) {
      pvcapath <- paste0(opt$funcdir,"pvca.R")
      filename=NA; exp_design <- covar
      if(file.exists(pvcapath)) {
        for(ii in names(pvca4plot)) {
          pvcadata <- pvca4plot[[ii]]
          source(pvcapath)
          pvca4plot[[ii]] <- pvcaObj
        }
        #save pdf
        pdf(paste0("plot/",opt$output,plotpdfname,"_pvca.pdf"),width = 12,height = 6)
        for(i in names(pvca4plot)) {
          randomEffectsMatrixWtAveProp <- pvca4plot[[i]]
          bp <- barplot(randomEffectsMatrixWtAveProp, xlab = "Effects",
                        ylab = "Weighted average proportion variance",
                        ylim= c(0,1.1),col = c("blue"), las=2,  
                        main=paste0("PVCA estimation (",i,")"))
          new_values = round(randomEffectsMatrixWtAveProp , 3)
          text(bp,randomEffectsMatrixWtAveProp,labels = new_values, pos=3, cex = 0.8)
        }
        dev.off()
        message("pvca plot was finished.") 
      } else warning("pvca.R is not exist. It will not draw pvca plot.")
    }
  }
} else  cat(" - no lme4 pcakge, not draw PVCA plot.- \n")

### pca 必须有covar #########
covar <- covar_batch; plotpdfname = "_batch_adj"
pca_list <- list(batch_adj_before =quant_raw,
                 batch_adj_after = quant_mega)
if(!is.null(covar)) {
  cat("- Draw PCA plot - \n")
  for(i in names(pca_list)) {
    if(!is.null(pca_list[[i]])) pca_list[[i]] <- pca_dat(pca_list[[i]])
  }
  message("PCA is calcluated.")
  save(covar,pca_list,file = paste0(opt$output,plotpdfname,"_pca.RData"))
  message("pca result was saved.")
  pdf(paste0("plot/",opt$output,plotpdfname,"_pcaplot.pdf"), width = 10, height = 8)
  for(i in names(pca_list)) {
    pcadata <- pca_list[[i]]
    if(!is.null(pcadata)) {
      if(inherits(pcadata,"prcomp"))
        pca_plot_var(pcadata, covar, title = i) else 
          warning(paste0(i, "step data run pca have error."))
    }
  }
  dev.off()
  message("pca plot was finished.") 
} else cat(" - covar is NULL, will not draw PCA plot.- \n")

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
#每个样本绘制密度图 标上颜色
density_plot(quant_raw, covar_batch, group=opt$batch_name, title = paste0(opt$output,"_batch_rm_before"))
density_plot(quant_mega, covar_batch, group=opt$batch_name, title = paste0(opt$output,"_batch_rm_after"))
#样本箱线图 （最多画50个）
boxplot(quant_1, main = "batch adj before data distribution in sample",ylab = "value");
boxplot(quant_2, main = "batch adj after data distribution in sample",ylab = "value");
#单个基因表达的密度图
try(densityplot(quant_raw,"batch adj before",savefile = FALSE))
try(densityplot(quant_mega,"batch adj after",savefile = FALSE))
dev.off()

cat("--- all finished ---\n")