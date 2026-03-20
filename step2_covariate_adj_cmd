#Version: 2.0
#created: Kefu liu(liukefu19@163.com)
#Maintainer: Kefu Liu
#Date: July 30, 2025
#Software: R
#250919 add PCA to find hidden covariate methods
#250919 add BIC, and evaluate_method parameter: pvca,pca,BIC or both

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
  make_option("--output", action="store", default=NA,
              type='character',help="output prefix name"),
  make_option("--funcdir", action="store", default=NA,
              type='character',help="preprocess function file path"),
  make_option("--keepvar", action="store", default=NA,
              type='character',help="keeped var column name, can be a vector, NULL means No var will keep. default is diagnose."),
  make_option("--vartype", action="store", default="value",
              type='character',help="variation type: value or factor, default is value."),
  make_option("--adjcov", action="store", default=NA,
              type='character',help="adjust covariants: obserbed or hidden, default is obs&hid."),
  make_option("--method", action="store", default=NA,
              type='character',help="estimate hidden factor method: sva, peer or pca, default is sva."),
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
  make_option("--tol", action="store", default= 1e-20,
              type='numeric',help="the tolerance for detecting linear dependencies, default is 1e-20.")
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
if(!is.null(opt$keepvar)){
  if(all(is.na(opt$keepvar))) {
    opt$keepvar = "diagnose"; message("Setting keepvar as default: diagnose.")}
}
opt$vartype <- match.arg(opt$vartype,c("value","factor"));
if(is.na(opt$adjcov)) {
  opt$adjcov = "obs&hid"; message("Setting adjcov as default: obs&hid.")}
opt$adjcov <- match.arg(opt$adjcov,c("obs&hid","obs","hid","obs+hid"));
if(is.na(opt$method)) {
  opt$method = "sva"; message("Setting method as default: sva.")}
opt$method <- match.arg(opt$method,c("sva","peer","pca"));
opt$svamethod <- match.arg(opt$svamethod,c("be","leek"));
if(opt$method == "peer"){
  if(!"peer" %in% existpck)
    stop("no peer package, peer only install in linux system.")
}
if(opt$method == "sva") {
  if(!"sva" %in% existpck)
    stop("no sva package, please install the R pakcage before run me.")
}
if(is.na(opt$nPCAchoose)) {
  opt$nPCAchoose = "BE"; message("Setting evaluate_method as default: BE.")}
opt$nPCAchoose <- match.arg(opt$nPCAchoose,c("BE","Elbow"));
if(is.na(opt$evaluate_method)) {
  opt$evaluate_method = "pvca&pca"; message("Setting evaluate_method as default: pvca&pca.")}
if(is.na(opt$BIC_cores)) {
  opt$BIC_cores = 30; message("Setting BIC_cores as default: 30.")}

# opt$nfactor will check in peer section
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

### 0.3 obj check (covar can be NULL) ###########
pct_threshold <- 0.6
tol = opt$tol;
.Machine$double.eps = 2.220446e-16

if(!"covar" %in% ls()) stop ("no covar object, please check the data.");
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
### 0.4 parameter check #######
#check keepvar
if(!is.null(covar) & !is.null(opt$keepvar)) {
  opt$keepvar <- unlist(strsplit(opt$keepvar,","));
  if(any(!opt$keepvar %in% colnames(covar))) {
    stop("keepvar is not in covar. if do not need keepvar, please setting to NULL.")
  }
  if(opt$vartype == "factor") {
    for(i in 1:length(opt$keepvar))
      covar[,opt$keepvar[i]] <- as.numeric(as.factor(covar[,opt$keepvar[i]]))
  }
  if(ncol(covar) != length(opt$keepvar)) {
    pos = match(opt$keepvar,colnames(covar));
    covar <- covar[,c(pos,(1:ncol(covar))[-pos])]
  }
} 
if(!"./plot" %in% list.dirs()) {dir.create("plot")}
### 1. obs cov adjust #########
cat('-----data preprocess 2 covariate adjust----- \n')
#quant_raw  quant1(after step1) covar1(after step1)
quant_raw <- quant; quant1 <- quant; covar1 = covar;
#no covar will setting to hid
if(is.null(covar)) {
  #250919 no covar, force to set to pca.
  if(opt$adjcov != "hid") {
    opt$adjcov = "hid";
    opt$keepvar = NULL
    if(opt$method == "sva") {
      opt$method = "pca";
      warning("covar is NULL. Force to Setting adjcov to hid and method to pca.")
    } else warning("covar is NULL. Force to Setting adjcov to hid.")
  }
}
#keep all obs, setting to hid.
if(!is.null(covar)) {
  if(opt$adjcov %in% c("obs","obs+hid") & ncol(covar) == length(opt$keepvar)) {
    message("covar info is the same with keepvar. No obs will be adjusted.")
    message("Setting adjcov to hid.")
    opt$adjcov = "hid"
  }
}
#check whether solve beta by remove one var.
if(!is.null(covar)) {
  varname <- colnames(covar)
  fmla <- as.formula(paste("~ ", paste(varname, collapse= "+")));
  X = model.matrix(fmla, data=covar);
  solveX <- try(solve(t(X) %*% X),silent = TRUE)
  if(inherits(beta,"try-error")) {
    if(!is.null(opt$keepvar)) {
      varname <- varname[!varname %in% opt$keepvar]
    }
    for(i in 1:length(varname)) {
      fmla <- as.formula(paste("~ ", paste(c(varname[-i],opt$keepvar), collapse= "+")));
      X = model.matrix(fmla, data=covar);
      solveX <- try(solve(t(X) %*% X),silent = TRUE)
      if(!inherits(beta,"try-error")) {
        message(varname[i],"remove will make the solve sucess.")
        warning("Run covariate adjust by fix covar.")
        covar <- covar[,!colnames(covar) %in% varname[i]]
        break}
    }
    if(inherits(beta,"try-error")) 
      stop("cannot solve beta, check the covar info.")
  }
}
#adj obs first, quant1 is obs adj res, covar1 is keepvar covar.
if(!is.null(covar)) {
  if(opt$adjcov %in% c("obs","obs+hid") & ncol(covar) > length(opt$keepvar)) {
    fmla <- as.formula(paste("~ ", paste(colnames(covar), collapse= "+")));
    X = model.matrix(fmla, data=covar); #lm model(modified)
    Y = quant;
    beta = try((solve(t(X)%*%X)%*%t(X))%*%t(Y))
    if(inherits(beta,"try-error")) {
      stop(paste0("cannot solve obs beta, check the covar info."))
    }
    #no keepvar, only keep the intercept
    if(is.null(opt$keepvar))  {
      covar1 = NULL;
      if(ncol(X) == 2) 
        to_regress <- as.matrix(X[,-1]) %*% matrix(beta[-1,],nrow = 1) else
          to_regress = as.matrix(X[,-1]) %*% (as.matrix(beta[-1,]))
    }
    #have keepvar, keep the intercept and keepvar
    if(!is.null(opt$keepvar)) {
      covar1 = covar[opt$keepvar];
      keepN = NULL;
      for(i in 1:length(opt$keepvar)) keepN = c(keepN,grep(opt$keepvar[i],colnames(X)))
      if(ncol(X)-length(keepN) == 1) 
        to_regress = 0 else if(ncol(X)-length(keepN) == 2) 
          to_regress <- as.matrix(X[,c(-1,-keepN)]) %*% matrix(beta[c(-1,-keepN),],nrow = 1) else
            to_regress = as.matrix(X[,c(-1,-keepN)]) %*% (as.matrix(beta[c(-1,-keepN),])) 
    }
    quant_adj_obs <- quant - t(to_regress);
    
    record = if(is.null(opt$keepvar)) paste0("all covar are adjusted.\n ") else paste0("covar ",paste0(opt$keepvar,collapse = ",")," is keeped.\n ")
    loginfo1 = list(record = paste0("adjust by linear regression, ",record,
                                    "covar correlation is plotted. \n ",
                                    "pvca before and after obs cov adjust and pca after obs adjust is plotted. \n "),
                    input = if(is.na(opt$covarfile)) opt$quantfile else paste(opt$quantfile,"and",opt$covarfile),
                    output = paste0(opt$output,"_adj_obs.RData"),
                    sessionInfo = sessionInfo(),
                    srcipt = "step2_covariate_adj v2.0", parameter = opt,
                    running_time = Sys.time())
    loginfo = c(loginfo, step2_cov_adj_obs = list(loginfo1))
    #obs save it as quant_final else quant_adj_obs
    if(opt$adjcov == "obs") { 
      quant = quant_adj_obs; 
      loginfo$step2_cov_adj_obs$output <- paste0(opt$output,"_adj_final.RData")
      save(quant, covar,loginfo,list = otherobj, file = paste0(opt$output,"_adj_final.RData")) 
    } else save(quant,covar,loginfo,list = otherobj,file = paste0(opt$output,"_adj_obs.RData")) 
    quant1 <- quant_adj_obs;
  }
}

### 2. predict hidden covariate and adjust ######
if(opt$adjcov != "obs" ) {
  #must have keepvar
  if(opt$method == "sva") {
    if(is.null(opt$keepvar)) stop("No keepvar, cannot use sva hidden factor estimate.")
    suppressPackageStartupMessages(library(sva))
    #hid: mod is keep covar; mod0 is 1
    #obs+hid: mod is keep covar; mod0 is 1
    if(opt$adjcov %in% c("obs+hid","hid")) {
      fmla <- as.formula(paste("~ ", paste(opt$keepvar, collapse= "+")));
      mod = model.matrix(fmla, data=covar1);
      mod0 = model.matrix(~1,data=covar1);
    }
    #obs&hid: mod is all covar; mod0 is keepvar
    if(opt$adjcov == "obs&hid") {
      fmla <- as.formula(paste("~ ", paste(colnames(covar1), collapse= "+")));
      mod = model.matrix(fmla,data=covar1);
      fmla <- as.formula(paste("~ ", paste(colnames(covar1)[!colnames(covar1) %in% opt$keepvar], collapse= "+")));
      mod0 = model.matrix(fmla, data=covar1);
    }
    n.sv = try(num.sv(quant1, mod, method = opt$svamethod));
    if(inherits(n.sv,"try-error")) {
      stop(paste0("cannot calculate num.sv, stop the process."))
    }
    if(n.sv > 0) {
      cat(paste0("--- sva: ", opt$svamethod ," method found ",n.sv," hidden covariance ---\n")) 
      svobj = try(sva(as.matrix(quant1),mod,mod0,n.sv=n.sv));
      if(inherits(svobj,"try-error")) {
        stop(paste0("sva fail, stop the process, please check covar or quant."))
      }} else cat(paste0("--- sva: ", opt$svamethod ," method found zero hidden covariance ---\n")) 
    #covar_hid
    if(opt$adjcov %in% c("obs+hid","hid")) {
      if(n.sv == 0) covar_hid = covar1[opt$keepvar] else
        covar_hid = data.frame(svobj$sv,covar1[opt$keepvar])
    } else {
      if(n.sv == 0) covar_hid = covar1 else {
        covar_hid = data.frame(svobj$sv, covar1);
      }
    }
  }
  #usual used for QTL analysis
  if(opt$method == "peer") {
    if(is.na(opt$nFactor))  {
      if(ncol(quant1) < 150) 
        nFactor = 15 else if (ncol(quant1) < 250)
          nFactor = 30 else if (ncol(quant1) < 350)
            nFactor = 45 else nFactor = 60
            message("Setting nFactor depend on sample size (GTEx suggestion): ", nFactor)
    } else nFactor = opt$nFactor
    if(nFactor < 1) stop("wrong nFactor.")
    suppressPackageStartupMessages(library(peer))
    expr = t(as.matrix(quant1)) 
    model = PEER()  # create the model object
    invisible(PEER_setPhenoMean(model,expr))# set the observed data
    invisible(PEER_setNk(model,nFactor)) # gradient number of factors
    invisible(PEER_getNk(model))
    invisible(PEER_setAdd_mean(model, TRUE))  # include an additional factor (covariate) to account for the mean expression
    #covar1 = NULL 1.obs+hid and no keep; 2.covar is NULL
    if(!is.null(covar1)) {
      invisible(PEER_setCovariates(model, as.matrix(covar1)))  # adding covariates has no effect on the model?
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
    if(is.null(covar1))
      covar_hid <- data.frame(factors) else
        covar_hid <- data.frame(factors,covar1)
  }
  #usual used for QTL analysis
  if(opt$method == "pca"){
    expr = t(as.matrix(quant1)) 
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
      covar_hid <-data.frame(PCsTop,covar[,colnames(covar) %in% colnames(knowncovFiltered)])
    } else covar_hid <- scale(PCsTop)
  }
  fmla <- as.formula(paste("~ ", paste(colnames(covar_hid), collapse= "+")));
  X = model.matrix(fmla, data=covar_hid)
  Y = quant1
  beta = try((solve(t(X)%*%X, tol = tol)%*%t(X))%*%t(Y))
  record2 = if(opt$method == "sva") 
    paste("using sva to find the hidden covar and find",n.sv,"hidden covariance.\n ") else 
      paste("using peer to calculate", opt$nFactor,"hidden covariance info.\n ")
  
  if(inherits(beta,"try-error")) {
    loginfoerror = list(record = paste0(record2, "adjust by linear regression, but cannot solve the beta, may tolerance is too small.\n "),
                        input = if(is.na(opt$covarfile)) opt$quantfile else paste(opt$quantfile,"and",opt$covarfile),
                        output = "only save the covar_hid",
                        packageVersion =if(opt$method == "sva") paste("sva",packageVersion("sva")) else paste("peer",packageVersion("peer")),
                        sessionInfo = sessionInfo(),
                        srcipt = "step2_covariate_adj v2.0", parameter = opt,
                        running_time = Sys.time())
    loginfo = c(loginfo, step2_cov_adj_error = list(loginfoerror))
    save(covar_hid,loginfo, file = paste0("covar_hid_",opt$output,".RData"))
    stop(paste0("cannot solve beta, just save hidden factor in covar_hid_",opt$output,".RData"))
  }
  if(!is.null(opt$keepvar)) {
    keepN = NULL;
    for(i in 1:length(opt$keepvar)) keepN = c(keepN,grep(opt$keepvar[i],colnames(X)))
    if(ncol(X)-length(keepN) == 1) 
      to_regress = 0 else if(ncol(X)-length(keepN) == 2)
        to_regress <- as.matrix(X[,c(-1,-keepN)]) %*% matrix(beta[c(-1,-keepN),],nrow = 1) else
          to_regress = as.matrix(X[,c(-1,-keepN)]) %*% (as.matrix(beta[c(-1,-keepN),]))
  }
  if(is.null(opt$keepvar)) {
    if (ncol(X) == 2) to_regress <- as.matrix(X[,-1]) %*% matrix(beta[-1,],nrow = 1) else
      to_regress = as.matrix(X[,-1]) %*% (as.matrix(beta[-1,]))
  }
  
  quant <- quant1 - t(to_regress);
  
  record = if(is.null(opt$keepvar)) paste0("all covar are adjusted.\n ") else paste0("covar ",paste0(opt$keepvar,collapse = ",")," is keeped.\n ")
  loginfo2 = list(record = paste0(record2, "adjust by linear regression, ",record,
                                  "covar included hidden covariance correlation is plotted. \n ",
                                  "pvca before and after obs cov adjust and pca after obs adjust is plotted. \n "),
                  input = if(is.na(opt$covarfile)) opt$quantfile else paste(opt$quantfile,"and",opt$covarfile),
                  output = paste0(opt$output,"_adj_final.RData"),
                  packageVersion =if(opt$method == "sva") paste("sva",packageVersion("sva")) else if(opt$method == "peer") paste("peer",packageVersion("peer")),
                  sessionInfo = sessionInfo(),
                  srcipt = "step2_covariate_adj v2.0", parameter = opt,
                  running_time = Sys.time())
  loginfo = c(loginfo, step2_cov_adj = list(loginfo2))
  
  save(quant, covar_hid, covar,loginfo,list = otherobj,
       file = paste0(opt$output,"_adj_final.RData")) #230201 updated
}

###  covariants correlation #######
if(!is.null(covar)){
  if(ncol(covar) != 1) {
    covar2 <- covar;
    for (i in 1:ncol(covar2)) {
      if(!class(covar2[,i]) %in% c("numeric","integer"))
        covar2[,i] <- as.numeric(as.factor(covar2[,i]))
    }
  } }
if(!grepl("hid",opt$adjcov)) covar_hid = NULL
if(!is.null(covar_hid)){
  if(ncol(covar_hid) != 1) {
    covar_hid2 <- covar_hid;
    for (i in 1:ncol(covar_hid2)) {
      if(!class(covar_hid2[,i]) %in% c("numeric","integer"))
        covar_hid2[,i] <- as.numeric(as.factor(covar_hid2[,i]))
    }
  }
}
if(!is.null(covar) | !is.null(covar_hid)) {
  cat("- Draw covar correlation plot - \n")
  #plot 
  suppressPackageStartupMessages(require(corrplot))
  col2 <- colorRampPalette(c("#053061","#2166AC","#4393C3","#92C5DE",
                             "#D1E5F0", "#FFFFFF","#FDDBC7","#F4A582",
                             "#D6604D","#B2182B","#67001F"))
  pdf(file = paste0("plot/",opt$output,"_",opt$method,"_cov_cor.pdf"),
      width = 10, height = 10)
  if(!is.null(covar))
    if(ncol(covar) != 1)
      corrplot(cor(covar2),p.mat = cor.mtest(covar2)$p, col = col2(200),
               insig = "p-value", order = "original",method = "ellipse",
               type = "lower", tl.pos = "d")
  if(!is.null(covar_hid))
    if(ncol(covar_hid) != 1)
      corrplot(cor(covar_hid2),p.mat = cor.mtest(covar_hid2)$p, col = col2(200),
               insig = "p-value", order = "original",method = "ellipse",
               type = "lower", tl.pos = "d")
  dev.off()
}
###  pvca  ########
if(grepl("pvca|PVCA",opt$evaluate_method)) {
  plotpdfname = paste0("-adj_",opt$adjcov)
  pvcalist <- list(adjust_before =quant_raw,
                   adjust_observed = quant1,
                   adjust_hidden = quant)
  if(opt$adjcov == "obs") pvcalist <- pvcalist[-3]
  if(opt$adjcov == "obs&hid") pvcalist <- pvcalist[-2] 
  if("lme4" %in% existpck) { #pvca
    if(!is.null(covar)) {
      cat("- Draw pvca plot - \n")
      if(ncol(covar) > 0) {
        pdf(paste0("plot/",opt$output,plotpdfname,"_pvca.pdf"),width = 12,height = 6)
        for(ii in names(pvcalist)) {
          suppressMessages(pvcaresult <- pvca_evaluate(pvcalist[[ii]],covar, pct_threshold = 0.6,output = ii))
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

###  pca #########
if(grepl("pca|PCA",opt$evaluate_method)) {
  plotpdfname = paste0("-adj_",opt$adjcov)
  pca_list <- list(adjust_before =quant_raw,
                   adjust_observed = quant1,
                   adjust_hidden = quant)
  if(opt$adjcov == "obs") pca_list <- pca_list[-3]
  if(opt$adjcov == "obs&hid") pca_list <- pca_list[-2]
  if(!is.null(covar)) {
    cat("- Draw PCA plot - \n")
    for(i in names(pca_list)) {
      if(opt$method == "pca" & i == "adj_before") {
        prcompResult$rotation <- prcompResult$rotation[,1:max(ncol(PCsTop),20)]
        prcompResult$x <- prcompResult$x[,1:max(ncol(PCsTop),20)]
        pca_list[[i]] <- prcompResult
      } else pca_list[[i]] <- pca_dat(pca_list[[i]])
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
} 

### boxplot ############
cat("- Draw boxplot/density plot - \n")
# random pick 50 samples to draw boxplot, rmOutlier will not darw
if(ncol(quant) > 50) {message("Samples large than 50, random pick 50 sample to draw intensity distribution.")
  pos <- sample(1:ncol(quant),50)
  quant_1 <- quant_raw[,pos]
  if(opt$adjcov %in% c("obs","obs+hid"))  quant1_1<- quant1[,pos]
  if(grepl("hid",opt$adjcov)) quant2_1 <- quant[,pos]
  width = ceiling(50/5);
} else {
  width = 10
  quant_1 <- quant_raw
  if(opt$adjcov %in% c("obs","obs+hid"))  quant1_1<- quant1
  if(grepl("hid",opt$adjcov)) quant2_1 <- quant
}
pdf(paste0("plot/",opt$output,"_adj_sample&gene distribution.pdf"),width = width,height = 6)
boxplot(quant_1, main = "adj before data distribution in sample",ylab = "log2 value");
if(opt$adjcov %in% c("obs","obs+hid")) boxplot(quant1_1,main = "adjust observed data distribution in sample",ylab = "log2 value");
if(grepl("hid",opt$adjcov)) boxplot(quant2_1, main = "adjust final data distribution in sample",ylab = "log2 value");
try(densityplot(quant_raw,"adj before",savefile = FALSE))
if(opt$adjcov %in% c("obs","obs+hid")) try(densityplot(quant1,"adj observed",savefile = FALSE))
if(grepl("hid",opt$adjcov)) try(densityplot(quant,"adj final",savefile = FALSE))
dev.off()

### BIC #######
if(grepl("bic|BIC",opt$evaluate_method)) {
  plotpdfname = paste0("-adj_",opt$adjcov)
  BIClist <- list(adjust_before =quant_raw,
                  adjust_observed = quant1,
                  adjust_hidden = quant)
  if(opt$adjcov == "obs") BIClist <- BIClist[-3]
  if(opt$adjcov == "obs&hid") BIClist <- BIClist[-2] 
  
  cat("- Draw BIC plot - \n")
  pdf(paste0("plot/",opt$output,plotpdfname,"_bic.pdf"),width = 12,height = 6)
  for(i in names(BIClist)) {
    BICresult <- BIC_evaluate(BIClist[[i]],covar,num_cores = opt$BIC_cores,output = paste0("BIC estimation (",i,")"))
    BIClist[[i]] <- BICresult$BICres
    print(BICresult$pic)
  }
  dev.off()
  save(BIClist,file = paste0(opt$output,plotpdfname,"_BIC.RData") )
  message("BIC plot was finished.") 
  
}
### finish report ####
cat("--- all finished ---\n")
cat("----step2_covariate_adj finished ------\n\n")