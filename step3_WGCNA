#Version: 2.0
#created: Kefu liu(liukefu19@163.com)
#Maintainer: Kefu Liu
#Date: July 30, 2025
#Software: R
#250720 add networkType and TOMType parameter

###0.package check ########
existpck <- rownames(installed.packages());
if(!"DDPNA" %in% existpck)
  stop("no DDPNA package, please install the R pakcage before run me.")
#if(!"Biobase" %in% existpck)
#  stop("no Biobase package, please install the R pakcage before run me.")
if(!"WGCNA" %in% existpck)
  stop("no WGCNA package, please install the R pakcage before run me.")

# object check##################
allobj <- ls();
if(!"input" %in% allobj) 
  input = NA;
if(!"output" %in% allobj) 
  output = NA;
if(!"funcdir" %in% allobj) 
  funcdir = NA;
if(!"onlytune" %in% allobj) onlytune = NA; 
if(!"tuning" %in% allobj) tuning = NA; 
if(!"power" %in% allobj) power = NA;
if(!"deepSplit" %in% allobj) deepSplit = NA;
if(!"minModuleSize" %in% allobj) minModuleSize = NA;
if(!"minKMEtoStay" %in% allobj) minKMEtoStay = NA;
if(!"detectCutHeight" %in% allobj) detectCutHeight = NA;
if(!"reassignThreshold" %in% allobj) reassignThreshold = NA;
if(!"mergeCutHeight" %in% allobj) mergeCutHeight = NA;
if(!"MaxMod0ratio" %in% allobj) MaxMod0ratio = NA;
if(!"networkType" %in% allobj) networkType = NA;
if(!"TOMType" %in% allobj) TOMType = NA;

opt = list(input= input,
           output = output, 
           funcdir = funcdir,
           onlytune = onlytune,
           power = power,
           tuning = tuning,
           deepSplit = deepSplit,
           minModuleSize = minModuleSize,
           minKMEtoStay = minKMEtoStay,
           detectCutHeight = detectCutHeight,
           reassignThreshold = reassignThreshold,
           mergeCutHeight = mergeCutHeight,
           MaxMod0ratio = MaxMod0ratio,
           networkType = networkType,
           TOMType = TOMType)
rmobj <- c(names(opt),"opt","option_list","allobj","allobj2","loginfo","existpck","funcpath")
allobj <- ls(); allobj <- allobj[!allobj %in% c("rmobj",rmobj)]

### 0. parameters evaluate ######
if(!dir.exists("plot")) dir.create("plot")
if(is.na(opt$input)) {
  if("quant" %in% ls()) {
    message("using quant obj to do the process.")
  } else if("wgcnadata" %in% ls()) {
    message("using wgcnadata obj to do the process.")
  } else
    stop ("no quantfile or quant object, see help info");
}
if(is.na(opt$output)) {
  opt$output = "data";
  message("no output, setting to 'data' as prefix")
}
if(is.na(opt$funcdir)) {
  warning('funcdir is setting to NA, means wgcna_plot.R in the default dir');
  if(!file.exists("wgcna_plot.R")) 
     stop("No wgcna_plot.R file in fundir path.")
} else {
  if(tail(unlist(strsplit(opt$funcdir,"")),1) == "/"){
    funcpath = paste0(opt$funcdir,"wgcna_plot.R");
  } else {
    opt$funcdir <-  paste0(opt$funcdir,"/")
    funcpath = paste0(opt$funcdir,"wgcna_plot.R");
  }
  if(!file.exists(funcpath)) 
      stop("No wgcna_plot.R file in fundir path.")
}
if(is.na(opt$tuning)) {
  opt$tuning = FALSE;
  message("tunning is setting to the default: FALSE.")}
if(is.na(opt$onlytune)) {
  opt$onlytune = FALSE;
  message("onlytune is setting to the default: FALSE.")}
if(is.na(opt$deepSplit)) {
  opt$deepSplit = 2;
  message("deepSplit is setting to the default: 2.")}
if(is.na(opt$minModuleSize)) {
  opt$minModuleSize = 50;
  message("minModuleSize is setting to the default: 50.")}
if(is.na(opt$minKMEtoStay)) {
  opt$minKMEtoStay = 0.3;
  message("minKMEtoStay is setting to the default: 0.3.")}
if(is.na(opt$detectCutHeight)) {
  opt$detectCutHeight = 0.995;
  message("detectCutHeight is setting to the default: 0.995.")}
if(is.na(opt$reassignThreshold)) {
  opt$reassignThreshold = 1e-4;
  message("reassignThreshold is setting to the default: 1e-4.")}
if(is.na(opt$mergeCutHeight)) {
  opt$mergeCutHeight = 0.15;
  message("mergeCutHeight is setting to the default: 0.15.")}
if(is.na(opt$MaxMod0ratio)) {
  opt$MaxMod0ratio = 0.5;
  message("MaxMod0ratio is setting to the default: 0.5.")}
if(is.na(opt$networkType)) {
  message("networkType is setting to the default: unsigned")}
if(is.na(opt$TOMType)) {
  opt$TOMType = "unsigned";
  message("TOMType is setting to the default: unsigned")}

### 0. input check ######
if(!is.na(opt$input)) {
  if(grepl("\\.RData$",opt$input)) {
    if(file.exists(opt$input)) 
      load(opt$input) else 
        stop(paste0("No ",opt$input," file."))
    name <- ls(); name <- name[grepl("quant",name)];
    if(length(name) == 1) {
      eval(parse(text = paste0("quant <- ", name)));
    }
  } else stop("input should be .RData format file.")
}
otherobj <- ls(); 
otherobj <- otherobj[!otherobj %in% c("rmobj",rmobj)]
otherobj <- otherobj[!sapply(otherobj, function(x) is.function(get(x)))]
if(!"wgcnadata" %in% ls()) wgcnadata = t(quant);

##auto_calculate_power#######
library(DDPNA)
library(WGCNA)

if(is.na(opt$power)) {
  message("power will calculated auto.")
  sftthreshould = 0.85
  SoftThresholdScaleGraph <- function(data,
                                      xlab = "Soft Threshold (power)",
                                      ylab = "Scale Free Topology Model Fit, signed R^2",
                                      main = "Scale independence",
                                      filename = NULL,...){#220604
    if (!requireNamespace("WGCNA", quietly = TRUE)) {
      stop ("WGCNA package is needed. Please install it.",
            call. = FALSE)}
    powers <- c(c(1:10), seq(from = 12, to = 24, by = 2));
    sft <- WGCNA::pickSoftThreshold(data,networkType = opt$networkType,powerVector =powers,...);
    if (is.character(filename) & length(filename) == 1){
      #if (!dir.exists("plot")) dir.create("plot");
      #pdf(paste0("plot/", "SoftThershold ", filename, ".pdf"))
      pdf(paste0( "SoftThershold ", filename, ".pdf")) #200703
    }
    plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[ ,2],
         xlab = xlab, ylab = ylab, type = "n",
         main = main);
    text(sft$fitIndices[ ,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
         labels=powers,cex=0.9,col="red");
    abline(h=0.90,col="red");
    if (is.character(filename) & length(filename) == 1) dev.off()
    sft;
  }
  sft = SoftThresholdScaleGraph(wgcnadata,filename = output,
                                corFnc = bicor)
  sft$powerEstimate = min(sft$fitIndices$Power[sft$fitIndices$SFT.R.sq > sftthreshould])
  if(is.null(sft$powerEstimate)) stop("No SFT reach the setting parameter.")
  if(is.infinite(sft$powerEstimate)) stop("No SFT reach the setting parameter.")
  cat("--the best WGCNA construct power is ", sft$powerEstimate,"--")
  opt$power = sft$powerEstimate
}

##tuning #####
if(opt$tuning) {
  message("tuning the WGCNA parameter.\n")
  #calculate cutheight
  {num <- nrow(wgcnadata)
    p = 0; r <- 0.6
    while (p < 0.05) {
      r <- r - 0.01
      p <- pt(-r * sqrt((num - 2)/(1 - r^2)), num - 2) * 2
    }
    while (p > 0.05) {
      r <- r + 0.001
      p <- pt(-r * sqrt((num - 2)/(1 - r^2)), num - 2) * 2
    }
    while (p < 0.05) {
      r <- r - 1e-04
      p <- pt(-r * sqrt((num - 2)/(1 - r^2)), num - 2) * 2
    }
    while (p > 0.05) {
      r <- r + 1e-05
      p <- pt(-r * sqrt((num - 2)/(1 - r^2)), num - 2) * 2
    }
    opt$detectCutHeight = 1 - r^opt$power
    opt$detectCutHeight = trunc(opt$detectCutHeight*1000000,5)/1000000
    if (opt$detectCutHeight < 0.995) 
      opt$detectCutHeight = 0.995
    rm(r, p, num)}
  
  WGCNAadjust <- try(wgcnatest(wgcnadata, power = opt$power,
                           corType = "bicor", 
                           networkType = opt$networkType,
                           TOMType = opt$TOMType,
                           maxBlockSize = 20000,
                           minModNum = 8, maxModNum = 30,
                           MaxMod0ratio = opt$MaxMod0ratio,
                           detectCutHeight = opt$detectCutHeight,
                           reassignThreshold=opt$reassignThreshold,
                           mergeCutHeight = opt$mergeCutHeight,
                           minModSize = opt$minModuleSize))
  if(!inherits(WGCNAadjust,"try-error")){
    write.csv(WGCNAadjust,file = paste0(opt$output,"_WGCNAadjust.csv"))
    if(ncol(WGCNAadjust) == 0) stop("No tuning to satisfy.")
    pos = which(rownames(WGCNAadjust) == 0)
    WGCNAadjust2 <- WGCNAadjust[,unlist(is.na(WGCNAadjust[pos+17,])) & unlist(!is.na(WGCNAadjust[10+pos,]))]
    if(ncol(WGCNAadjust2) == 0) {
      message("No tuning to make the ModNum between 10-16.")
      WGCNAadjust2 = WGCNAadjust} 
    posstay = which(rownames(WGCNAadjust2) == "minKMEtoStay")
    if(sum(unlist(WGCNAadjust2[posstay,] == 0.3))>0) {
      WGCNAadjust2 = WGCNAadjust2[which(WGCNAadjust2[posstay,] == 0.3)]
      message("pick minKMEtoStay is 0.3.")
      opt$minKMEtoStay = 0.3
    } else opt$minKMEtoStay = 0.2
    opt$deepSplit = WGCNAadjust2[1,1]
    record1 = "run tunning, and auto pick a parameter.\n "
  } else {
    record1 = "tuning is failed.\n "
    message(record1)}
}
if(!opt$tuning)  record1 = NULL
if(is.na(opt$minKMEtoStay)) opt$minKMEtoStay = 0.3
##WGCNA###################
if(!opt$onlytune) {
  record1 = paste0(record1,"use bicor, unsigned TOMtype, pamRespectsDendro is FALSE, other parameter is setted by manual or auto.\n ",
                   "power:",opt$power,", deepSplit:",opt$deepSplit,
                   ", minModuleSize:",opt$minModuleSize,
                   ", minKMEtoStay:",opt$minKMEtoStay,
                   ", detectCutHeight:",opt$detectCutHeight,
                   ", reassignThreshold:",opt$reassignThreshold,
                   ", mergeCutHeight:",opt$mergeCutHeight,"\n ")
  cat(record1)
  CoExpNet = blockwiseModules(wgcnadata, power = opt$power, 
                              maxBlockSize = 20000, corType = "bicor",
                              networkType = opt$networkType,
                              TOMType = opt$TOMType,
                              deepSplit = opt$deepSplit,
                              detectCutHeight =opt$detectCutHeight,
                              minModuleSize = opt$minModuleSize,
                              minKMEtoStay = opt$minKMEtoStay,
                              reassignThreshold = opt$reassignThreshold, 
                              mergeCutHeight = opt$mergeCutHeight,
                              numericLabels = TRUE, pamRespectsDendro = FALSE,
                              verbose = 3)
  
  if(!"inf" %in% ls()) inf = data.frame(GN =rownames(quant))
  moduleinf <- Module_inf(CoExpNet,inf,inftype = "none",IDname = "GN")
  MEinf <- ME_inf(CoExpNet$MEs,quant,intensity.type = "none")
  
  loginfo1 = list(record = paste0(record1,
                                  "do WGCNA and draw regular WGCNA plot.\n"),
                  input = opt$input,
                  output = paste0(opt$output,"_WGCNAnet.Rdata"),
                  sessionInfo = sessionInfo(),
                  srcipt = "step3_WGCNA_R v1.1", parameter = opt,
                  running_time = Sys.time())
  
  loginfo = c(loginfo, step3_WGCNA_R = list(loginfo1))
  
  save(CoExpNet, MEinf, moduleinf,loginfo,list=otherobj, file = paste0(opt$output,"_WGCNAnet.Rdata")) 
  
  rm(input)
  if(length(unique(moduleinf$moduleNum)) != 1) {
    output = opt$output
    source(funcpath)
  } else message("No module detect, No WGCNA plot.")
  print(table(moduleinf$moduleNum))
} else {
  message("onlytune, see SoftThershold_xx.pdf or xx_WGCNAadjust.csv to see the parameter.")
}


  
