#Version: 2.0
#created: Kefu liu(liukefu19@163.com)
#Maintainer: Kefu Liu
#Date: July 30, 2025
#Software: R

#log2FCcut = 1; adjPcut = 0.05; 
#varequal=TRUE; do_ttest = FALSE; small_size = TRUE
#quant; covar; grpname = "group" ;   groupB groupA
# object check ##################
allobj <- ls(); 

if(!"input" %in% allobj) 
  input = NA;
if(!"output" %in% allobj) 
  output = NA;
if(!"grpname" %in% allobj) 
  grpname = NA;
if(!"groupA" %in% allobj) 
  groupA = NA;
if(!"groupB" %in% allobj) 
  groupB = NA;
if(!"do_ttest" %in% allobj) 
  do_ttest = NA;
if(!"varequal" %in% allobj) 
  varequal = NA;
if(!"log2FCcut" %in% allobj) 
  log2FCcut = NA;
if(!"adjPcut" %in% allobj) 
  adjPcut = NA;

opt = list(input = input, output = output,
           grpname = grpname, groupA= groupA,groupB = groupB,
           do_ttest=do_ttest,varequal=varequal,
           log2FCcut=log2FCcut,adjPcut=adjPcut)

rmobj <- c(names(opt),"opt","option_list","allobj","allobj2","quant","loginfo","existpck","funcpath")
allobj <- ls(); allobj <- allobj[!allobj %in% c("rmobj",rmobj)]

### 0 parameters evaluate&read input ######
if(is.na(opt$output)) {
  opt$output = "data";
  message("no output, setting to 'data' as prefix")
}
if(is.na(opt$grpname)) {
  opt$grpname = "group"; message("Setting grpname as default: group.")}
if(is.na(opt$log2FCcut)) {
  opt$log2FCcut = 1; message("Setting log2FCcut as default: 1.")}
if(is.na(opt$adjPcut)) {
  opt$adjPcut = 0.05; message("Setting adjPcut as default: 0.05.")}
if(is.na(opt$varequal)) {
  opt$varequal = TRUE; message("Setting varequal as default: TRUE.")}

if(!is.na(opt$input)) {
  if(grepl("\\.RData$",opt$input)) {
    if(file.exists(opt$input)) 
      load(opt$input) else 
        stop(paste0("No ",opt$input," file."))
  } else stop("input should be .RData format file.")
}
if(!"quant" %in% ls()) 
  stop ("no quant object, see help info.");
if(!"covar" %in% ls()) stop("No covar, see help info.");
otherobj <- ls(); 
otherobj <- otherobj[!otherobj %in% c("rmobj",rmobj)]
otherobj <- otherobj[!sapply(otherobj, function(x) is.function(get(x)))]
if(!"loginfo" %in% ls()) loginfo = list();

### covar check #############
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
if(any(!opt$grpname %in% colnames(covar)))
  stop("grpname is not in covar, please check it.")
group_list <- covar[,opt$grpname];
group = unique(group_list)

#groupA groupB
if(length(group) != 2) {
  if(all(is.na(opt$groupA))) stop("No groupA, see readme document.")
  if(all(is.na(opt$groupB))) stop("No groupB, see readme document.")
  opt$groupA = unlist(strsplit(opt$groupA,","));
  opt$groupB = unlist(strsplit(opt$groupB,",")); 
  if(!all(opt$groupA %in% group)) stop("Not all groupA belongs to group.")
  if(!all(opt$groupB %in% group)) stop("Not all groupB belongs to group.")
  if(length(opt$groupA) != length(opt$groupB)) stop("groupA and groupB should have the same length.")
} else {
  if(all(is.na(opt$groupA)) | all(is.na(opt$groupA))) {
    opt$groupA = group[1]; opt$groupB = group[2];
    message(paste0("setting groupA as ",group[1],"; groupB as ",group[2]))
  }
  opt$groupA = unlist(strsplit(opt$groupA,","));
  opt$groupB = unlist(strsplit(opt$groupB,",")); 
  if(!all(opt$groupA %in% group)) stop("Not all groupA belongs to group.")
  if(!all(opt$groupB %in% group)) stop("Not all groupB belongs to group.")
  if(length(groupA) != length(groupB)) stop("groupA and groupB should have the same length.")
}

### samplesize check##########################
cat("---  differential test in one variable   ---\n")
if(min(table(group_list)) < 10) {
  small_size = TRUE;
  message("min group size is less than 10. Consider as small size sample. ")
} else {
  small_size = FALSE;
  message("min group size is large than 10.")
}
if(small_size) {
  opt$do_ttest = TRUE;
  cat("Only do t test.\n")
}
if(is.na(opt$do_ttest)) {
  opt$do_ttest = FALSE;
  message("large sample size. Setting do_ttest to FALSE.")}


## t #####
recordnum = NULL
if(opt$do_ttest) {
  if(opt$varequal) message("using student t.test.") else
    message("using welch t.test.")
  ttest_result=list()
  ttest_result_all=list()
  for(i in seq_along(opt$groupB)){
    ctl = which(group_list == opt$groupA[i])
    ind=which(group_list== opt$groupB[i])
    grp <- factor(group_list[c(ctl,ind)],levels = c(opt$groupA[i],opt$groupB[i]))
    p = log2FC = NULL
    for(j in 1:nrow(quant)){
      if(sum(!is.na(quant[j,ctl])) < 3 | sum(!is.na(quant[j,ind])) < 3)
        {
        p <- c(p,NA);log2FC <- c(log2FC,NA)
      } else {
        ctlsd <- sd(quant[j,ctl],na.rm = TRUE)/mean(unlist(quant[j,ctl]));
        indsd <- sd(quant[j,ind],na.rm = TRUE)/mean(unlist(quant[j,ind]))
        if( abs(ctlsd) < 1e-6| abs(indsd)< 1e-6) {
          p <- c(p,NA);log2FC <- c(log2FC,NA)} else{
          t = t.test(unlist(quant[j,c(ctl,ind)])~grp, var.equal = opt$varequal);
          log2FC <- c(log2FC,t$estimate[2]-t$estimate[1]);
          p <- c(p,t$p.value);
        }
      }
      
      if(j == trunc(0.25*nrow(quant)) & j > 2000)
        print(paste0("....t.test:groupset:",i, " 25% completed ......"))
      if(j == trunc(0.5*nrow(quant)) & j > 4000)
        print(paste0("....t.test:groupset:",i, " 50% completed ......")) 
      if(j == trunc(0.75*nrow(quant)) & j > 4000)
        print(paste0("....t.test:groupset:",i, " 75% completed ......"))
      if(j == nrow(quant) & j > 6000)
        print(paste0("....t.test:groupset:",i, " 100% completed ......"))
    }
    p.adj <- p.adjust(p, method = "fdr")
    result <- data.frame(SYMBOL = rownames(quant),log2FC,p.value = p,p.adj)
    ttest_result_all[[i]]=result;
    ind2=which(abs(result$log2FC)> opt$log2FCcut & result$p.adj < opt$adjPcut)
    com=unique(group_list[c(ctl,ind)])
    comp=paste(com,collapse=" vs ")
    message(paste("found",length(ind2), "DEG in" ,comp))
    recordnum = paste(recordnum, "found",length(ind2), "DEG in" ,comp," t.test\n ")
    if(length(ind2) == 0) 
      ttest_result[[i]]=NA else 
        ttest_result[[i]]=result[ind2,]
    names(ttest_result_all)[i] = paste0(opt$groupA[i],"vs",opt$groupB[i])
    names(ttest_result)[i] = paste0(opt$groupA[i],"vs",opt$groupB[i])
  }
}

## wilcox #######
if(!small_size) {
  message("big sample size, using Wilcoxon test.")
  wilcox_result=list()
  wilcox_result_all=list()
  
  for(i in seq_along(opt$groupB)){
    ctl = which(group_list == opt$groupA[i])
    ind=which(group_list==opt$groupB[i])
    grp <- factor(group_list[c(ctl,ind)],levels = c(opt$groupA[i],opt$groupB[i]))
    p = log2FC = NULL
    for(j in 1:nrow(quant)){
      w = wilcox.test(unlist(quant[j,c(ctl,ind)])~grp);
      meanvalue = tapply(unlist(quant[j,c(ctl,ind)]),grp,median)
      log2FC <- c(log2FC,meanvalue[2]-meanvalue[1]);
      p <- c(p,w$p.value);
      if(j == trunc(0.25*nrow(quant)) & j > 2000)
        print(paste0("....Wilcoxon: groupset:",i, " 25% completed ......"))
      if(j == trunc(0.5*nrow(quant)) & j > 4000)
        print(paste0("....Wilcoxon: groupset:",i, " 50% completed ......")) 
      if(j == trunc(0.75*nrow(quant)) & j > 4000)
        print(paste0("....Wilcoxon:groupset:",i, " 75% completed ......"))
      if(j == nrow(quant) & j > 6000)
        print(paste0("....Wilcoxon:groupset:",i, " 100% completed ......"))
    }
    p.adj <- p.adjust(p, method = "fdr")
    result <- data.frame(SYMBOL = rownames(quant),log2FC,p.value = p,p.adj)
    
    wilcox_result_all[[i]]=result;
    ind2=which(abs(result$log2FC)> opt$log2FCcut & result$p.adj< opt$adjPcut)
    com=unique(group_list[c(ctl,ind)])
    comp=paste(com,collapse=" vs ")
    message(paste("found",length(ind2), "DEG in" ,comp))
    recordnum = paste(recordnum, "found",length(ind2), "DEG in" ,comp," wilcox test\n ")
    
    if(length(ind2) == 0) 
      wilcox_result[[i]]=NA else 
        wilcox_result[[i]]=result[ind2,]
    names(wilcox_result_all)[i] = paste0(opt$groupA[i],"vs",opt$groupB[i])
    names(wilcox_result)[i] = paste0(opt$groupA[i],"vs",opt$groupB[i])
  }
}

##save ########
#do_ttest can be setting as T,F,NA
#small_size T,F depend on sample size
#if small_size is T, force to setting do_ttest T, and only do ttest
#if small_size is F, (do_ttest T, both do t & wilcoxn)|(only do wilcoxn)
if(small_size) {
  record = if(opt$varequal) paste0("smallsize, only do student t test.\n ") else paste0("smallsize, only do welch t test.\n ") 
} else if(opt$do_ttest) {
  record = if(opt$varequal) paste0("Do student t test and wilcox test.\n ") else paste0("Do welch t test and wilcox test.\n ") 
} else record = "Do wilcox test.\n "
record2 = if("ttest_result_all" %in% ls()) paste(names(ttest_result_all),collapse = " and ") else paste(names(wilcox_result_all),collapse = " and ")
loginfo1 = list(record = paste0(record,
                                record2, "is compared.\n ",
                                recordnum),
                input = opt$input,
                output = paste0(opt$output,"_diffstat_classic.RData"),
                sessionInfo = sessionInfo(),
                srcipt = "step3_DEG_classic v2.0", parameter = opt,
                running_time = Sys.time())
loginfo = c(loginfo, DEG_classic = list(loginfo1))

if(small_size) {
  save(ttest_result,loginfo,file = paste0(opt$output,"_DEG_classic.RData"))
  save(ttest_result_all,loginfo,
       file = paste0(opt$output,"_diffstat_classic.RData"))
  cat("--Finish DEG analysis--\n")
} else {
    if(opt$do_ttest) {
      save(ttest_result,wilcox_result,loginfo,
           file = paste0(opt$output,"_DEG_classic.RData"))
      save(ttest_result_all,wilcox_result_all,loginfo,
           file = paste0(opt$output,"_diffstat_classic.RData"))
      cat("--Finish DEG analysis, do t test and wilcoxon test--\n")
      } else {
      save(wilcox_result,loginfo,
           file = paste0(opt$output,"_DEG_classic.RData"))
      save(wilcox_result_all,loginfo,
           file = paste0(opt$output,"_diffstat_classic.RData"))
      cat("--Finish DEG analysis, only do wilcoxon test--\n")
      }
}
cat("--- saved in ",opt$output,"_DEG_classic.RData and ",opt$output,"_diffstat_classic.RData ---\n",sep = "")