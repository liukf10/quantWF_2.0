#Version: 2.0
#created: Kefu liu(liukefu19@163.com)
#Maintainer: Kefu Liu
#Date: July 25, 2025
#Software: R

#input \t seperate
### package check ########
existpck <- rownames(installed.packages());
if(!"data.table" %in% existpck) 
  stop("no data.table package, please install the R package before run me.")



## object check ############
allobj <- ls();
if(!"quantfile" %in% allobj) 
  quantfile = NA;
if(!"covarfile" %in% allobj) 
  covarfile = NA;
if(!"output" %in% allobj) 
  output = NA;
if(!"IDinfCol" %in% allobj) 
  IDinfCol = NA;
if(!"NAvalue" %in% allobj)
  NAvalue = NA

opt = list(quantfile = quantfile, covarfile = covarfile, output = output,
           IDinfCol = IDinfCol, NAvalue =NAvalue)
rmobj <- c(names(opt),"opt","option_list","allobj","allobj2","quant","loginfo","existpck","funcpath")
allobj <- ls(); allobj <- allobj[!allobj %in% c("rmobj",rmobj)]

### 0.1 parameters evaluate ######
if(is.na(opt$output)) {
  opt$output = "data";
  message("no output, setting to 'data' as prefix")
}
if(all(is.na(opt$IDinfCol))) {
  opt$IDinfCol <- NA
  saveother = TRUE
  message("no IDinfCol, it will put all col without intensity into IDinf.")
} else {
  saveother = FALSE
  opt$IDinfCol <- unlist(strsplit(opt$IDinfCol,","));
  if(!any(is.na(as.numeric(opt$IDinfCol)))) opt$IDinfCol <- as.numeric(opt$IDinfCol)
}
if(is.na(opt$NAvalue)) {
  message("NAvalue setting to NA.")
}

### 0.2 read quant & covar info ####
#no quantfile will stop
if(is.na(opt$quantfile)) {
    stop ("no quantfile, see help info");
} else {
  if(!file.exists(opt$quantfile)) 
    stop(opt$quantfile, " is not exist in the dir.")
}

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

### 1.data extract ######
cat('-----raw MS data extract----- \n')
require(data.table)
rawdata <- fread(opt$quantfile, data.table = FALSE)
pos <- grep("(I|i)ntensity",colnames(rawdata))
if(length(pos) == 0) {
  pos <- grep("(A|a)bundance",colnames(rawdata))
  if(length(pos) == 0) {stop("No intensity or abundance in input colname.")}
  message("using abundance to extract the data.")
}
intensity <- rawdata[,pos]
if(!is.numeric(as.matrix(intensity))) 
  warning("intensity are not all numeric, please check before the later process.")
colnames(intensity) <- gsub("(I|i)ntensity","",colnames(intensity))
colnames(intensity) <- gsub("(A|a)bundance","",colnames(intensity))

#extract inf and other
if(!saveother) {
  #should be vector
  if(is.vector(opt$IDinfCol)) {
    # if is numeric, should be column number
    if(is.numeric(opt$IDinfCol)){
      #should not overlapped with intensity column.
      #or not overlapped with column name.
      if(!all(opt$IDinfCol %in% (1:ncol(rawdata))[-pos])) 
        stop("IDinfCol is numeric, but not all match input column number or overlapped with intensity column.")
    } else {
      if(!all(opt$IDinfCol %in% colnames(rawdata)[-pos])) 
        stop("IDinfCol is col name, but not all match input column name.")
      opt$IDinfCol <- which(colnames(rawdata) %in% opt$IDinfCol)
    }
  } else stop("IDinfCol should be a vector.")
  inf <- rawdata[,opt$IDinfCol]
  other <- rawdata[,-c(pos,opt$IDinfCol)]
} else {
  inf <- rawdata[,-pos]
  other <- NULL
}
if(!is.na(opt$NAvalue)) intensity[intensity %in% opt$NAvalue] <- NA;

IDinf <- data.frame(ID = inf[,1],inf);
quant <- intensity;

#put the first column as ID as quant rowname.
IDinf$ID <- make.names(IDinf$ID)
#modified the replication name as _1,_2... and make.names to make it readable
if(sum(duplicated(IDinf$ID)) > 0) {
  IDo<- data.frame(table(IDinf$ID))
  IDo <- subset(IDo, Freq > 1)
  ID <- IDinf$ID
  for(i in 1:nrow(IDo)) {
    var <- IDo$Var1[i]
    ID[IDinf$ID == var] <- paste0(ID[IDinf$ID == var],"_",1:IDo$Freq[i])
  }
  IDinf$ID <- ID
}
rownames(quant) <- IDinf$ID

loginfo1 = list(record = paste0("Extract the MS data.\n"),
                input = if(is.na(opt$covarfile)) opt$quantfile else paste(opt$quantfile,"and",opt$covarfile),
                output = paste0(opt$output,"_raw.RData"),
                packageVersion =paste("data.table",packageVersion("data.table")),
                sessionInfo = sessionInfo(),
                srcipt = "step0_MS v2.0", parameter = opt,
                running_time = Sys.time())
loginfo = c(loginfo, data_extract = list(loginfo1))

if(is.null(other)) 
  save(IDinf,quant,loginfo,list=otherobj, file = paste0(opt$output,"_raw.RData")) else
    save(IDinf,quant,other,loginfo,list=otherobj, file = paste0(opt$output,"_raw.RData"))

cat('---rawdata is extracted--- \n')
n1 = ncol(quant); nr1 <- nrow(quant);
cat(paste0("raw data have ", n1, " sample X ", nr1," rows \n"))
