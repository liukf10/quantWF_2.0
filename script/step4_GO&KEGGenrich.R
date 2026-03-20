
func_enrich <- function(gene,organism,keyType,universe=NULL,filename = "1",...) {
  #if(!"organism" %in% ls()) {stop("No organism. It must be human, mouse or rat.")}
  if(missing(organism)) {organism = "human"; 
  warning("organism is not setting, setting to human")}
  organism = match.arg(organism, c("human", "mouse", "rat"))
  if(missing(keyType)) {keyType = "SYMBOL"; 
  warning("keyType is not setting, setting to SYMBOL")}
  if(organism == "human") {
    library(org.Hs.eg.db);  OrgDb = "org.Hs.eg.db"; orgdb2 <- "hsa";
  } 
  if(organism == "mouse") {
    library(org.Mm.eg.db);  OrgDb = "org.Mm.eg.db"; orgdb2 <- "mmu";
  }
  if(organism == "rat") {
    library(org.Rn.eg.db);  OrgDb = "org.Rn.eg.db"; orgdb2 <- "rno";
  }
  library(clusterProfiler)
  if(!keyType %in% idType(OrgDb = OrgDb)) stop("keyType is wrong, it can be search by idType.")
  if(is.null(universe)) {
    BP <- enrichGO(gene,OrgDb,keyType = keyType, ont = "BP",
                   readable = TRUE,...)
    CC <- enrichGO(gene,OrgDb,keyType = keyType, ont = "CC",
                   readable = TRUE,...)
    MF <- enrichGO(gene,OrgDb,keyType = keyType, ont = "MF",
                   readable = TRUE,...)
    if(keyType %in% "ENTREZID") {
      KEGG <- enrichKEGG(gene, organism = orgdb2, keyType = "ncbi-geneid",...)
    } else if(keyType %in% "UNIPROT") {
      KEGG <- enrichKEGG(gene, organism = orgdb2, keyType = "uniprot",...)
    } else {
      idtransfer1 <- bitr(gene, keyType, "ENTREZID",OrgDb);
      KEGG <- enrichKEGG(idtransfer1$ENTREZID, organism = orgdb2, keyType = "ncbi-geneid",...)
    }
  } else{
    BP <- enrichGO(gene,OrgDb,keyType = keyType, ont = "BP",universe = universe,
                   readable = TRUE,...)
    CC <- enrichGO(gene,OrgDb,keyType = keyType, ont = "CC",universe = universe,
                   readable = TRUE,...)
    MF <- enrichGO(gene,OrgDb,keyType = keyType, ont = "MF",universe = universe,
                   readable = TRUE,...)
    if(keyType %in% "ENTREZID") {
      KEGG <- enrichKEGG(gene, organism = orgdb2, keyType = "ncbi-geneid",universe = universe,...)
    } else if(keyType %in% "UNIPROT") {
      KEGG <- enrichKEGG(gene, organism = orgdb2, keyType = "uniprot",universe = universe,...)
    } else {
      idtransfer1 <- bitr(gene, keyType, "ENTREZID",OrgDb)
      idtransfer2 <- bitr(universe, keyType, "ENTREZID",OrgDb)
      KEGG <- enrichKEGG(idtransfer1$ENTREZID, organism = orgdb2, 
                         keyType = "ncbi-geneid",universe = idtransfer2$ENTREZID,...)
    }
  }
  if(any(nrow(as.data.frame(BP)) > 0 | nrow(as.data.frame(CC)) > 0 | nrow(as.data.frame(MF)) > 0 | nrow(as.data.frame(KEGG)) > 0)){
    cat("draw dotplot and save in ",paste0("func_enrich_",filename,".pdf"))
    pdf(paste0("func_enrich_",filename,".pdf"))
    if(nrow(as.data.frame(BP)) > 0) { print(dotplot(BP));}
    if(nrow(as.data.frame(CC)) > 0) { print(dotplot(CC));}
    if(nrow(as.data.frame(MF)) > 0) { print(dotplot(MF));}
    if(nrow(as.data.frame(KEGG)) > 0) { print(dotplot(KEGG));}
    dev.off()
  }
  enrichlist <- list(BP = BP, CC = CC, MF = MF,KEGG = KEGG)
}

func_enrich_list <- function(genelist, TERM=NULL,organism, keyType, universe = NULL, filename = "",...) {
  if(!all(c("TERM","gene") %in% colnames(genelist))) 
    stop("genelist should have TERM and gene column.\n It can generate by term2gene function.")
  BP=CC=MF=KEGG=NULL
  if(length(unique(genelist$TERM)) < 2) stop("genelist should have at least 2 TERM.")
  if(is.null(TERM)) TERM <- unique(genelist$TERM)
  for(i in TERM) {
    x <- subset(genelist,TERM == i)$gene
    efunc <- func_enrich(x,organism=organism,keyType=keyType,universe=universe,
                         filename = paste0(filename,"_",i),...)
    BP <- c(BP,efunc$BP)
    CC <- c(CC,efunc$CC)
    MF <- c(MF,efunc$MF)
    KEGG <- c(KEGG,efunc$KEGG)
  }
  names(BP) <- names(CC) <- names(MF) <- names(KEGG) <- TERM
  enrichlist <- list(BP = BP, CC = CC, MF = MF,KEGG = KEGG)
}



gsea_enrich <- function(genelist,organism,keyType,filename = "1",...){
  if(missing(organism)) {organism = "human"; 
  warning("organism is not setting, setting to human")}
  if(missing(keyType)) {keyType = "SYMBOL"; 
  warning("keyType is not setting, setting to SYMBOL")}
  if(organism == "human") {
    library(org.Hs.eg.db);  OrgDb = "org.Hs.eg.db"; orgdb2 <- "hsa";
  } 
  if(organism == "mouse") {
    library(org.Mm.eg.db);  OrgDb = "org.Mm.eg.db"; orgdb2 <- "mmu";
  }
  if(organism == "rat") {
    library(org.Rn.eg.db);  OrgDb = "org.Rn.eg.db"; orgdb2 <- "rno";
  }
  library(clusterProfiler)
  if(!keyType %in% idType(OrgDb = OrgDb)) stop("keyType is wrong, it can be search by idType.")
  
  BP <- gseGO(genelist,OrgDb,keyType = keyType, ont = "BP",...)
  CC <- gseGO(genelist,OrgDb,keyType = keyType, ont = "CC",...)
  MF <- gseGO(genelist,OrgDb,keyType = keyType, ont = "MF",...)
  if(keyType %in% "ENTREZID") {
    KEGG <- gseKEGG(genelist, organism = orgdb2, keyType = "ncbi-geneid",...)
  } else if(keyType %in% "UNIPROT") {
    KEGG <- gseKEGG(genelist, organism = orgdb2, keyType = "uniprot",...)
  } else {
    idtransfer1 <- bitr(names(genelist), keyType, "ENTREZID",OrgDb);
    genelist2 <- genelist[names(genelist) %in% idtransfer1[,1]]
    names(genelist2) <- idtransfer1$ENTREZID[match(names(genelist2),idtransfer1[,1])]
    KEGG <- gseKEGG(genelist2, organism = orgdb2, keyType = "ncbi-geneid",...)
  }
  
  if(any(nrow(as.data.frame(BP)) > 0 | nrow(as.data.frame(CC)) > 0 | nrow(as.data.frame(MF)) > 0 | nrow(as.data.frame(KEGG)) > 0)){
    cat("draw dotplot and save in ",paste0("func_enrich_",filename,".pdf"))
    pdf(paste0("gsea_enrich_",filename,".pdf"))
    if(nrow(as.data.frame(BP)) > 0) { print(dotplot(BP));}
    if(nrow(as.data.frame(CC)) > 0) { print(dotplot(CC));}
    if(nrow(as.data.frame(MF)) > 0) { print(dotplot(MF));}
    if(nrow(as.data.frame(KEGG)) > 0) { print(dotplot(KEGG));}
    dev.off()
  }
  enrichlist <- list(BP = BP, CC = CC, MF = MF,KEGG = KEGG)
}


cat("gsea analysis: genelist is a descending order vector, name is ID, value is rank score.")


gsea_1_step <- function(stat, id = "ori.ID", log2fc = "log2FC", pvalue = "pvalue", 
                        organism = "human",keyType = "SYMBOL",filename = "1",valname,...) {
  if(organism == "human") {
    library(org.Hs.eg.db);  OrgDb = "org.Hs.eg.db"; orgdb2 <- "hsa";
  } 
  if(organism == "mouse") {
    library(org.Mm.eg.db);  OrgDb = "org.Mm.eg.db"; orgdb2 <- "mmu";
  }
  if(organism == "rat") {
    library(org.Rn.eg.db);  OrgDb = "org.Rn.eg.db"; orgdb2 <- "rno";
  }
  if(!id %in%  colnames(stat)) {
    stop("wrong id colname setting.")
  }
  if(missing(valname)) {
    if(!all(c(log2fc,pvalue) %in%  colnames(stat)))
    {stop("wrong log2fc or pvalue colname setting.")}} else {
      if(!valname %in% colnames(stat))
      {stop("wrong valname colname setting.")}
    }
  #id 转换成symbol
  if(!keyType %in% c("ENTREZID","SYMBOL")) {
    idtransfer <- bitr(stat[,id], keyType, "SYMBOL",OrgDb)
    keyType = "SYMBOL"
    stat <- stat[stat[,id] %in% idtransfer[,1],]
    stat[,id] <- idtransfer$SYMBOL[match(stat[,id],idtransfer[,1])]
  }
  stat = stat[!duplicated(stat[,id]),]
  
  if(!dir.exists("gsea")) dir.create("gsea")
  path = getwd()
  setwd("gsea")
  if(missing(valname)) {
    stat <- stat[!is.na(stat[,log2fc]),];
    stat <- stat[!is.na(stat[,pvalue]),];
    stat$Zscore <- -log10(stat[,pvalue])
    stat$deltaScore <- stat[,log2fc]*stat$Zscore
    stat$Zscore[stat[,log2fc]<0] <- -stat$Zscore[stat[,log2fc]<0]
    
    stat <- stat[order(stat[,log2fc], decreasing= TRUE),]
    rank_log2fc <- stat[,log2fc]; names(rank_log2fc) <- stat[,id]
    stat <- stat[order(stat$Zscore, decreasing= TRUE),]
    rank_Zscore <- stat$Zscore; names(rank_Zscore) <- stat[,id]
    stat <- stat[order(stat$deltaScore, decreasing= TRUE),]
    rank_deltaScore <- stat$deltaScore; names(rank_deltaScore) <- stat[,id]
    
    gsea_log2fc <- gsea_enrich(rank_log2fc,organism,keyType,filename = paste0(filename,"_log2fc"),...)
    gsea_Zscore <- gsea_enrich(rank_Zscore,organism,keyType,filename = paste0(filename,"_pvalue"),...)
    gsea_deltaScore <- gsea_enrich(rank_deltaScore,organism,keyType,filename = paste0(filename,"_deltascore"),...)
    
    output = list(rank = list(log2fc = rank_log2fc,
                     pvalue = rank_Zscore,
                     deltascore = rank_deltaScore),
         gsea = list(log2fc = gsea_log2fc,
                     pvalue = gsea_Zscore,
                     deltascore = gsea_deltaScore))
    #ENTREZID 转换成symbol
    if(keyType %in% "ENTREZID") {
      idtransfer <- bitr(stat[,id], keyType, "SYMBOL",OrgDb)
      
      rank_log2fc <- rank_log2fc[names(rank_log2fc) %in% idtransfer[,1]]
      names(rank_log2fc) <- idtransfer$SYMBOL[match(names(rank_log2fc),idtransfer[,1])]
      rank_log2fc = rank_log2fc[!duplicated(names(rank_log2fc))]
      rank_Zscore <- rank_Zscore[names(rank_Zscore) %in% idtransfer[,1]]
      names(rank_Zscore) <- idtransfer$SYMBOL[match(names(rank_Zscore),idtransfer[,1])]
      rank_Zscore = rank_Zscore[!duplicated(names(rank_Zscore))]
      rank_deltaScore <- rank_deltaScore[names(rank_deltaScore) %in% idtransfer[,1]]
      names(rank_deltaScore) <- idtransfer$SYMBOL[match(names(rank_deltaScore),idtransfer[,1])]
      rank_deltaScore = rank_deltaScore[!duplicated(names(rank_deltaScore))]
    }
    write.table(rank_log2fc,file = paste0(filename ,"_log2fc.rnk"),
                quote = FALSE, col.names = FALSE,sep = "\t")
    write.table(rank_Zscore,file = paste0(filename ,"_pvalue.rnk"),
                quote = FALSE, col.names = FALSE,sep = "\t")
    write.table(rank_deltaScore,file = paste0(filename ,"_deltascore.rnk"),
                quote = FALSE, col.names = FALSE,sep = "\t")
    
  } else {
    stat <- stat[order(stat[,valname], decreasing= TRUE),]
    rank_val <- stat[,valname]; names(rank_val) <- stat[,id];
    gsea_val <- gsea_enrich(rank_val,organism,keyType,filename = paste0(filename),...)
    output =  list(rank = rank_val,
                   gsea = gsea_val)
    
    #ENTREZID 转换成symbol
    if(keyType %in% "ENTREZID") {
      idtransfer <- bitr(stat[,id], keyType, "SYMBOL",OrgDb)
      rank_val <- rank_val[names(rank_log2fc) %in% idtransfer[,1]]
      names(rank_val) <- idtransfer$SYMBOL[match(names(rank_val),idtransfer[,1])]
      rank_val = rank_val[!duplicated(names(rank_val))]
    }
    
    write.table(rank_val,file = paste0(filename ,".rnk"),
                quote = FALSE, col.names = FALSE,sep = "\t")
    
  }
  setwd(path)
  output
  
}

