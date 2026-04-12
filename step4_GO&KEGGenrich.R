
func_enrich <- function(gene,organism,keyType,universe=NULL,filename = "1",pvalueCutoff = 0.05, qvalueCutoff = 0.2, ...) {
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
                   readable = TRUE,pvalueCutoff = pvalueCutoff, qvalueCutoff = qvalueCutoff, ...)
    CC <- enrichGO(gene,OrgDb,keyType = keyType, ont = "CC",
                   readable = TRUE,pvalueCutoff = pvalueCutoff, qvalueCutoff = qvalueCutoff, ...)
    MF <- enrichGO(gene,OrgDb,keyType = keyType, ont = "MF",
                   readable = TRUE,pvalueCutoff = pvalueCutoff, qvalueCutoff = qvalueCutoff, ...)
    if(keyType %in% "ENTREZID") {
      KEGG <- enrichKEGG(gene, organism = orgdb2, keyType = "ncbi-geneid",pvalueCutoff = pvalueCutoff, qvalueCutoff = qvalueCutoff, ...)
    } else if(keyType %in% "UNIPROT") {
      KEGG <- enrichKEGG(gene, organism = orgdb2, keyType = "uniprot",pvalueCutoff = pvalueCutoff, qvalueCutoff = qvalueCutoff, ...)
    } else {
      idtransfer1 <- bitr(gene, keyType, "ENTREZID",OrgDb);
      KEGG <- enrichKEGG(idtransfer1$ENTREZID, organism = orgdb2, keyType = "ncbi-geneid",pvalueCutoff = pvalueCutoff, qvalueCutoff = qvalueCutoff, ...)
    }
  } else{
    BP <- enrichGO(gene,OrgDb,keyType = keyType, ont = "BP",universe = universe,
                   readable = TRUE,pvalueCutoff = pvalueCutoff, qvalueCutoff = qvalueCutoff, ...)
    CC <- enrichGO(gene,OrgDb,keyType = keyType, ont = "CC",universe = universe,
                   readable = TRUE,pvalueCutoff = pvalueCutoff, qvalueCutoff = qvalueCutoff, ...)
    MF <- enrichGO(gene,OrgDb,keyType = keyType, ont = "MF",universe = universe,
                   readable = TRUE,pvalueCutoff = pvalueCutoff, qvalueCutoff = qvalueCutoff, ...)
    if(keyType %in% "ENTREZID") {
      KEGG <- enrichKEGG(gene, organism = orgdb2, keyType = "ncbi-geneid",universe = universe,pvalueCutoff = pvalueCutoff, qvalueCutoff = qvalueCutoff, ...)
    } else if(keyType %in% "UNIPROT") {
      KEGG <- enrichKEGG(gene, organism = orgdb2, keyType = "uniprot",universe = universe,pvalueCutoff = pvalueCutoff, qvalueCutoff = qvalueCutoff, ...)
    } else {
      idtransfer1 <- bitr(gene, keyType, "ENTREZID",OrgDb)
      idtransfer2 <- bitr(universe, keyType, "ENTREZID",OrgDb)
      KEGG <- enrichKEGG(idtransfer1$ENTREZID, organism = orgdb2, 
                         keyType = "ncbi-geneid",universe = idtransfer2$ENTREZID,pvalueCutoff = pvalueCutoff, qvalueCutoff = qvalueCutoff, ...)
    }
  }
  
  # 统计结果数量
  bp_count <- ifelse(is.null(BP) || nrow(as.data.frame(BP)) == 0, 0, nrow(as.data.frame(BP)))
  cc_count <- ifelse(is.null(CC) || nrow(as.data.frame(CC)) == 0, 0, nrow(as.data.frame(CC)))
  mf_count <- ifelse(is.null(MF) || nrow(as.data.frame(MF)) == 0, 0, nrow(as.data.frame(MF)))
  kegg_count <- ifelse(is.null(KEGG) || nrow(as.data.frame(KEGG)) == 0, 0, nrow(as.data.frame(KEGG)))
  
  has_result <- any(bp_count > 0 | cc_count > 0 | mf_count > 0 | kegg_count > 0)
  
  if(has_result){
    # 创建 plot 目录
    if(!dir.exists("plot")) dir.create("plot")
    
    cat("draw dotplot and save in ",paste0("plot/func_enrich_",filename,".pdf and .png\n"))
    
    # 准备 Markdown 报告内容
    md_content <- paste0("# Step 4: GO/KEGG Enrichment Analysis Report\n\n")
    md_content <- paste0(md_content, "## Input Information\n")
    md_content <- paste0(md_content, "- Organism: ", organism, "\n")
    md_content <- paste0(md_content, "- Key type: ", keyType, "\n")
    md_content <- paste0(md_content, "- Number of input genes: ", length(gene), "\n")
    if(!is.null(universe)) {
      md_content <- paste0(md_content, "- Universe size: ", length(universe), "\n")
    }
    md_content <- paste0(md_content, "\n")
    
    md_content <- paste0(md_content, "## Parameters\n")
    md_content <- paste0(md_content, "- p-value cutoff: ", pvalueCutoff, "\n")
    md_content <- paste0(md_content, "- q-value cutoff: ", qvalueCutoff, "\n")
    md_content <- paste0(md_content, "\n")
    
    md_content <- paste0(md_content, "## Results Summary\n")
    md_content <- paste0(md_content, "- GO BP: ", bp_count, " enriched terms\n")
    md_content <- paste0(md_content, "- GO CC: ", cc_count, " enriched terms\n")
    md_content <- paste0(md_content, "- GO MF: ", mf_count, " enriched terms\n")
    md_content <- paste0(md_content, "- KEGG: ", kegg_count, " enriched pathways\n")
    md_content <- paste0(md_content, "\n")
    
    md_content <- paste0(md_content, "## Plots\n\n")
    
    # 保存 PDF
    pdf(paste0("plot/func_enrich_",filename,".pdf"))
    if(bp_count > 0) { print(dotplot(BP, title = "GO BP Enrichment"));}
    if(cc_count > 0) { print(dotplot(CC, title = "GO CC Enrichment"));}
    if(mf_count > 0) { print(dotplot(MF, title = "GO MF Enrichment"));}
    if(kegg_count > 0) { print(dotplot(KEGG, title = "KEGG Enrichment"));}
    dev.off()
    
    # 保存 PNG
    if(bp_count > 0) {
      png(paste0("plot/func_enrich_",filename,"_BP.png"), width = 800, height = 600, res = 100)
      print(dotplot(BP, title = "GO BP Enrichment"))
      dev.off()
      md_content <- paste0(md_content, "### GO BP Enrichment\n")
      md_content <- paste0(md_content, "![GO BP](plot/func_enrich_",filename,"_BP.png)\n\n")
    }
    if(cc_count > 0) {
      png(paste0("plot/func_enrich_",filename,"_CC.png"), width = 800, height = 600, res = 100)
      print(dotplot(CC, title = "GO CC Enrichment"))
      dev.off()
      md_content <- paste0(md_content, "### GO CC Enrichment\n")
      md_content <- paste0(md_content, "![GO CC](plot/func_enrich_",filename,"_CC.png)\n\n")
    }
    if(mf_count > 0) {
      png(paste0("plot/func_enrich_",filename,"_MF.png"), width = 800, height = 600, res = 100)
      print(dotplot(MF, title = "GO MF Enrichment"))
      dev.off()
      md_content <- paste0(md_content, "### GO MF Enrichment\n")
      md_content <- paste0(md_content, "![GO MF](plot/func_enrich_",filename,"_MF.png)\n\n")
    }
    if(kegg_count > 0) {
      png(paste0("plot/func_enrich_",filename,"_KEGG.png"), width = 800, height = 600, res = 100)
      print(dotplot(KEGG, title = "KEGG Enrichment"))
      dev.off()
      md_content <- paste0(md_content, "### KEGG Enrichment\n")
      md_content <- paste0(md_content, "![KEGG](plot/func_enrich_",filename,"_KEGG.png)\n\n")
    }
    
    md_content <- paste0(md_content, "## Output Files\n")
    md_content <- paste0(md_content, "- Enrichment results: list containing BP, CC, MF, KEGG\n")
    md_content <- paste0(md_content, "- PDF plots: plot/func_enrich_",filename,".pdf\n")
    md_content <- paste0(md_content, "- PNG plots: plot/func_enrich_",filename,"_*.png\n")
    
    # 写入 Markdown 报告
    writeLines(md_content, paste0("func_enrich_",filename,".md"))
    cat("Markdown report saved to ", paste0("func_enrich_",filename,".md\n"))
  }
  enrichlist <- list(BP = BP, CC = CC, MF = MF,KEGG = KEGG)
}

func_enrich_list <- function(genelist, TERM=NULL,organism, keyType, universe = NULL, filename = "",pvalueCutoff = 0.05, qvalueCutoff = 0.2, ...) {
  if(!all(c("TERM","gene") %in% colnames(genelist))) 
    stop("genelist should have TERM and gene column.\n It can generate by term2gene function.")
  BP=CC=MF=KEGG=NULL
  if(length(unique(genelist$TERM)) < 2) stop("genelist should have at least 2 TERM.")
  if(is.null(TERM)) TERM <- unique(genelist$TERM)
  for(i in TERM) {
    x <- subset(genelist,TERM == i)$gene
    efunc <- func_enrich(x,organism=organism,keyType=keyType,universe=universe,
                         filename = paste0(filename,"_",i),pvalueCutoff = pvalueCutoff, qvalueCutoff = qvalueCutoff, ...)
    BP <- c(BP,efunc$BP)
    CC <- c(CC,efunc$CC)
    MF <- c(MF,efunc$MF)
    KEGG <- c(KEGG,efunc$KEGG)
  }
  names(BP) <- names(CC) <- names(MF) <- names(KEGG) <- TERM
  enrichlist <- list(BP = BP, CC = CC, MF = MF,KEGG = KEGG)
}



gsea_enrich <- function(genelist,organism,keyType,filename = "1",pvalueCutoff = 0.05, ...){
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
  
  BP <- gseGO(genelist,OrgDb,keyType = keyType, ont = "BP",pvalueCutoff = pvalueCutoff, ...)
  CC <- gseGO(genelist,OrgDb,keyType = keyType, ont = "CC",pvalueCutoff = pvalueCutoff, ...)
  MF <- gseGO(genelist,OrgDb,keyType = keyType, ont = "MF",pvalueCutoff = pvalueCutoff, ...)
  if(keyType %in% "ENTREZID") {
    KEGG <- gseKEGG(genelist, organism = orgdb2, keyType = "ncbi-geneid",pvalueCutoff = pvalueCutoff, ...)
  } else if(keyType %in% "UNIPROT") {
    KEGG <- gseKEGG(genelist, organism = orgdb2, keyType = "uniprot",pvalueCutoff = pvalueCutoff, ...)
  } else {
    idtransfer1 <- bitr(names(genelist), keyType, "ENTREZID",OrgDb);
    genelist2 <- genelist[names(genelist) %in% idtransfer1[,1]]
    names(genelist2) <- idtransfer1$ENTREZID[match(names(genelist2),idtransfer1[,1])]
    KEGG <- gseKEGG(genelist2, organism = orgdb2, keyType = "ncbi-geneid",pvalueCutoff = pvalueCutoff, ...)
  }
  
  # 统计结果数量
  bp_count <- ifelse(is.null(BP) || nrow(as.data.frame(BP)) == 0, 0, nrow(as.data.frame(BP)))
  cc_count <- ifelse(is.null(CC) || nrow(as.data.frame(CC)) == 0, 0, nrow(as.data.frame(CC)))
  mf_count <- ifelse(is.null(MF) || nrow(as.data.frame(MF)) == 0, 0, nrow(as.data.frame(MF)))
  kegg_count <- ifelse(is.null(KEGG) || nrow(as.data.frame(KEGG)) == 0, 0, nrow(as.data.frame(KEGG)))
  total_count <- bp_count + cc_count + mf_count + kegg_count
  
  has_result <- any(bp_count > 0 | cc_count > 0 | mf_count > 0 | kegg_count > 0)
  
  if(has_result){
    # 创建 plot 目录
    if(!dir.exists("plot")) dir.create("plot")
    
    cat("draw dotplot and save in ",paste0("plot/gsea_enrich_",filename,".pdf and .png\n"))
    
    # 准备 Markdown 报告内容
    md_content <- paste0("# Step 4: GSEA Enrichment Analysis Report\n\n")
    md_content <- paste0(md_content, "## Input Information\n")
    md_content <- paste0(md_content, "- Organism: ", organism, "\n")
    md_content <- paste0(md_content, "- Key type: ", keyType, "\n")
    md_content <- paste0(md_content, "- Number of input genes: ", length(genelist), "\n")
    md_content <- paste0(md_content, "\n")
    
    md_content <- paste0(md_content, "## Parameters\n")
    md_content <- paste0(md_content, "- p-value cutoff: ", pvalueCutoff, "\n")
    md_content <- paste0(md_content, "\n")
    
    md_content <- paste0(md_content, "## Results Summary\n")
    md_content <- paste0(md_content, "- GO BP: ", bp_count, " enriched terms\n")
    md_content <- paste0(md_content, "- GO CC: ", cc_count, " enriched terms\n")
    md_content <- paste0(md_content, "- GO MF: ", mf_count, " enriched terms\n")
    md_content <- paste0(md_content, "- KEGG: ", kegg_count, " enriched pathways\n")
    md_content <- paste0(md_content, "- Total enriched: ", total_count, "\n")
    md_content <- paste0(md_content, "\n")
    
    md_content <- paste0(md_content, "## Plots\n\n")
    
    # 保存 PDF
    pdf(paste0("plot/gsea_enrich_",filename,".pdf"))
    if(bp_count > 0) { print(dotplot(BP, title = "GSEA GO BP"));}
    if(cc_count > 0) { print(dotplot(CC, title = "GSEA GO CC"));}
    if(mf_count > 0) { print(dotplot(MF, title = "GSEA GO MF"));}
    if(kegg_count > 0) { print(dotplot(KEGG, title = "GSEA KEGG"));}
    dev.off()
    
    # 保存 PNG
    if(bp_count > 0) {
      png(paste0("plot/gsea_enrich_",filename,"_BP.png"), width = 800, height = 600, res = 100)
      print(dotplot(BP, title = "GSEA GO BP"))
      dev.off()
      md_content <- paste0(md_content, "### GSEA GO BP\n")
      md_content <- paste0(md_content, "![GSEA BP](plot/gsea_enrich_",filename,"_BP.png)\n\n")
    }
    if(cc_count > 0) {
      png(paste0("plot/gsea_enrich_",filename,"_CC.png"), width = 800, height = 600, res = 100)
      print(dotplot(CC, title = "GSEA GO CC"))
      dev.off()
      md_content <- paste0(md_content, "### GSEA GO CC\n")
      md_content <- paste0(md_content, "![GSEA CC](plot/gsea_enrich_",filename,"_CC.png)\n\n")
    }
    if(mf_count > 0) {
      png(paste0("plot/gsea_enrich_",filename,"_MF.png"), width = 800, height = 600, res = 100)
      print(dotplot(MF, title = "GSEA GO MF"))
      dev.off()
      md_content <- paste0(md_content, "### GSEA GO MF\n")
      md_content <- paste0(md_content, "![GSEA MF](plot/gsea_enrich_",filename,"_MF.png)\n\n")
    }
    if(kegg_count > 0) {
      png(paste0("plot/gsea_enrich_",filename,"_KEGG.png"), width = 800, height = 600, res = 100)
      print(dotplot(KEGG, title = "GSEA KEGG"))
      dev.off()
      md_content <- paste0(md_content, "### GSEA KEGG\n")
      md_content <- paste0(md_content, "![GSEA KEGG](plot/gsea_enrich_",filename,"_KEGG.png)\n\n")
    }
    
    md_content <- paste0(md_content, "## Output Files\n")
    md_content <- paste0(md_content, "- GSEA results: list containing BP, CC, MF, KEGG\n")
    md_content <- paste0(md_content, "- PDF plots: plot/gsea_enrich_",filename,".pdf\n")
    md_content <- paste0(md_content, "- PNG plots: plot/gsea_enrich_",filename,"_*.png\n")
    
    # 写入 Markdown 报告
    writeLines(md_content, paste0("gsea_enrich_",filename,".md"))
    cat("Markdown report saved to ", paste0("gsea_enrich_",filename,".md\n"))
  }
  enrichlist <- list(BP = BP, CC = CC, MF = MF,KEGG = KEGG)
}


cat("gsea analysis: genelist is a descending order vector, name is ID, value is rank score.")


gsea_1_step <- function(stat, id = "ori.ID", log2fc = "log2FC", pvalue = "pvalue", 
                        organism = "human",keyType = "SYMBOL",filename = "1",valname,pvalueCutoff = 0.05, ...) {
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
  if(!dir.exists("plot")) dir.create("plot")
  
  # 准备 Markdown 报告内容
  md_content <- paste0("# Step 4: GSEA Analysis Report\n\n")
  md_content <- paste0(md_content, "## Input Information\n")
  md_content <- paste0(md_content, "- Input file: ", filename, "\n")
  md_content <- paste0(md_content, "- Organism: ", organism, "\n")
  md_content <- paste0(md_content, "- Key type: ", keyType, "\n")
  md_content <- paste0(md_content, "- Total genes: ", nrow(stat), "\n")
  md_content <- paste0(md_content, "\n")
  
  md_content <- paste0(md_content, "## Parameters\n")
  md_content <- paste0(md_content, "- p-value cutoff: ", pvalueCutoff, "\n")
  md_content <- paste0(md_content, "\n")
  
  result_summary <- ""
  
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
    
    gsea_log2fc <- gsea_enrich(rank_log2fc,organism,keyType,filename = paste0(filename,"_log2fc"),pvalueCutoff = pvalueCutoff, ...)
    gsea_Zscore <- gsea_enrich(rank_Zscore,organism,keyType,filename = paste0(filename,"_pvalue"),pvalueCutoff = pvalueCutoff, ...)
    gsea_deltaScore <- gsea_enrich(rank_deltaScore,organism,keyType,filename = paste0(filename,"_deltascore"),pvalueCutoff = pvalueCutoff, ...)
    
    output = list(rank = list(log2fc = rank_log2fc,
                     pvalue = rank_Zscore,
                     deltascore = rank_deltaScore),
         gsea = list(log2fc = gsea_log2fc,
                     pvalue = gsea_Zscore,
                     deltascore = gsea_deltaScore))
    
    # 统计结果数量
    total_enriched <- sum(
      sum(sapply(gsea_log2fc, function(x) ifelse(is.null(x) || nrow(as.data.frame(x)) == 0, 0, nrow(as.data.frame(x))))),
      sum(sapply(gsea_Zscore, function(x) ifelse(is.null(x) || nrow(as.data.frame(x)) == 0, 0, nrow(as.data.frame(x))))),
      sum(sapply(gsea_deltaScore, function(x) ifelse(is.null(x) || nrow(as.data.frame(x)) == 0, 0, nrow(as.data.frame(x)))))
    )
    result_summary <- paste0("- Enriched pathways (log2FC rank): ", 
                            sum(sapply(gsea_log2fc, function(x) ifelse(is.null(x) || nrow(as.data.frame(x)) == 0, 0, nrow(as.data.frame(x))))), "\n")
    result_summary <- paste0(result_summary, "- Enriched pathways (p-value rank): ", 
                            sum(sapply(gsea_Zscore, function(x) ifelse(is.null(x) || nrow(as.data.frame(x)) == 0, 0, nrow(as.data.frame(x))))), "\n")
    result_summary <- paste0(result_summary, "- Enriched pathways (deltaScore rank): ", 
                            sum(sapply(gsea_deltaScore, function(x) ifelse(is.null(x) || nrow(as.data.frame(x)) == 0, 0, nrow(as.data.frame(x))))), "\n")
    result_summary <- paste0(result_summary, "- Total enriched pathways: ", total_enriched, "\n")
    
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
    gsea_val <- gsea_enrich(rank_val,organism,keyType,filename = paste0(filename),pvalueCutoff = pvalueCutoff, ...)
    output =  list(rank = rank_val,
                   gsea = gsea_val)
    
    # 统计结果数量
    total_enriched <- sum(sapply(gsea_val, function(x) ifelse(is.null(x) || nrow(as.data.frame(x)) == 0, 0, nrow(as.data.frame(x)))))
    result_summary <- paste0("- Enriched pathways: ", total_enriched, "\n")
    
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
  
  # 完成 Markdown 报告
  md_content <- paste0(md_content, "## Results Summary\n")
  md_content <- paste0(md_content, result_summary)
  md_content <- paste0(md_content, "\n")
  
  md_content <- paste0(md_content, "## Plots\n")
  md_content <- paste0(md_content, "Plots are saved in the plot/ directory\n\n")
  
  md_content <- paste0(md_content, "## Output Files\n")
  md_content <- paste0(md_content, "- GSEA results: gsea/ directory\n")
  md_content <- paste0(md_content, "- Rank files: *.rnk\n")
  md_content <- paste0(md_content, "- PDF plots: plot/gsea_enrich_",filename,"*.pdf\n")
  md_content <- paste0(md_content, "- PNG plots: plot/gsea_enrich_",filename,"*.png\n")
  
  # 写入 Markdown 报告
  writeLines(md_content, paste0("gsea_1_step_",filename,".md"))
  cat("Markdown report saved to ", paste0("gsea_1_step_",filename,".md\n"))
  
  setwd(path)
  output
  
}
