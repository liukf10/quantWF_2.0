

library(clusterProfiler)
term2gene <- function(IDlist,universe = NULL) {
  if(is.list(IDlist)) {
    if(is.null(names(IDlist))) names(IDlist) = 1:length(IDlist)
    genelist = NULL
    for(i in 1:length(IDlist)) {
      if(is.null(genelist)) genelist =  data.frame(TERM = names(IDlist)[i], gene = IDlist[[i]]) else
        genelist <- rbind(genelist,
                          data.frame(TERM = names(IDlist)[i], gene = IDlist[[i]]))
    }
  } else if(is.vector(IDlist)) {
    genelist =  data.frame(TERM = "list1", gene = IDlist)} else
      stop("wrong IDlist class.")
  if(!is.null(universe)) {
    genelist <- rbind(genelist,
                      data.frame(TERM = "universe", gene = universe))
  }
  genelist
}
GS2GSenrich <- function(genelist1,genelist2,universe=NULL) {
  genelist1 = genelist1[!is.na(genelist1$gene),]
  genelist1 = genelist1[genelist1$gene != "",]
  genelist2 = genelist2[!is.na(genelist2$gene),]
  genelist2 = genelist2[genelist2$gene != "",]
  if("universe" %in% genelist2$TERM) {
    universe = genelist2$gene[genelist2$TERM == "universe"]
    maxGSSize = max(table(genelist2$TERM[genelist2$TERM != "universe"]))
    term2 = unique(genelist2$TERM)[unique(genelist2$TERM) != "universe"]
  } else {
    maxGSSize = max(table(genelist2$TERM));
    term2 = unique(genelist2$TERM)}
  Z = data.frame(matrix(0,ncol = length(term2),nrow = nlevels(as.factor(genelist1$TERM))))
  colnames(Z) <- term2
  rownames(Z) <- unique(genelist1$TERM)
  count = p = enrichfold = Z
  for(i in rownames(Z)) {
    eOther <- enricher(gene = subset(genelist1,TERM == i)$gene,
                       universe = universe,
                       TERM2GENE=genelist2, 
                       maxGSSize =maxGSSize )@result
    eOther$enrichFold = as.numeric(gsub("/.*","",eOther$GeneRatio))/as.numeric(gsub(".*/","",eOther$GeneRatio))/(as.numeric(gsub("/.*","",eOther$BgRatio))/as.numeric(gsub(".*/","",eOther$BgRatio)))
    eOther$Z = -log10(eOther$pvalue)
    Z[i,match(eOther$ID,colnames(Z))] <- eOther$Z
    enrichfold[i,match(eOther$ID,colnames(Z))] <- eOther$enrichFold
    p[i,match(eOther$ID,colnames(Z))] <- eOther$pvalue
    count[i,match(eOther$ID,colnames(Z))] <- eOther$Count
  }
  list(enrichfold=enrichfold, Z=Z,
       p = p, count= count)
}

HeatMapMe = function(GS2GS, transpose = FALSE, pcutoff = 0.05, 
                     filename = NULL,width = 6, height = 7,
                     cex.text = 1, xLabelsAngle = 0,xLabelsAdj = 0.5,
                     mar =  c(5, 4, 4, 2) + 0.1, plotdir = "plot"){
  if(transpose) {
    p <- data.frame(t(GS2GS$p))
    count = t(GS2GS$count);
    enrichfold = t(GS2GS$enrichfold)
    
  } else {
    p = GS2GS$p;  count = GS2GS$count;
    enrichfold = GS2GS$enrichfold
    }
  
  p[p==0] = 2; pos = p > pcutoff
  count[pos] = "";
  p <- signif(p, digits = 2); p <- format(p, scientific = TRUE)
  p[pos] = "";
  textMatrix = paste(as.matrix(count), "\n(", as.matrix(p), ")", sep = "")
  dim(textMatrix) = dim(count)
  textMatrix[textMatrix == "\n()"] <- ""
  
  # 创建 plot 目录
  if(!is.null(filename)) {
    if(!dir.exists(plotdir)) dir.create(plotdir)
    
    # 保存 PDF
    pdf(file.path(plotdir, paste0(filename, ".pdf")), width = width, height = height)
    par(mar = mar)
    WGCNA::labeledHeatmap(Matrix = enrichfold, xLabels = colnames(p), 
                          yLabels = rownames(p), cex.lab = 1, colorLabels = TRUE, 
                          colors = WGCNA::blueWhiteRed(100)[51:100], textMatrix = textMatrix, 
                          setStdMargins = FALSE, cex.text = cex.text, xLabelsAngle = xLabelsAngle,xLabelsAdj = xLabelsAdj)
    dev.off()
    
    # 保存 PNG
    png(file.path(plotdir, paste0(filename, ".png")), width = width * 100, height = height * 100, res = 100)
    par(mar = mar)
    WGCNA::labeledHeatmap(Matrix = enrichfold, xLabels = colnames(p), 
                          yLabels = rownames(p), cex.lab = 1, colorLabels = TRUE, 
                          colors = WGCNA::blueWhiteRed(100)[51:100], textMatrix = textMatrix, 
                          setStdMargins = FALSE, cex.text = cex.text, xLabelsAngle = xLabelsAngle,xLabelsAdj = xLabelsAdj)
    dev.off()
    
    cat("Heatmap saved to ", file.path(plotdir, paste0(filename, ".pdf")), " and ", 
        file.path(plotdir, paste0(filename, ".png")), "\n")
  } else {
    par(mar = mar)
    WGCNA::labeledHeatmap(Matrix = enrichfold, xLabels = colnames(p), 
                          yLabels = rownames(p), cex.lab = 1, colorLabels = TRUE, 
                          colors = WGCNA::blueWhiteRed(100)[51:100], textMatrix = textMatrix, 
                          setStdMargins = FALSE, cex.text = cex.text, xLabelsAngle = xLabelsAngle,xLabelsAdj = xLabelsAdj)
  }
}

G2GSenrich <- function(gene,genelist2,universe=NULL) {
  genelist2 = genelist2[!is.na(genelist2$gene),]
  genelist2 = genelist2[genelist2$gene != "",]
  if("universe" %in% genelist2$TERM) {
    universe = genelist2$gene[genelist2$TERM == "universe"]
    maxGSSize = max(table(genelist2$TERM[genelist2$TERM != "universe"]))
    term2 = unique(genelist2$TERM)[unique(genelist2$TERM) != "universe"]
  } else {
    maxGSSize = max(table(genelist2$TERM));
    term2 = unique(genelist2$TERM)}
  eOther <- enricher(gene = gene,universe = universe,
                     TERM2GENE=genelist2, 
                     maxGSSize =maxGSSize)
}


compareEnrichRes <- function(enrichlist,path = path) {
  load(paste0(path,"compareCluster_template.RData"))
  temp <- xx
  if(!is.list(enrichlist)) stop("enrichlist should be a list.")
  if(is.null(names(enrichlist))) stop("enrichlist should have name.")
  if(length(enrichlist) < 2) stop("enrichlist only 1 element.")
  multicluster = geneCluster =  NULL
  for(i in names(enrichlist)) {
    if(is.null(multicluster)) {
      multicluster <- data.frame(Cluster = i,as.data.frame(enrichlist[[i]]))
      geneCluster <- list(enrichlist[[i]]@gene)} else {
        multicluster <- rbind(multicluster, data.frame(Cluster = i,as.data.frame(enrichlist[[i]])))
        geneCluster <- c(geneCluster,list(enrichlist[[i]]@gene))  
    }
  }
  temp@compareClusterResult <- multicluster
  temp@compareClusterResult$Cluster <- as.factor(temp@compareClusterResult$Cluster)
  temp@geneClusters <- geneCluster
  names(temp@geneClusters) <- names(enrichlist)
  temp
}

# 新增：生成 Gene Set Enrichment 分析报告的函数
gs_enrich_report <- function(GS2GS_result, filename = "GS_enrich_report", 
                              genelist1_name = "Gene Set 1", 
                              genelist2_name = "Gene Set 2",
                              pcutoff = 0.05) {
  if(!dir.exists("plot")) dir.create("plot")
  
  # 生成热图
  HeatMapMe(GS2GS_result, filename = filename, pcutoff = pcutoff)
  
  # 统计显著富集的数量
  p_matrix <- GS2GS_result$p
  sig_count <- sum(p_matrix < pcutoff, na.rm = TRUE)
  total_tests <- length(p_matrix)
  
  # 准备 Markdown 报告内容
  md_content <- paste0("# Step 4: Gene Set to Gene Set Enrichment Analysis Report\n\n")
  
  md_content <- paste0(md_content, "## Input Information\n")
  md_content <- paste0(md_content, "- Gene Set 1: ", genelist1_name, "\n")
  md_content <- paste0(md_content, "- Gene Set 2: ", genelist2_name, "\n")
  md_content <- paste0(md_content, "\n")
  
  md_content <- paste0(md_content, "## Parameters\n")
  md_content <- paste0(md_content, "- p-value cutoff: ", pcutoff, "\n")
  md_content <- paste0(md_content, "\n")
  
  md_content <- paste0(md_content, "## Results Summary\n")
  md_content <- paste0(md_content, "- Total tests: ", total_tests, "\n")
  md_content <- paste0(md_content, "- Significant enrichments (p < ", pcutoff, "): ", sig_count, "\n")
  md_content <- paste0(md_content, "- Significance rate: ", round(sig_count/total_tests * 100, 2), "%\n")
  md_content <- paste0(md_content, "\n")
  
  md_content <- paste0(md_content, "## Plots\n\n")
  md_content <- paste0(md_content, "### Enrichment Heatmap\n")
  md_content <- paste0(md_content, "![Enrichment Heatmap](plot/", filename, ".png)\n\n")
  
  md_content <- paste0(md_content, "## Output Files\n")
  md_content <- paste0(md_content, "- PDF plot: plot/", filename, ".pdf\n")
  md_content <- paste0(md_content, "- PNG plot: plot/", filename, ".png\n")
  
  # 写入 Markdown 报告
  writeLines(md_content, paste0(filename, ".md"))
  cat("Markdown report saved to ", paste0(filename, ".md\n"))
}
