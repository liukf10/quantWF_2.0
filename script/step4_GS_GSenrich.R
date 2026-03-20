

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
                     mar =  c(5, 4, 4, 2) + 0.1){
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
  if(!is.null(filename)) pdf(filename, width = width, height = height)
  par(mar = mar) #bottom  left top right 
  WGCNA::labeledHeatmap(Matrix = enrichfold, xLabels = colnames(p), 
                        yLabels = rownames(p), cex.lab = 1, colorLabels = TRUE, 
                        colors = WGCNA::blueWhiteRed(100)[51:100], textMatrix = textMatrix, 
                        setStdMargins = FALSE, cex.text = cex.text, xLabelsAngle = xLabelsAngle,xLabelsAdj = xLabelsAdj)
  if(!is.null(filename))  dev.off()
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

#260201 fix the bug when have no enrich res.
compareEnrichRes <- function(enrichlist,path = path) {
  load(paste0(path,"compareCluster_template.RData"))
  temp <- xx
  if(!is.list(enrichlist)) stop("enrichlist should be a list.")
  if(is.null(names(enrichlist))) stop("enrichlist should have name.")
  if(length(enrichlist) < 2) stop("enrichlist only 1 element.")
  namesterm = multicluster = geneCluster =  NULL
  for(i in names(enrichlist)) {
    if(is.null(multicluster)) {
      enres <- as.data.frame(enrichlist[[i]])
      if(nrow(enres) == 0) warning(i," have no enrich result.") else { #if no enrich result, will skip. 260201
        multicluster <- data.frame(Cluster = i,as.data.frame(enrichlist[[i]]))
        geneCluster <- list(enrichlist[[i]]@gene)
        namesterm <- i
      }
    } else {
      enres <- as.data.frame(enrichlist[[i]])
      if(nrow(enres) == 0) warning(i," have no enrich result.") else { #260201
        multicluster <- rbind(multicluster, data.frame(Cluster = i,as.data.frame(enrichlist[[i]])))
        geneCluster <- c(geneCluster,list(enrichlist[[i]]@gene)) 
        namesterm <- c(namesterm,i) #260201
      }
    }
  }
  temp@compareClusterResult <- multicluster
  temp@compareClusterResult$Cluster <- as.factor(temp@compareClusterResult$Cluster)
  temp@geneClusters <- geneCluster
  names(temp@geneClusters) <- namesterm #260201
  temp
}

