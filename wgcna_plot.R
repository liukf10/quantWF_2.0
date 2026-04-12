#Version: 2.1
#created: Kefu liu(liukefu19@163.com)
#Maintainer: Kefu Liu
#Date: March 25, 2026
#Software: R
#V2.1: 添加outputdir参数支持指定输出目录，同时生成PDF和PNG格式图表


objs <- ls();

if("input" %in% objs) {
  load(paste0(input,"_WGCNAnet.Rdata"))
} else if(!all(c("CoExpNet","MEinf","moduleinf") %in% objs))
  stop("No one of CoExpNet, MEinf or moduleinf.")

#V2.1: 使用outputdir参数创建plot目录
outputdir <- ifelse(exists("outputdir") && !is.na(outputdir), outputdir, ".")
plot_dir <- file.path(outputdir, "plot")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

#plot1
moduleLabels = CoExpNet$colors
moduleColors = labels2colors(moduleLabels)

#plot1.2
MCol <- data.frame(moduleColors, unmerged = labels2colors(CoExpNet$unmergedColors))


#plot2
library(stringr)
MEs_col = MEinf
colnames(MEs_col) = paste0("ME", labels2colors(
  as.numeric(str_replace_all(colnames(MEinf),"ME",""))))
#MEs_col <- MEs_col[,c(13,11,6,3,1,8,5,7,12,4,9,2,10)]
MEs =  orderMEs(MEinf)


if("covar" %in% objs) {meta <- covar; message("using covar to draw module_phenotype_relationship.")}
if("meta" %in% ls()) {
  for (i in 1:(ncol(meta)+1)) {
    if(i == ncol(meta)+1) {
      if(all(rownames(MEinf) %in% rownames(meta)))
        break else i = 0;
    }
    if(all(rownames(MEinf) %in% meta[,i]))
      break
  }
  if(i == 0) stop("Not all sample match with covarfile.")
  if(i == (ncol(meta)+1)) 
    pos <- match(rownames(MEinf),rownames(meta)) else {
      pos <- match(rownames(MEinf),meta[,i]);
      meta <- meta[pos,]; meta <- meta[-i]; 
      rownames(meta) <- rownames(MEinf);
    }
  message("Extract meta info, and draw module_pheotype_relationship plot.")
  #plot3
  #Module and meta relationship
  for(i in 1:ncol(meta))  meta[,i]<-as.numeric(factor(meta[,i]));
  #MEs <- MEinf[,c(13,14,8,4,9,10,1,2,11,3,7,5,12,6,15)];
  modmodCor = cor(MEs, meta, use = "p")
  modmodP = corPvalueStudent(modmodCor, nrow(MEs));
  modPadj = matrix(p.adjust(modmodP,method = "fdr"),
                   nrow = nrow(modmodP), ncol = ncol(modmodP));
  rownames(modPadj) <- rownames(modmodP); colnames(modPadj) <- colnames(modmodP);
  
  textMatrix = paste(signif(modmodCor, 2), "\n(", signif(modPadj, 1), ")", sep = "")
  dim(textMatrix) = dim(modmodCor)
  textMatrix[modPadj>0.05]<-""
} else {
  message("No meta info, and no module_pheotype_relationship plot.")
  }



# Save as PDF
pdf(file = file.path(plot_dir, paste0(output,"_WGCNA.pdf")),
    width = 20, height = 10)
plotDendroAndColors(CoExpNet$dendrograms[[1]], 
                    moduleColors[CoExpNet$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
plotDendroAndColors(CoExpNet$dendrograms[[1]], 
                    MCol[CoExpNet$blockGenes[[1]],],
                    groupLabels = NULL,
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
plotEigengeneNetworks(MEs_col, "Eigengene adjacency heatmap", 
                      marDendro = c(3,3,2,4), greyLabel = 0,
                      marHeatmap = c(3,4,2,2), plotDendrograms = F, 
                      xLabelsAngle = 90)
plotEigengeneNetworks(MEs, "Eigengene adjacency heatmap", 
                      marDendro = c(3,3,2,4),
                      marHeatmap = c(3,4,2,2), plotDendrograms = T, 
                      xLabelsAngle = 90)
if("meta" %in% ls()) {
  par(mar = c(5, 6, 3, 4));
  labeledHeatmap(Matrix = t(modmodCor), xLabels = colnames(MEs), 
                 yLabels = colnames(meta), 
                 cex.lab = 0.7, xLabelsAngle = 0, xLabelsAdj = 0.5,
                 ySymbols = colnames(meta), colorLabels = FALSE, 
                 colors = blueWhiteRed(50), 
                 textMatrix = t(textMatrix), setStdMargins = FALSE, 
                 cex.text = 0.6, zlim = c(-1,1),
                 main = paste("Module-trait relationships"))
}
graphics.off()

# Save as PNG - each plot as a separate file
# 1. Module dendrogram and colors
png(file = file.path(plot_dir, paste0(output,"_WGCNA_module_dendro.png")),
    width = 12, height = 8, units = "in", res = 300)
plotDendroAndColors(CoExpNet$dendrograms[[1]], 
                    moduleColors[CoExpNet$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
graphics.off()

# 2. Module dendrogram with merged colors
png(file = file.path(plot_dir, paste0(output,"_WGCNA_merged_colors.png")),
    width = 12, height = 8, units = "in", res = 300)
plotDendroAndColors(CoExpNet$dendrograms[[1]], 
                    MCol[CoExpNet$blockGenes[[1]],],
                    groupLabels = NULL,
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
graphics.off()

# 3. Eigengene adjacency heatmap (no dendrograms)
png(file = file.path(plot_dir, paste0(output,"_WGCNA_eigengene_heatmap.png")),
    width = 10, height = 8, units = "in", res = 300)
plotEigengeneNetworks(MEs_col, "Eigengene adjacency heatmap", 
                      marDendro = c(3,3,2,4), greyLabel = 0,
                      marHeatmap = c(3,4,2,2), plotDendrograms = F, 
                      xLabelsAngle = 90)
graphics.off()

# 4. Eigengene adjacency heatmap (with dendrograms)
png(file = file.path(plot_dir, paste0(output,"_WGCNA_eigengene_dendro.png")),
    width = 10, height = 8, units = "in", res = 300)
plotEigengeneNetworks(MEs, "Eigengene adjacency heatmap", 
                      marDendro = c(3,3,2,4),
                      marHeatmap = c(3,4,2,2), plotDendrograms = T, 
                      xLabelsAngle = 90)
graphics.off()

# 5. Module-trait relationships (if meta exists)
if("meta" %in% ls()) {
  png(file = file.path(plot_dir, paste0(output,"_WGCNA_module_trait.png")),
      width = 12, height = 10, units = "in", res = 300)
  par(mar = c(5, 6, 3, 4));
  labeledHeatmap(Matrix = t(modmodCor), xLabels = colnames(MEs), 
                 yLabels = colnames(meta), 
                 cex.lab = 0.7, xLabelsAngle = 0, xLabelsAdj = 0.5,
                 ySymbols = colnames(meta), colorLabels = FALSE, 
                 colors = blueWhiteRed(50), 
                 textMatrix = t(textMatrix), setStdMargins = FALSE, 
                 cex.text = 0.6, zlim = c(-1,1),
                 main = paste("Module-trait relationships"))
  graphics.off()
}

# Also save individual module-trait heatmap as separate files if meta exists
if("meta" %in% ls()) {
  # Save module-trait heatmap as separate PDF
  pdf(file = file.path(plot_dir, paste0("Module_trait_heatmap_", output, ".pdf")),
      width = 12, height = 10)
  par(mar = c(5, 8, 3, 4));
  labeledHeatmap(Matrix = t(modmodCor), xLabels = colnames(MEs), 
                 yLabels = colnames(meta), 
                 cex.lab = 0.8, xLabelsAngle = 0, xLabelsAdj = 0.5,
                 ySymbols = colnames(meta), colorLabels = FALSE, 
                 colors = blueWhiteRed(50), 
                 textMatrix = t(textMatrix), setStdMargins = FALSE, 
                 cex.text = 0.7, zlim = c(-1,1),
                 main = paste("Module-trait relationships"))
  graphics.off()
  
  # Save module-trait heatmap as separate PNG
  png(file = file.path(plot_dir, paste0("Module_trait_heatmap_", output, ".png")),
      width = 12, height = 10, units = "in", res = 300)
  par(mar = c(5, 8, 3, 4));
  labeledHeatmap(Matrix = t(modmodCor), xLabels = colnames(MEs), 
                 yLabels = colnames(meta), 
                 cex.lab = 0.8, xLabelsAngle = 0, xLabelsAdj = 0.5,
                 ySymbols = colnames(meta), colorLabels = FALSE, 
                 colors = blueWhiteRed(50), 
                 textMatrix = t(textMatrix), setStdMargins = FALSE, 
                 cex.text = 0.7, zlim = c(-1,1),
                 main = paste("Module-trait relationships"))
  graphics.off()
}

message(paste0("WGCNA plots saved to: ", plot_dir))
message("Generated files:")
message(paste0("  - ", output, "_WGCNA.pdf"))
message(paste0("  - ", output, "_WGCNA_module_dendro.png"))
message(paste0("  - ", output, "_WGCNA_merged_colors.png"))
message(paste0("  - ", output, "_WGCNA_eigengene_heatmap.png"))
message(paste0("  - ", output, "_WGCNA_eigengene_dendro.png"))
if("meta" %in% ls()) {
  message(paste0("  - ", output, "_WGCNA_module_trait.png"))
  message(paste0("  - Module_trait_heatmap_", output, ".pdf"))
  message(paste0("  - Module_trait_heatmap_", output, ".png"))
}
