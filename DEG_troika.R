#Version: 2.1
#created: Kefu liu(liukefu19@163.com)
#Maintainer: Kefu Liu
#Date: March 25, 2026
#Software: R
#20250717 fixed limma and edgeR have ctl/ind log2FC
#20250918 add loginfo and message in last, and fix DESeq2 DEG extract error
#V2.1: 添加outputdir参数支持指定输出目录
#v2.1: 修改DESeq2输出结果，保证是用groupB/groupA的结果
#v2.1: 添加RRHO分析、火山图绘制和report.md生成功能

#######
#log2FCcut = 1; adjPcut = 0.05;
#counts; covar;grpname = "group"

#说明############
#padj为NA
#Note on p-values set to NA: some values in the results table can be set to NA for one of the following reasons:
#1. 所有样品中为0, baseMean为0,log2FC,p,padj为NA
#2.有一个样品有异常值,p,padj设置成NA. Cooks's distance 来检测outlier.
#3. 如果基因没有通过DESeq2的自动做的independent filtering(有低的mean normalized count),只有padj会设置成NA.

# object check##################
allobj <- ls();
if(!"counts" %in% allobj) stop("No counts obj, see readme document.");
if(!"covar" %in% allobj) stop("No covar obj, see readme document.");
if(!"output" %in% allobj) {
  output = "data";
  message("no output, setting to 'data' as prefix")
}
if(!"outputdir" %in% allobj) {
  outputdir = ".";
  message("no outputdir, setting to current directory")
}
if(!"funcdir" %in% allobj) {
  funcdir = NA;
  message("no funcdir, setting to NA")
}
if(!"grpname" %in% allobj) {
  grpname = "group"; message("Setting grpname as default: group.")}
if(is.vector(covar)) {
  if(length(covar) != ncol(counts)) stop("covar length is not equal with counts column number.")
  message("covar is a vector, considered in the same order with counts.")
  covar <- data.frame(samplename = colnames(counts),grpname = covar);
  colnames(covar)[2] = grpname;
}
for (i in 1:(ncol(covar)+1)) {
  if(i == ncol(covar)+1) {
    if(all(colnames(counts) %in% rownames(covar)))
      break else i = 0;
  }
  if(all(colnames(counts) %in% covar[,i]))
    break
}
if(i == 0) stop("Not all sample match with covar.")
if(i == (ncol(covar)+1)) 
  pos <- match(colnames(counts),rownames(covar)) else {
    pos <- match(colnames(counts),covar[,i]);
    covar <- covar[pos,]; covar <- covar[-i]; 
    rownames(covar) <- colnames(counts);
  } #covar & counts are matched
if(any(!grpname %in% colnames(covar)))
  stop("grpname is not in covar, please check it.")
group_list <- covar[,grpname];
group = unique(group_list)


if(!"log2FCcut" %in% allobj) {
  log2FCcut = 1; message("Setting log2FCcut as default: 1.")}
if(!"adjPcut" %in% allobj) {
  adjPcut = 0.05; message("Setting adjPcut as default: 0.05.")}
if(!"idinf_col" %in% allobj) {
  idinf_col = "SYMBOL"; message("Setting idinf_col as default: SYMBOL.")}

if(length(group) != 2) {
  if(!"groupA" %in% allobj) stop("No groupA, see readme document.")
  if(!"groupB" %in% allobj) stop("No groupB, see readme document.")
  if(!all(groupA %in% group)) stop("Not all groupA belongs to group.")
  if(!all(groupB %in% group)) stop("Not all groupB belongs to group.")
  if(length(groupA) != length(groupB)) stop("groupA and groupB should have the same length.")
} else {
  if(!"groupA" %in% allobj | !"groupB" %in% allobj) {
    groupA = group[1]; groupB = group[2];
    message(paste0("setting groupA as ",group[1],"; groupB as ",group[2]))
  }
  if(!all(groupA %in% group)) stop("Not all groupA belongs to group.")
  if(!all(groupB %in% group)) stop("Not all groupB belongs to group.")
  if(length(groupA) != length(groupB)) stop("groupA and groupB should have the same length.")
}

### Initialize directories and load functions ########
# Create plot and assets directories
plot_dir <- file.path(outputdir, "plot")
assets_dir <- file.path(outputdir, "assets")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)
if(!dir.exists(assets_dir)) dir.create(assets_dir, recursive = TRUE)

# Load preprocess_func.R for utility functions
if(is.na(funcdir)) {
  if(file.exists("preprocess_func.R")) {
    source("preprocess_func.R")
    message("Loaded preprocess_func.R from current directory")
  } else {
    message("No preprocess_func.R file in current dir, some functions may not be available.")
  }
} else {
  if(tail(unlist(strsplit(funcdir,"")),1) == "/"){
    funcpath = paste0(funcdir,"preprocess_func.R");
  } else {
    funcdir <-  paste0(funcdir,"/")
    funcpath = paste0(funcdir,"preprocess_func.R");}
  if(file.exists(funcpath)) {
    source(funcpath)
    message("Loaded preprocess_func.R from: ", normalizePath(funcpath))
  } else {
    message("No preprocess_func.R file in funcdir path, some functions may not be available.")
  }
}

###0.package check ########
existpck <- rownames(installed.packages());
if(!"limma" %in% existpck)
  stop("no limma package, please install the R pakcage before run me.")
if(!"edgeR" %in% existpck)
  stop("no edgeR package, please install the R pakcage before run me.")
if(!"DESeq2" %in% existpck)
  stop("no DESeq2 package, please install the R pakcage before run me.")

# Check for RRHO package
rrho_available <- "RRHO" %in% existpck
if(rrho_available) {
  message("RRHO package found, RRHO analysis will be performed.")
} else {
  message("RRHO package not found, skipping RRHO analysis. Install with: BiocManager::install('RRHO')")
}

##limma #####
library(limma)
library(edgeR)
limma_result=list()
limma_result_all=list()

for(i in seq_along(groupB)){
  ctl = which(group_list == groupA[i])
  ind=which(group_list==groupB[i])
  grp <- factor(group_list[c(ctl,ind)],levels = c(groupA[i],groupB[i]))
  design <- model.matrix(~0+grp)
  colnames(design) <- c(groupA[i],groupB[i])
  rownames(design) <- colnames(counts)[c(ctl,ind)]
  dge <- DGEList(counts=counts[,c(ctl,ind)])
  # calculate Norm Factor. default is TMM
  dge <- calcNormFactors(dge)
  # add 3 before logCPM, it can reduce low expression gene variation
  logCPM <- cpm(dge, log=TRUE, prior.count=3)
  v <- voom(dge,design, plot= TRUE, normalize="quantile")
  fit <- lmFit(v, design)
  com=unique(group_list[c(ind,ctl)])#250717 修改 使得上调为实验组上调
  comp=paste(com,collapse="-")
  cont.matrix <- makeContrasts(contrasts=c(comp),levels = design)
  fit2 <- contrasts.fit(fit,cont.matrix)
  fit2 <- eBayes(fit2)
  tmp <- topTable(fit2, coef=comp, n=Inf)
  DEG_limma_voom <- na.omit(tmp)
  #head(DEG_limma_voom)
  result=cbind(rownames(DEG_limma_voom),DEG_limma_voom)
  colnames(result)[1]="SYMBOL"
  limma_result_all[[i]]=result;
  ind2=which(abs(DEG_limma_voom$logFC)> log2FCcut & DEG_limma_voom$adj.P.Val< adjPcut)
  message(paste("found",length(ind2), "DEG in" ,paste(com,collapse=" vs "), "by limma."))
  loginfo <- paste0("log2FC cutoff: ",log2FCcut, "; adjust Pvalue cutoff: ",adjPcut,
                   "\n found ",length(ind2), " DEG in" ,paste(com,collapse=" vs "), " by limma.\n")
  
  if(length(ind2) == 0) 
    limma_result[[i]]=NA else 
      limma_result[[i]]=result[ind2,]
  names(limma_result_all)[i] = paste0(groupA[i],"vs",groupB[i])
  names(limma_result)[i] = paste0(groupA[i],"vs",groupB[i])
  
}

##edgeR #######
#### 第一步，构建edgeR的DGEList对象，并过滤
library(edgeR)
edgeR_result=list()
edgeR_result_all=list()

for(i in seq_along(groupB)){
  ctl = which(group_list == groupA[i])
  ind=which(group_list==groupB[i])
  grp <- factor(group_list[c(ctl,ind)],levels = c(groupA[i],groupB[i]))
  design <- model.matrix(~0+grp)
  colnames(design) <- c(groupA[i],groupB[i])
  rownames(design) <- colnames(counts)[c(ctl,ind)]
  DEG <- DGEList(counts=counts[,c(ctl,ind)])
  keep <- filterByExpr(DEG,c(ctl,ind))
  table(keep)
  DEG <- DEG[keep ,keep.lib.sizes=FALSE]
  #DEG$samples$lib.size <- colSums(DEG$counts)
  DEG <- calcNormFactors(DEG)
  #### 第二步，差异表达分析
  DEG <- estimateGLMCommonDisp(DEG,design)
  DEG <- estimateGLMTrendedDisp(DEG, design)
  DEG <- estimateGLMTagwiseDisp(DEG, design)
  # 拟合线性模型
  fit <- glmFit(DEG, design)
  # 进行差异分析
  # 1,-1意味着前比后
  com=unique(group_list[c(ind,ctl)])
  comp=paste(com,collapse="-")
  cont.matrix <- makeContrasts(contrasts=c(comp),levels = design)
  lrt <- glmLRT(fit, contrast=cont.matrix) 
  #### 第三步，提取过滤差异分析结果
  edgeR_DEG <- topTags(lrt, n=nrow(DEG))
  edgeR_DEG <- as.data.frame(edgeR_DEG)
  result=cbind(rownames(edgeR_DEG),edgeR_DEG)
  colnames(result)[1]="SYMBOL"
  edgeR_result_all[[i]]=result
  ind2=which(abs(edgeR_DEG$logFC)> log2FCcut &edgeR_DEG$FDR < adjPcut);
  message(paste("found",length(ind2), "DEG in" ,paste(com,collapse=" vs "),"by edgeR."))
  loginfo <- paste(loginfo,
                   "found",length(ind2), "DEG in" ,paste(com,collapse=" vs "), "by edgeR.\n")
  
  if(length(ind2) == 0) 
    edgeR_result[[i]]=NA else
      edgeR_result[[i]]=result[ind2,]
  names(edgeR_result_all)[i] = paste0(groupA[i],"vs",groupB[i])
  names(edgeR_result)[i] = paste0(groupA[i],"vs",groupB[i])
  
}



##DESeq2 #########
suppressMessages(library(DESeq2))
#### 第一步，构建DESeq2的DESeq对象
counts_round=round(counts)
DESeq2_result=list()
DESeq2_result_all=list()
for(i in seq_along(groupB)){
  ctl = which(group_list == groupA[i])
  ind=which(group_list==groupB[i])
  com=unique(group_list[c(ind,ctl)])
  comp=paste(com,collapse="-")
  condition=factor(c(group_list[c(ctl,ind)]),levels = c(groupA[i],groupB[i]))
  colData <- data.frame(row.names=colnames(counts[,c(ctl,ind)]),group_list=condition)
  
  dds <- DESeqDataSetFromMatrix(countData = counts_round[,c(ctl,ind)],colData = colData,
                                design = ~group_list)
  #### 第二步，进行差异表达分析
  
  dds2 <- DESeq(dds)
  # 保存差异表达分析结果
  #### 第二步，提取差异分析结果
  res <- results(dds2)
  res <- res[order(res$padj),]
  DEG <- as.data.frame(res)
  result=cbind(rownames(DEG),DEG)
  colnames(result)[1]="SYMBOL";
  DESeq2_result_all[[i]]=result
  ind2=which(abs(DEG$log2FoldChange)>log2FCcut & DEG$padj < adjPcut & !is.na(DEG$padj)); #250918 old is ind
  message(paste("found",length(ind2), "DEG in" ,paste(com,collapse=" vs "),"by DESeq2."))
  loginfo <- paste(loginfo,
                   "found",length(ind2), "DEG in" ,paste(com,collapse=" vs "), "by DESeq2.\n")
  
  if(length(ind2) == 0) 
    DESeq2_result[[i]]=NA else 
      DESeq2_result[[i]]=result[ind2,]
  names(DESeq2_result_all)[i] = paste0(groupA[i],"vs",groupB[i])
  names(DESeq2_result)[i] = paste0(groupA[i],"vs",groupB[i])
}

## Function to create volcano plot for troika methods
create_volcano_plot_troika <- function(result, method, comp_name, log2FCcut, adjPcut, id_col = "SYMBOL") {
  # Load required libraries
  if(!requireNamespace("ggplot2", quietly = TRUE) || 
     !requireNamespace("ggrepel", quietly = TRUE)) {
    message("ggplot2 or ggrepel not available, skipping volcano plot.")
    return(NULL)
  }
  
  library(ggplot2)
  library(ggrepel)
  
  # Determine column names based on method
  if(method == "limma") {
    logfc_col <- "logFC"
    padj_col <- "adj.P.Val"
  } else if(method == "edgeR") {
    logfc_col <- "logFC"
    padj_col <- "FDR"
  } else if(method == "DESeq2") {
    logfc_col <- "log2FoldChange"
    padj_col <- "padj"
  } else {
    stop("Unknown method: ", method)
  }
  
  # Check if required columns exist
  if(!logfc_col %in% colnames(result)) {
    message(paste("Column", logfc_col, "not found in result for method", method))
    return(NULL)
  }
  if(!padj_col %in% colnames(result)) {
    message(paste("Column", padj_col, "not found in result for method", method))
    return(NULL)
  }
  
  # Clean data for plotting - remove NA/Inf values
  valid_rows <- !is.na(result[[logfc_col]]) & !is.na(result[[padj_col]]) &
                !is.infinite(result[[logfc_col]]) & !is.infinite(result[[padj_col]]) &
                result[[padj_col]] > 0
  result_clean <- result[valid_rows, ]
  
  if(nrow(result_clean) == 0) {
    message("No valid data for volcano plot after removing NA/Inf values")
    return(NULL)
  }
  
  # Calculate significance
  result_clean$significance <- "Not Sig"
  result_clean$significance[result_clean[[logfc_col]] > log2FCcut & result_clean[[padj_col]] < adjPcut] <- "Up"
  result_clean$significance[result_clean[[logfc_col]] < -log2FCcut & result_clean[[padj_col]] < adjPcut] <- "Down"
  result_clean$significance <- factor(result_clean$significance, levels = c("Up", "Down", "Not Sig"))
  
  # Calculate y-axis limit
  valid_pvals <- result_clean[[padj_col]][!is.na(result_clean[[padj_col]]) & result_clean[[padj_col]] > 0]
  if(length(valid_pvals) == 0) {
    y_limit <- 10
  } else {
    y_max <- max(-log10(valid_pvals), na.rm = TRUE)
    y_limit <- min(y_max * 1.1, 300)
  }
  
  # Create plot
  p <- ggplot(result_clean, aes(x = .data[[logfc_col]], y = -log10(.data[[padj_col]]), color = significance)) +
    geom_point(alpha = 0.6, size = 1.5) +
    scale_color_manual(values = c("Up" = "#E41A1C", "Down" = "#377EB8", "Not Sig" = "grey50"),
                       labels = c("Up" = paste0("Up (", sum(result_clean$significance == "Up", na.rm = TRUE), ")"),
                                  "Down" = paste0("Down (", sum(result_clean$significance == "Down", na.rm = TRUE), ")"),
                                  "Not Sig" = "Not Significant")) +
    geom_vline(xintercept = c(-log2FCcut, log2FCcut), linetype = "dashed", color = "grey40") +
    geom_hline(yintercept = -log10(adjPcut), linetype = "dashed", color = "grey40") +
    coord_cartesian(ylim = c(0, y_limit)) +
    labs(title = paste0("Volcano Plot: ", comp_name, " (", method, ")"),
         subtitle = paste0("log2FC cutoff: ", log2FCcut, ", adj.P cutoff: ", adjPcut),
         x = "Log2 Fold Change",
         y = "-Log10 Adjusted P-value",
         color = "Significance") +
    theme_bw() +
    theme(legend.position = "right",
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = 10))
  
  # Add labels for significant genes (top 20 by p-value)
  sig_genes <- result_clean[result_clean$significance %in% c("Up", "Down"), ]
  if(nrow(sig_genes) > 0) {
    sig_genes <- sig_genes[order(sig_genes[[padj_col]]), ]
    top_genes <- head(sig_genes, 20)
    
    if(id_col %in% colnames(top_genes) && any(!is.na(top_genes[[id_col]]))) {
      p <- p + geom_label_repel(
        data = top_genes,
        aes(label = .data[[id_col]]),
        size = 3,
        max.overlaps = 15,
        box.padding = 0.3,
        point.padding = 0.3,
        segment.size = 0.2,
        na.rm = TRUE
      )
    }
  }
  
  return(p)
}

# Function to save plots as both PDF and PNG
save_plot_both <- function(plot, filename_base, width = 8, height = 6, dpi = 300) {
  pdf_file <- file.path(plot_dir, paste0(filename_base, ".pdf"))
  png_file <- file.path(assets_dir, paste0(filename_base, ".png"))
  
  # Save as PDF
  pdf(pdf_file, width = width, height = height)
  print(plot)
  dev.off()
  
  # Save as PNG
  png(png_file, width = width * dpi, height = height * dpi, res = dpi)
  print(plot)
  dev.off()
  
  return(list(pdf = pdf_file, png = png_file))
}

# Generate volcano plots for all methods
message("\n--- Generating volcano plots ---")
volcano_plots <- list()

for(i in seq_along(groupB)) {
  comp_name <- paste0(groupA[i], "vs", groupB[i])
  safe_name <- gsub("[^a-zA-Z0-9]", "_", comp_name)
  
  # limma volcano
  volc_limma <- create_volcano_plot_troika(limma_result_all[[i]], "limma", comp_name, log2FCcut, adjPcut, idinf_col)
  if(!is.null(volc_limma)) {
    plot_files <- save_plot_both(volc_limma, paste0("step3_troika_volcano_limma_", safe_name), width = 10, height = 8)
    volcano_plots[[paste0(comp_name, "_limma")]] <- plot_files$png
    message(paste0("  Volcano plot (limma) saved: ", basename(plot_files$png)))
  }
  
  # edgeR volcano
  volc_edgeR <- create_volcano_plot_troika(edgeR_result_all[[i]], "edgeR", comp_name, log2FCcut, adjPcut, idinf_col)
  if(!is.null(volc_edgeR)) {
    plot_files <- save_plot_both(volc_edgeR, paste0("step3_troika_volcano_edgeR_", safe_name), width = 10, height = 8)
    volcano_plots[[paste0(comp_name, "_edgeR")]] <- plot_files$png
    message(paste0("  Volcano plot (edgeR) saved: ", basename(plot_files$png)))
  }
  
  # DESeq2 volcano
  volc_deseq <- create_volcano_plot_troika(DESeq2_result_all[[i]], "DESeq2", comp_name, log2FCcut, adjPcut, idinf_col)
  if(!is.null(volc_deseq)) {
    plot_files <- save_plot_both(volc_deseq, paste0("step3_troika_volcano_DESeq2_", safe_name), width = 10, height = 8)
    volcano_plots[[paste0(comp_name, "_DESeq2")]] <- plot_files$png
    message(paste0("  Volcano plot (DESeq2) saved: ", basename(plot_files$png)))
  }
}
message("--- Volcano plot generation completed ---\n")

## RRHO Analysis ########
rrho_results <- list()

if(rrho_available) {
  message("\n--- Performing RRHO analysis ---")
  library(RRHO)
  
  for(i in seq_along(groupB)) {
    comp_name <- names(limma_result_all)[i]
    safe_name <- gsub("[^a-zA-Z0-9]", "_", comp_name)
    message(paste0("  Processing comparison: ", comp_name))
    
    # Get data for each method
    limma_data <- limma_result_all[[i]]
    edgeR_data <- edgeR_result_all[[i]]
    deseq_data <- DESeq2_result_all[[i]]
    
    # Function to prepare RRHO data
    prepare_rrho_data <- function(data, logfc_col, padj_col) {
      # Create ranking based on -log10(padj) * sign(logFC)
      data$rank <- -log10(data[[padj_col]]) * sign(data[[logfc_col]])
      data <- data[order(data$rank, decreasing = TRUE), ]
      result <- data[, c("SYMBOL", "rank")]
      colnames(result) <- c("SYMBOL", "rank")
      return(result)
    }
    
    # Prepare data for each method
    limma_rrho <- prepare_rrho_data(limma_data, "logFC", "adj.P.Val")
    edgeR_rrho <- prepare_rrho_data(edgeR_data, "logFC", "FDR")
    deseq_rrho <- prepare_rrho_data(deseq_data, "log2FoldChange", "padj")
    
    # Get common genes
    common_genes <- Reduce(intersect, list(limma_rrho$SYMBOL, edgeR_rrho$SYMBOL, deseq_rrho$SYMBOL))
    message(paste0("    Common genes across all methods: ", length(common_genes)))
    
    if(length(common_genes) < 100) {
      message("    Too few common genes, skipping RRHO for this comparison.")
      next
    }
    
    # Filter to common genes
    limma_rrho <- limma_rrho[limma_rrho$SYMBOL %in% common_genes, ]
    edgeR_rrho <- edgeR_rrho[edgeR_rrho$SYMBOL %in% common_genes, ]
    deseq_rrho <- deseq_rrho[deseq_rrho$SYMBOL %in% common_genes, ]
    
    # Perform pairwise RRHO comparisons
    comparisons <- list(
      list(name1 = "limma", name2 = "edgeR", data1 = limma_rrho, data2 = edgeR_rrho),
      list(name1 = "limma", name2 = "DESeq2", data1 = limma_rrho, data2 = deseq_rrho),
      list(name1 = "edgeR", name2 = "DESeq2", data1 = edgeR_rrho, data2 = deseq_rrho)
    )
    
    for(comp in comparisons) {
      tryCatch({
        message(paste0("    RRHO: ", comp$name1, " vs ", comp$name2))
        
        # Run RRHO
        rrho_res <- RRHO(comp$data1, comp$data2, 
                         stepsize = max(1, floor(length(common_genes) / 100)),
                         labels = c(comp$name1, comp$name2),
                         log10.ind = TRUE)
        
        # Save result
        rrho_key <- paste0(comp_name, "_", comp$name1, "_vs_", comp$name2)
        rrho_results[[rrho_key]] <- rrho_res
        
        # Generate and save heatmap
        pdf_file <- file.path(plot_dir, paste0("step3_troika_rrho_", safe_name, "_", comp$name1, "_", comp$name2, ".pdf"))
        png_file <- file.path(assets_dir, paste0("step3_troika_rrho_", safe_name, "_", comp$name1, "_", comp$name2, ".png"))
        
        # Save PDF
        pdf(pdf_file, width = 8, height = 8)
        tryCatch({
          RRHO::heatmapRRHO(rrho_res, 
                           alternative = "two.sided",
                           color.gradient = c("blue", "white", "red"),
                           main = paste0("RRHO: ", comp$name1, " vs ", comp$name2, "\n", comp_name))
        }, error = function(e) {
          message(paste0("      Warning: RRHO heatmap error for ", comp$name1, " vs ", comp$name2, ": ", e$message))
        })
        dev.off()
        
        # Save PNG
        png(png_file, width = 2400, height = 2400, res = 300)
        tryCatch({
          RRHO::heatmapRRHO(rrho_res, 
                           alternative = "two.sided",
                           color.gradient = c("blue", "white", "red"),
                           main = paste0("RRHO: ", comp$name1, " vs ", comp$name2, "\n", comp_name))
        }, error = function(e) {
          message(paste0("      Warning: RRHO heatmap error for ", comp$name1, " vs ", comp$name2, ": ", e$message))
        })
        dev.off()
        
        message(paste0("      RRHO heatmap saved: ", basename(png_file)))
        
      }, error = function(e) {
        message(paste0("    Error in RRHO analysis for ", comp$name1, " vs ", comp$name2, ": ", e$message))
      })
    }
  }
  message("--- RRHO analysis completed ---\n")
} else {
  message("RRHO package not available, skipping RRHO analysis.")
}

##save ########
save(limma_result,edgeR_result,DESeq2_result,
     file = file.path(outputdir, paste0(output,"_DEG_troika.RData")))
save(limma_result_all,edgeR_result_all,DESeq2_result_all,rrho_results,
     file = file.path(outputdir, paste0(output,"_stat_troika.RData")))

loginfo <- paste(loginfo,
                 "   saved in",file.path(outputdir, paste0(output,"_stat_troika.RData")),"\n")
cat("************************************\n")
cat(loginfo)
cat("************************************\n")

## Generate Report Markdown ####
report_file <- file.path(outputdir, paste0(output, "_DEG_troika_report.md"))

# Initialize report
report_lines <- c()
report_lines <- c(report_lines, "# DEG Troika Analysis Report")
report_lines <- c(report_lines, "")

# Analysis Overview
report_lines <- c(report_lines, "## Analysis Overview")
report_lines <- c(report_lines, "")
report_lines <- c(report_lines, "| Item | Value |")
report_lines <- c(report_lines, "|------|-------|")
report_lines <- c(report_lines, sprintf("| Date | %s |", Sys.time()))
report_lines <- c(report_lines, sprintf("| Output Prefix | %s |", output))
report_lines <- c(report_lines, sprintf("| Output Directory | %s |", outputdir))
report_lines <- c(report_lines, sprintf("| Log2FC cutoff | %s |", log2FCcut))
report_lines <- c(report_lines, sprintf("| adjP cutoff | %s |", adjPcut))
report_lines <- c(report_lines, sprintf("| Group Variable | %s |", grpname))
report_lines <- c(report_lines, sprintf("| Group A | %s |", paste(groupA, collapse = ", ")))
report_lines <- c(report_lines, sprintf("| Group B | %s |", paste(groupB, collapse = ", ")))
report_lines <- c(report_lines, sprintf("| ID info column | %s |", idinf_col))
report_lines <- c(report_lines, sprintf("| RRHO Analysis | %s |", ifelse(rrho_available, "Performed", "Skipped (package not installed)")))
report_lines <- c(report_lines, "")

# Input Information
report_lines <- c(report_lines, "## Input Information")
report_lines <- c(report_lines, "")
report_lines <- c(report_lines, sprintf("| Sample count | %d |", ncol(counts)))
report_lines <- c(report_lines, sprintf("| Feature count | %d |", nrow(counts)))
report_lines <- c(report_lines, "")

# Group Information
report_lines <- c(report_lines, "## Group Information")
report_lines <- c(report_lines, "")
report_lines <- c(report_lines, "| Group | Sample Count |")
report_lines <- c(report_lines, "|-------|-------------|")
all_groups <- unique(c(groupA, groupB))
for(grp in all_groups) {
  grp_count <- length(which(group_list == grp))
  report_lines <- c(report_lines, sprintf("| %s | %d |", grp, grp_count))
}
report_lines <- c(report_lines, "")

# DEG Summary by Method
report_lines <- c(report_lines, "## DEG Summary by Method")
report_lines <- c(report_lines, "")
report_lines <- c(report_lines, "| Comparison | Method | Total DEGs | Up | Down |")
report_lines <- c(report_lines, "|------------|--------|------------|-----|------|")

for(i in seq_along(groupB)) {
  comp_name <- paste0(groupA[i], "vs", groupB[i])
  
  # limma stats
  limma_data <- limma_result_all[[i]]
  limma_up <- sum(limma_data$logFC > log2FCcut & limma_data$adj.P.Val < adjPcut, na.rm = TRUE)
  limma_down <- sum(limma_data$logFC < -log2FCcut & limma_data$adj.P.Val < adjPcut, na.rm = TRUE)
  limma_total <- limma_up + limma_down
  report_lines <- c(report_lines, sprintf("| %s | limma | %d | %d | %d |", comp_name, limma_total, limma_up, limma_down))
  
  # edgeR stats
  edgeR_data <- edgeR_result_all[[i]]
  edgeR_up <- sum(edgeR_data$logFC > log2FCcut & edgeR_data$FDR < adjPcut, na.rm = TRUE)
  edgeR_down <- sum(edgeR_data$logFC < -log2FCcut & edgeR_data$FDR < adjPcut, na.rm = TRUE)
  edgeR_total <- edgeR_up + edgeR_down
  report_lines <- c(report_lines, sprintf("| %s | edgeR | %d | %d | %d |", comp_name, edgeR_total, edgeR_up, edgeR_down))
  
  # DESeq2 stats
  deseq_data <- DESeq2_result_all[[i]]
  deseq_up <- sum(deseq_data$log2FoldChange > log2FCcut & deseq_data$padj < adjPcut & !is.na(deseq_data$padj), na.rm = TRUE)
  deseq_down <- sum(deseq_data$log2FoldChange < -log2FCcut & deseq_data$padj < adjPcut & !is.na(deseq_data$padj), na.rm = TRUE)
  deseq_total <- deseq_up + deseq_down
  report_lines <- c(report_lines, sprintf("| %s | DESeq2 | %d | %d | %d |", comp_name, deseq_total, deseq_up, deseq_down))
}
report_lines <- c(report_lines, "")

# RRHO Analysis section
if(rrho_available && length(rrho_results) > 0) {
  report_lines <- c(report_lines, "## RRHO Analysis")
  report_lines <- c(report_lines, "")
  report_lines <- c(report_lines, "RRHO (Rank Rank Hypergeometric Overlap) analysis was performed to assess the consistency between different DEG methods.")
  report_lines <- c(report_lines, "")
  report_lines <- c(report_lines, "### RRHO Comparisons")
  report_lines <- c(report_lines, "")
  report_lines <- c(report_lines, "| Comparison | Methods Compared |")
  report_lines <- c(report_lines, "|------------|------------------|")
  
  for(rrho_key in names(rrho_results)) {
    parts <- strsplit(rrho_key, "_")[[1]]
    comp <- paste(parts[1:(length(parts)-3)], collapse = "_")
    method1 <- parts[length(parts)-2]
    method2 <- parts[length(parts)]
    report_lines <- c(report_lines, sprintf("| %s | %s vs %s |", comp, method1, method2))
  }
  report_lines <- c(report_lines, "")
  
  # RRHO Heatmaps
  report_lines <- c(report_lines, "### RRHO Heatmaps")
  report_lines <- c(report_lines, "")
  report_lines <- c(report_lines, "Heatmaps show the rank-rank hypergeometric overlap between methods. Red indicates agreement, blue indicates disagreement.")
  report_lines <- c(report_lines, "")
  
  for(rrho_key in names(rrho_results)) {
    parts <- strsplit(rrho_key, "_")[[1]]
    comp <- paste(parts[1:(length(parts)-3)], collapse = "_")
    safe_comp <- gsub("[^a-zA-Z0-9]", "_", comp)
    method1 <- parts[length(parts)-2]
    method2 <- parts[length(parts)]
    
    report_lines <- c(report_lines, sprintf("**%s: %s vs %s**", comp, method1, method2))
    report_lines <- c(report_lines, "")
    report_lines <- c(report_lines, sprintf("![RRHO](assets/step3_troika_rrho_%s_%s_%s.png)", safe_comp, method1, method2))
    report_lines <- c(report_lines, "")
  }
} else if(rrho_available) {
  report_lines <- c(report_lines, "## RRHO Analysis")
  report_lines <- c(report_lines, "")
  report_lines <- c(report_lines, "RRHO analysis was attempted but no valid results were generated (possibly due to insufficient common genes).")
  report_lines <- c(report_lines, "")
} else {
  report_lines <- c(report_lines, "## RRHO Analysis")
  report_lines <- c(report_lines, "")
  report_lines <- c(report_lines, "RRHO analysis was skipped because the RRHO package is not installed.")
  report_lines <- c(report_lines, "")
  report_lines <- c(report_lines, "To install RRHO, run: `BiocManager::install('RRHO')`")
  report_lines <- c(report_lines, "")
}

# Visualization section
report_lines <- c(report_lines, "## Visualization")
report_lines <- c(report_lines, "")

# Volcano Plots
report_lines <- c(report_lines, "### Volcano Plots")
report_lines <- c(report_lines, "")

for(i in seq_along(groupB)) {
  comp_name <- paste0(groupA[i], "vs", groupB[i])
  safe_name <- gsub("[^a-zA-Z0-9]", "_", comp_name)
  
  report_lines <- c(report_lines, sprintf("**%s**", comp_name))
  report_lines <- c(report_lines, "")
  
  report_lines <- c(report_lines, "| Method | Plot |")
  report_lines <- c(report_lines, "|--------|------|")
  report_lines <- c(report_lines, sprintf("| limma | ![Volcano](assets/step3_troika_volcano_limma_%s.png) |", safe_name))
  report_lines <- c(report_lines, sprintf("| edgeR | ![Volcano](assets/step3_troika_volcano_edgeR_%s.png) |", safe_name))
  report_lines <- c(report_lines, sprintf("| DESeq2 | ![Volcano](assets/step3_troika_volcano_DESeq2_%s.png) |", safe_name))
  report_lines <- c(report_lines, "")
}

# Output Files
report_lines <- c(report_lines, "## Output Files")
report_lines <- c(report_lines, "")
report_lines <- c(report_lines, "| File | Description |")
report_lines <- c(report_lines, "|------|-------------|")
report_lines <- c(report_lines, sprintf("| `%s_DEG_troika.RData` | Filtered DEG results (significant genes only) |", output))
report_lines <- c(report_lines, sprintf("| `%s_stat_troika.RData` | Complete statistical results (all genes) |", output))
report_lines <- c(report_lines, sprintf("| `plot/step3_troika_volcano_*.pdf` | Volcano plots (PDF format) |"))
report_lines <- c(report_lines, sprintf("| `assets/step3_troika_volcano_*.png` | Volcano plots (PNG format for reports) |"))
if(rrho_available && length(rrho_results) > 0) {
  report_lines <- c(report_lines, sprintf("| `plot/step3_troika_rrho_*.pdf` | RRHO heatmaps (PDF format) |"))
  report_lines <- c(report_lines, sprintf("| `assets/step3_troika_rrho_*.png` | RRHO heatmaps (PNG format for reports) |"))
}
report_lines <- c(report_lines, "")

# Processing Log
report_lines <- c(report_lines, "## Processing Log")
report_lines <- c(report_lines, "")
report_lines <- c(report_lines, "```")
report_lines <- c(report_lines, loginfo)
report_lines <- c(report_lines, "```")
report_lines <- c(report_lines, "")

# Footer
report_lines <- c(report_lines, "---")
report_lines <- c(report_lines, "")
report_lines <- c(report_lines, paste0("*Analysis completed at:* ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
report_lines <- c(report_lines, "")
report_lines <- c(report_lines, "*Generated by DEG_troika v2.1*")

# Write report
writeLines(report_lines, report_file)
message("\nMarkdown report saved to: ", report_file)
