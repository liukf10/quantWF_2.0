rm(list = ls())
setwd("D:/func/quantWF_2.0")
x <- list.files()
needtofix <- x[grepl("\\.R$",x)]
needtofix <- needtofix[!needtofix %in% c("preprocess_func.R","pvca.R","wgcna_plot.R","wgcna_tuning.R",
                                         "step4_GO&KEGGenrich.R","step4_GS_GSenrich.R","DEG_troika.R")]
if(!is.null(needtofix)){
  changed = gsub("\\.R","",needtofix)
  needdel = changed[changed %in% x]
  for(i in needdel) file.remove(i)
  for(i in needtofix) file.rename(i,gsub("\\.R","",i))
}

