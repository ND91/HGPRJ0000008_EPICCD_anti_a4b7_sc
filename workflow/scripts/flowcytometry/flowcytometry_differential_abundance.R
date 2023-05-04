#!/usr/bin/env Rscript

# The goal of this script is to perform differential abundance testing on the celltypes measured through flow cytometry.

library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop(paste0("Script needs 2 arguments. Current input is:", args))
}

celltype_percentage_csv <- args[1]
dacs_csv <- args[2]

celltype_percentage <- read.csv(celltype_percentage_csv)

dacs <- celltype_percentage %>%
  dplyr::group_by(Celltype) %>% 
  do(w = wilcox.test(Percentage_relative_PBMCs~Response, data=., paired=FALSE)) %>% 
  summarise(Celltype, pvalue_wilcox = w$p.value)
  
write.csv(dacs, dacs_csv)

sessionInfo()