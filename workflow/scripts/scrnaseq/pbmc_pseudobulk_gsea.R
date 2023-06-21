#!/usr/bin/env Rscript

# The goal of this script is to perform gene set enrichment analyses against the KEGG genesets.

library(DESeq2)
library(dplyr)
library(fgsea)
library(msigdbr)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop(paste0("Script needs 2 arguments. Current input is:", args))
}

degs_list_rds <- args[1]
fgsea_list_rds <- args[2]

#gs_hs <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME")
gs_hs <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG")
#gs_hs <- msigdbr(species = "Homo sapiens", category = "H")

gs_hs_list <- split(x = gs_hs$gene_symbol, f = gs_hs$gs_name)

degs_list <- readRDS(degs_list_rds)

fgsea_list <- lapply(degs_list, function(deg_object){
  degs_df <- deg_object$degs %>%
    data.frame() %>%
    dplyr::filter(!is.na(gene))
  
  waldstats <- degs_df$stat
  names(waldstats) <- degs_df$gene
  
  fgsea_obj <- fgsea(pathways = gs_hs_list, 
                     stats    = waldstats,
                     minSize  = 15,
                     maxSize  = 500)
  
  fgsea_obj <- fgsea_obj[order(fgsea_obj$pval),]
  
  return(fgsea_obj)
})

saveRDS(fgsea_list, fgsea_list_rds)

sessionInfo()