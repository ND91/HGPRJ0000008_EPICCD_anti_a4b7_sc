#!/usr/bin/env bash

# The goal of this script is to perform pseudobulk differential expression analysis on the different celltypes at l3 level.

library(Seurat)
library(dplyr)
library(DESeq2)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop(paste0("Script needs 3 arguments. Current input is:", args))
}

seurat_rds <- args[1]
functions_r <- args[2]
degs_list_rds <- args[3]

source(functions_r)

seuratObject <- readRDS(seurat_rds)

pb_sample_metadata <- seuratObject@meta.data %>%
  dplyr::select(sampleID, status) %>%
  unique()
rownames(pb_sample_metadata) <- pb_sample_metadata$sampleID

degs_ivni <- seuratDE(seuratobj = seuratObject, 
                      cellsampleID = "sampleID", 
                      cellclusterID = "manual_l3", 
                      sampleinfo = pb_sample_metadata, 
                      design = "~status", 
                      contrast = c("status", "Involved", "Uninvolved"))

saveRDS(degs_ivni, degs_list_rds, compress = "gzip")

sessionInfo()