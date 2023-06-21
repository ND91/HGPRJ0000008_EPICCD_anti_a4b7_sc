#!/usr/bin/env bash

# The goal of this script is to perform pseudobulk differential expression analysis on the different celltypes.

library(Seurat)
library(dplyr)
library(DESeq2)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop(paste0("Script needs 4 arguments. Current input is:", args))
}

seurat_rds <- args[1]
level <- args[2]
functions_r <- args[3]
degs_list_rds <- args[4]

source(functions_r)

seuratObject <- readRDS(seurat_rds)

pb_sample_metadata <- seuratObject@meta.data %>%
  dplyr::select(SampleID, Response, Sex, Age) %>%
  unique() %>%
  dplyr::mutate(Response_coded = ifelse(Response == "Responder", "R", "NR"))
rownames(pb_sample_metadata) <- pb_sample_metadata$SampleID

if(level != "manual_l0"){
  degs_rvnr <- seuratDE(seuratobj = seuratObject, 
                        cellsampleID = "SampleID", 
                        cellclusterID = level, 
                        sampleinfo = pb_sample_metadata, 
                        design = "~Response_coded+Sex+Age", 
                        contrast = c("Response_coded", "R", "NR"))
} else{
  degs_rvnr <- seuratDE(seuratobj = seuratObject, 
                        cellsampleID = "SampleID", 
                        sampleinfo = pb_sample_metadata, 
                        design = "~Response_coded+Sex+Age", 
                        contrast = c("Response_coded", "R", "NR"))
}

degs_rvnr <- degs_rvnr[lapply(degs_rvnr, length) != 0]

saveRDS(degs_rvnr, degs_list_rds, compress = "gzip")

sessionInfo()