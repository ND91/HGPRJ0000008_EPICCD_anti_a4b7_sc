#!/usr/bin/env bash

# The goal of this script is to perform differential abundance analysis on the different celltypes relative to the total population of PBMCs.

library(Seurat)
library(speckle)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop(paste0("Script needs 2 arguments. Current input is:", args))
}

seurat_rds <- args[1]
dacs_csv <- args[2]

seuratObject <- readRDS(seurat_rds)

dacs <- propeller(clusters = seuratObject@meta.data$manual_l3, 
                  sample = seuratObject@meta.data$SampleID, 
                  group = seuratObject@meta.data$Response)

write.csv(dacs, dacs_csv)

sessionInfo()