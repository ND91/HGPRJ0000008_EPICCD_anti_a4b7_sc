#!/usr/bin/env bash

# The goal of this script is to perform differential abundance analysis on the different celltypes relative to the total population of PBMCs.

library(Seurat)
library(speckle)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop(paste0("Script needs 2 arguments. Current input is:", args))
}

sce_clustered_rds <- args[1]
dacs_csv <- args[2]

sce_clustered <- readRDS(sce_clustered_rds)

dacs <- speckle::propeller(clusters = colData(sce_clustered)$manual_l3, 
                           sample = colData(sce_clustered)$SampleID, 
                           group = colData(sce_clustered)$Response)

write.csv(dacs, dacs_csv)

sessionInfo()