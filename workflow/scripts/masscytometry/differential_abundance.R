#!/usr/bin/env bash

# The goal of this script is to perform differential abundance analysis on the different celltypes relative to the total population of PBMCs.

devtools::install_github("phipsonlab/speckle")

library(speckle)
library(dplyr)
library(SingleCellExperiment)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop(paste0("Script needs 2 arguments. Current input is:", args))
}

sce_clustered_rds <- args[1]
dacs_csv <- args[2]

sce_clustered <- readRDS(sce_clustered_rds)

overlapping_samples <- data.frame(colData(sce_clustered)) %>%
  dplyr::filter(PBMC_scrnaseq == "Yes")

sce_clustered_overlapping <- sce_clustered[,rownames(overlapping_samples)]

dacs <- speckle::propeller(clusters = colData(sce_clustered_overlapping)$manual_l3, 
                           sample = colData(sce_clustered_overlapping)$Sample_ID, 
                           group = colData(sce_clustered_overlapping)$Response)

write.csv(dacs, dacs_csv)

sessionInfo()