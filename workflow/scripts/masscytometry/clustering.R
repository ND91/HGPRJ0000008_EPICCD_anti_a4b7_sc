#!/usr/bin/env Rscript
# The goal of this script is to perform UMAP dimension reduction on the masscytometry data.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop(paste0("Script needs 2 arguments. Current input is:", args))
}

library(SingleCellExperiment)
library(scater)
library(dplyr)
library(uwot)

sce_annotated_rds <- args[1]
sce_clustered_rds <- args[2]

sce_annotated <- readRDS(sce_annotated_rds)

set.seed(796318)
sce_clustered <- sce_annotated
reducedDims(sce_clustered) <- list(UMAP = DataFrame(scater::calculateUMAP(sce_annotated, exprs_values = "ncounts")))

saveRDS(sce_clustered, sce_clustered_rds)

sessionInfo()
