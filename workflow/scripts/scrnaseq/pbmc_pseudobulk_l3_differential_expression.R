#!/usr/bin/env bash

# The goal of this script is to perform pseudobulk differential expression analysis on the different celltypes.

library(Seurat)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop(paste0("Script needs 2 arguments. Current input is:", args))
}

seurat_rds <- args[1]
degs_list_rds <- args[2]

seuratDE <- function(seuratobj, cellsampleID, cellclusterID = NULL, sampleinfo, design, contrast = NULL, name = NULL){
  # seuratDE is a function to perform pseudobulk differential expression on all celltypes in a column listed in one of the SeuratObject columns. It essentially creates a pseudobulk dataset and wraps around DESeq2 to perform DE analyses.
  #   input:
  #     seuratobj: The SeuratObject of interest.
  #     sampleinfo: A dataframe containing the sample metadata, which will be fed to DESeq2 as the colData.  Note that this is not the same as the cell metadata (albeit a reduced form thereof).
  #     cellsampleID: Column name in the SeuratObject metadata containing the sample identifier.  
  #     cellclusterID: Column name in the SeuratObject by which the cells should be aggregated. If left empty, function will perform pseudobulk DE on all cells.
  #     design: The design formula to be used in the analysis. Note that the columns should be present in \"sampleinfo\".
  #     contrast: The comparison of interest per the DESeq2 requirements. Either \"name\" or \"contrast\" should be provided.
  #     name: The name of the result you want from the DESeq2 output. Either \"name\" or \"contrast\" should be provided.
  #   output: A list of celltypes each of which includes a DESeq2 object and the results per the contrasts indicated previously ranked by p-value.
  
  require(DESeq2)
  
  if(!class(seuratobj) != "SeuratObj") stop("\"seuratobj\" must be of the class SeuratObject")
  if(!cellsampleID %in% colnames(seuratobj@meta.data)) stop("\"cellsampleID\" cannot be found among the cell metadata column names in the provided \"seuratobj\"")
  #if(!cellclusterID %in% colnames(seuratobj@meta.data)) warning("\"cellclusterID\" cannot be found among the cell metadata column names in the provided \"seuratobj\", will pool all cells together")
  if(is.null(contrast) & is.null(name)) stop("Either \"name\" or \"contrast\" should be provided, not both.")
  
  if(!is.null(cellclusterID)){
    seuratobj_list <- Seurat::SplitObject(object = seuratobj, split.by = cellclusterID)
  } else{
    seuratobj_list <- list(seuratobj)
  }
  
  de_list <- lapply(seuratobj_list, function(seuratentry){
    
    tryCatch(expr = {
      counts_mat <- GetAssayData(seuratentry, slot = "counts")
      
      sourceID <- paste0("~0+", cellsampleID)
      
      counts_sample <- counts_mat %*% model.matrix(as.formula(sourceID), data = seuratentry@meta.data)
      colnames(counts_sample) <- gsub(cellsampleID, "", colnames(counts_sample))
      
      counts_sample <- counts_sample[which(Matrix::rowSums(counts_sample) != 0),]
      
      dds <- DESeqDataSetFromMatrix(countData = counts_sample,
                                    colData = sampleinfo[colnames(counts_sample),],
                                    design = as.formula(design))
      dds <- DESeq(dds)
      
      if(is.null(name)){
        degs <- results(dds, contrast = contrast, independentFiltering = T)
      } else{
        degs <- results(dds, name = name, independentFiltering = T)
      }
      
      degs <- degs[order(degs$pvalue),]
      degs <- degs[!is.na(degs$padj),]
      degs$gene <- rownames(degs)
      return(list(degs = degs,
                  dds = dds))
    },
    error = function(cond) {
      message("Failed with message:")
      message(cond)
      
      return(NULL)
    })
  })
  return(de_list)
}

seuratObject <- readRDS(seurat_rds)

pb_sample_metadata <- seuratObject@meta.data %>%
  dplyr::select(SampleID, Response, Sex, Age) %>%
  unique() %>%
  dplyr::mutate(Response_coded = ifelse(Response == "Responder", "R", "NR"))
rownames(pb_sample_metadata) <- pb_sample_metadata$SampleID

degs_rvnr <- seuratDE(seuratobj = seuratObject, 
                      cellsampleID = "SampleID", 
                      #cellclusterID = "manual_l3", 
                      sampleinfo = pb_sample_metadata, 
                      design = "~Response_coded+Sex+Age", 
                      contrast = c("Response_coded", "R", "NR"))

degs_rvnr <- degs_rvnr[lapply(degs_rvnr, length) != 0]

saveRDS(degs_rvnr, degs_list_rds, compress = "gzip")

sessionInfo()