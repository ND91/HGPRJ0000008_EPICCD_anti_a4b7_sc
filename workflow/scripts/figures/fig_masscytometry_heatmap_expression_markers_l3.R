#!/usr/bin/env Rscript

# The goal of this script is to create a heatmap of the proteins for all celltypes (l3).

library(SingleCellExperiment)
library(dplyr)
library(ComplexHeatmap)
library(viridisLite)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop(paste0("Script needs 3 arguments. Current input is:", args))
}

sce_rds <- args[1]
manual_l3_order_xlsx <- args[2]
marker_proteins <- c("CD45", "ITGA4", "ITGAL", "CD3", "CD4", "HAVR2", "CD44", "IL2RA", "CD8A", "CCR7", "CD7", "IL7RA", "PDCD1", "CTLA4", "CD28", "CD69", "CD45RA", "CD45RO", "CXCR3", "CCR4", "CCR9", "CD5", "CCR5", "TNR4", "TNR6", "TNR9", "CD9", "CD2", "CD57", "KLRB1", "CD16", "CCR10", "CD27", "HLADR", "CES1", "CD14")
heatmap_pdf <- args[3]

sce <- readRDS(sce_rds)

manual_l3_order <- readxl::read_excel(manual_l3_order_xlsx) %>%
  dplyr::mutate(number = as.numeric(factor(Celltype, levels = Celltype)),
                celltype_number = paste0(number, ". ", Celltype),
                celltype_number = factor(celltype_number, levels = celltype_number))

marker_expr <- assay(sce)[marker_proteins,]

marker_expr_median <- data.frame(protein = rownames(marker_expr), as.matrix(marker_expr)) %>%
  tidyr::pivot_longer(-protein, names_to = "EventID", values_to = "expression") %>%
  dplyr::inner_join(data.frame(EventID = colnames(sce), 
                               Celltype = colData(sce)$manual_l3), 
                    by = "EventID") %>%
  dplyr::mutate(Celltype = factor(Celltype), 
                protein = factor(protein)) %>%
  dplyr::group_by(Celltype, protein, .drop = F) %>%
  dplyr::summarize(Median = mean(expression),
                   Percentage = mean(expression>0)*100) %>%
  dplyr::mutate(Median = ifelse(is.na(Median), 0, Median),
                Percentage = ifelse(is.na(Percentage), 0, Percentage),
                Celltype = factor(Celltype, manual_l3_order$Celltype),
                protein = factor(protein, marker_proteins))

marker_expr_median_wide <- marker_expr_median %>%
  tidyr::pivot_wider(-Percentage, names_from = Celltype, values_from = Median)
marker_expr_median_wide_df <- data.frame(marker_expr_median_wide[,-1], row.names = marker_expr_median_wide$protein)
colnames(marker_expr_median_wide_df) <- colnames(marker_expr_median_wide)[-1]

overlapping_celltypes <- manual_l3_order$Celltype[manual_l3_order$Celltype %in% colnames(marker_expr_median_wide_df)]

marker_expr_median_wide_df <- marker_expr_median_wide_df[,overlapping_celltypes]

heatmap_complexheatmap <- ComplexHeatmap::pheatmap(as.matrix(marker_expr_median_wide_df), 
                                                   scale = "row",
                                                   color = viridis(1000),
                                                   cluster_cols = F)

pdf(width = 7.5, height = 10, file = heatmap_pdf)
print(heatmap_complexheatmap)
dev.off()

sessionInfo()
