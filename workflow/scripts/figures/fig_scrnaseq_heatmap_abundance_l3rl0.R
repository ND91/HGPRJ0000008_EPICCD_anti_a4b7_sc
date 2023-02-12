#!/usr/bin/env Rscript

# The goal of this script is to create a boxplot of l3 relative to all PBMCs grouped by response and facetted by lineage.

library(Seurat)
library(dplyr)
library(ComplexHeatmap)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop(paste0("Script needs 2 arguments. Current input is:", args))
}

seurat_rds_path <- args[1]
heatmap_pdf <- args[2]

abundance_l3rl0_long <- seuratObject@meta.data %>%
  dplyr::mutate(manual_l3 = as.factor(manual_l3)) %>%
  dplyr::group_by(SampleID, Response, manual_l3, .drop = F) %>%
  dplyr::summarize(Nl3sample = n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(SampleID) %>%
  dplyr::mutate(Ncells = sum(Nl3sample),
                Ncellprop = Nl3sample/Ncells,
                Ncellprop = ifelse(is.na(Ncellprop), 0, Ncellprop))

abundance_l3rl0_wide <- abundance_l3rl0_long %>%
  dplyr::select(-c(Response, Nl3sample, Ncells)) %>%
  tidyr::pivot_wider(names_from = SampleID, values_from = Ncellprop) %>%
  data.frame(., row.names = .$manual_l3) %>%
  dplyr::select(-manual_l3)

abundance_l3rl0_wide_scaled <- t(apply(abundance_l3rl0_wide*100, 1, scale))
colnames(abundance_l3rl0_wide_scaled) <- colnames(abundance_l3rl0_wide)

heatmap_colannotation <- abundance_l3rl0_long %>%
  dplyr::select(SampleID, Response) %>%
  unique() %>%
  data.frame(., row.names = .$SampleID) %>%
  select(-SampleID) %>%
  dplyr::arrange(Response)

heatmap_rowannotation <- seuratObject@meta.data %>%
  dplyr::select(manual_l1, manual_l3) %>%
  dplyr::filter(manual_l3 %in% rownames(abundance_l3rl0_wide_scaled)) %>%
  unique() %>%
  dplyr::arrange(manual_l1, manual_l3) 

# Number of cells

plotobj <- Heatmap(abundance_l3rl0_wide_scaled[heatmap_rowannotation$manual_l3, rownames(heatmap_colannotation)], 
                   name = "Abundance relative to all PBMCs",
                   show_row_dend = F,
                   split = heatmap_rowannotation$manual_l1,
                   cluster_columns = F,
                   col = viridis(100),
                   #rect_gp = gpar(col = "black", lwd = 1),
                   top_annotation = HeatmapAnnotation(Response = heatmap_colannotation$Response,
                                                      gp = gpar(col = "black"),
                                                      annotation_name_side = "left",
                                                      col = list(Response = c("Responder" = "#00BFC4", "Non-responder" = "#F8766D"))))

pdf(width = 8.5, height = 10, file = heatmap_pdf)
print(plotobj)
dev.off()

png(width = 8.5, height = 10, units = "in", res = 120, file = heatmap_png)
print(plotobj)
dev.off()
