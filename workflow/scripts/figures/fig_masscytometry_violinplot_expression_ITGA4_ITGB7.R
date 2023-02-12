#!/usr/bin/env Rscript

# The goal of this script is to create a boxplot of l1 relative to all PBMCs grouped by response and facetted by lineage.

library(Seurat)
library(dplyr)
library(ggplot2)
library(ggrastr)
library(ggrepel)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop(paste0("Script needs 2 arguments. Current input is:", args))
}

seurat_rds_path <- args[1]
manual_l3_order_xlsx <- args[2]
violinplot_pdf <- args[3]

manual_l3_order <- readxl::read_excel(manual_l3_order_xlsx, col_names = F) %>%
  pull(`...1`)

marker_expr <- GetAssayData(seuratObject, assay = "RNA")[which(rownames(GetAssayData(seuratObject, assay = "RNA")) %in% c("ITGA4", "ITGB7")), ]

marker_expr <- data.frame(GeneID = rownames(marker_expr), marker_expr) %>%
  tidyr::pivot_longer(-GeneID, names_to = "CB", values_to = "nUMIs") %>%
  dplyr::mutate(CB = stringr::str_replace(CB, "\\.", "-")) %>%
  dplyr::inner_join(data.frame(CB = rownames(seuratObject@meta.data), 
                               Celltype = seuratObject@meta.data$manual_l3), 
                    by = "CB") %>%
  dplyr::mutate(Celltype = factor(Celltype), 
                GeneID = factor(GeneID),
                Celltype = factor(Celltype, manual_l3_order))

plotobj <- marker_expr %>% 
  ggplot(aes(x = Celltype, y = nUMIs)) +
  geom_jitter_rast(alpha = 0.25) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  #geom_violin() +
  facet_wrap(~GeneID, nrow = 2, ncol = 1) +
  theme_bw() +
  theme(legend.pos = "bottom", 
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())


pdf(width = 10, height = 5, file = dotplot_pdf)
print(plotobj)
dev.off()

png(width = 10, height = 5, units = "in", res = 120, file = dotplot_png)
print(plotobj)
dev.off()
