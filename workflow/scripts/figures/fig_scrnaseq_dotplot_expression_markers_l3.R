#!/usr/bin/env Rscript

# The goal of this script is to create a dotplot of the canonical markers for all celltypes (l3).

library(Seurat)
library(dplyr)
library(ggplot2)
library(ggrastr)
library(ggrepel)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop(paste0("Script needs 3 arguments. Current input is:", args))
}

seurat_rds <- args[1]
manual_l3_order_xlsx <- args[2]
marker_genes <- c("PTPRC", "MKI67", "CD3D", "CD4", "CD8A", "SELL", "CCR5", "S100A4", "FOXP3", "TIGIT", "CTLA4", "TRDC", "TRGC1", "GNLY", "KLRB1", "NCAM1", "NKG7", "GZMA", "GZMB", "GZMK", "GZMH", "IL7R", "CD2", "BANK1", "MS4A1", "CD27", "IGHD", "IGHM", "JCHAIN", "HLA-DRA", "CST3", "CD14", "FCGR3A", "CD1C", "FCER1A", "SIGLEC6", "AXL", "IRF4", "IRF8", "CLEC4C", "KLF4", "PPBP", "CD34")
dotplot_pdf <- args[3]

seuratObject <- readRDS(seurat_rds)

manual_l3_order <- readxl::read_excel(manual_l3_order_xlsx, col_names = T) %>%
  pull(Celltype)

marker_expr <- GetAssayData(seuratObject, assay = "RNA")[which(rownames(GetAssayData(seuratObject, assay = "RNA")) %in% marker_genes), ]

marker_expr_median <- data.frame(GeneID = rownames(marker_expr), marker_expr) %>%
  tidyr::pivot_longer(-GeneID, names_to = "CB", values_to = "nUMIs") %>%
  dplyr::mutate(CB = stringr::str_replace(CB, "\\.", "-")) %>%
  dplyr::inner_join(data.frame(CB = rownames(seuratObject@meta.data), 
                               Celltype = seuratObject@meta.data$manual_l3), 
                    by = "CB") %>%
  dplyr::mutate(Celltype = factor(Celltype), 
                GeneID = factor(GeneID)) %>%
  dplyr::group_by(Celltype, GeneID, .drop = F) %>%
  dplyr::summarize(Median = median(log1p(nUMIs)),
                   Percentage = mean(nUMIs>0)*100) %>%
  dplyr::mutate(Median = ifelse(is.na(Median), 0, Median),
                Percentage = ifelse(is.na(Percentage), 0, Percentage),
                Celltype = factor(Celltype, manual_l3_order),
                GeneID = factor(GeneID, marker_genes))

dotplot_ggplot2 <- marker_expr_median %>% 
  ggplot(aes(x = GeneID, y = Celltype, col = Median)) +
  geom_tile(alpha = 0, col = "grey") +
  geom_point(aes(size = Percentage)) +
  theme_bw() +
  theme(legend.pos = "bottom", 
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.spacing.x=unit(0, "lines"))


pdf(width = 10.5, height = 7.5, file = dotplot_pdf)
print(dotplot_ggplot2)
dev.off()

sessionInfo()
