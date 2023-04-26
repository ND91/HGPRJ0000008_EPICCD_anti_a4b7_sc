#!/usr/bin/env Rscript

# The goal of this script is to create a scatterplot of the first two UMAP dimensions of all PBMCs colored by APC markers CD86 and HLA-DRA, as well as classical monocyte marker CD14 and non-classical monocyte marker FCGR3A.

library(Seurat)
library(dplyr)
library(ggplot2)
library(ggrastr)
library(ggrepel)
library(viridis)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop(paste0("Script needs 4 arguments. Current input is:", args))
}

seurat_rds <- args[1]
manual_l3_order_xlsx <- args[2]
umap_pdf <- args[3]
boxplot_pdf <- args[4]

seuratObject <- readRDS(seurat_rds)
manual_l3_order <- readxl::read_excel(manual_l3_order_xlsx) %>%
  dplyr::mutate(number = as.numeric(factor(Celltype, levels = Celltype)),
                celltype_number = paste0(number, ". ", Celltype),
                celltype_number = factor(celltype_number, levels = celltype_number))

expr_df <- data.frame(CB = colnames(seuratObject),
                      Embeddings(seuratObject[["umap"]]),
                      expr = GetAssayData(seuratObject)[c("CD86"),],
                      gene = "CD86",
                      seuratObject@meta.data) %>%
  dplyr::bind_rows(data.frame(CB = colnames(seuratObject),
                              Embeddings(seuratObject[["umap"]]),
                              expr = GetAssayData(seuratObject)["HLA-DRA",],
                              gene = "HLA-DRA",
                              seuratObject@meta.data)) %>%
  dplyr::bind_rows(data.frame(CB = colnames(seuratObject),
                              Embeddings(seuratObject[["umap"]]),
                              expr = GetAssayData(seuratObject)["CD14",],
                              gene = "CD14",
                              seuratObject@meta.data)) %>%
  dplyr::bind_rows(data.frame(CB = colnames(seuratObject),
                              Embeddings(seuratObject[["umap"]]),
                              expr = GetAssayData(seuratObject)["FCGR3A",],
                              gene = "FCGR3A",
                              seuratObject@meta.data)) %>%
  dplyr::left_join(manual_l3_order, by = c("manual_l3" = "Celltype")) %>%
  dplyr::mutate(manual_l3 = factor(manual_l3, levels = manual_l3_order$Celltype),
                expr_rank = rank(expr, ties.method="first"),
                gene = factor(gene, levels = c("HLA-DRA", "CD86", "CD14", "FCGR3A")))

# UMAP

umap_ggobj <- ggplot(expr_df, aes(x = UMAP_1, y = UMAP_2, order = expr_rank)) +
  geom_point_rast(show.legend = F, size = 0.25, col = "#d3d3d3") +
  geom_point_rast(data = expr_df %>%
                    dplyr::filter(expr>0), aes(col = expr), show.legend = T, size = 0.25) +
  labs(y = "",
       x = "") +
  guides(colour = guide_legend(override.aes = list(size = 3),
                               title = "nUMIs")) +
  facet_wrap(~gene, nrow = 2, ncol = 2) +
  scale_color_viridis() +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.pos = "bottom",
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text.x = element_text(face = "bold"))

pdf(width = 8, height = 8, file = umap_pdf)
print(umap_ggobj)
dev.off()

sessionInfo()