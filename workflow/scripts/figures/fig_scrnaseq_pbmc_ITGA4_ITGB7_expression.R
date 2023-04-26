#!/usr/bin/env Rscript

# The goal of this script is to create a scatterplot of the first two UMAP dimensions of all PBMCs colored by ITGA4 and ITGB7 gene expression as obtained through scRNAseq.

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
                      expr = GetAssayData(seuratObject)[c("ITGA4"),],
                      gene = "ITGA4",
                      seuratObject@meta.data) %>%
  dplyr::bind_rows(data.frame(CB = colnames(seuratObject),
                              Embeddings(seuratObject[["umap"]]),
                              expr = GetAssayData(seuratObject)["ITGB7",],
                              gene = "ITGB7",
                              seuratObject@meta.data)) %>%
  dplyr::left_join(manual_l3_order, by = c("manual_l3" = "Celltype")) %>%
  dplyr::mutate(manual_l3 = factor(manual_l3, levels = manual_l3_order$Celltype),
                expr_rank = rank(expr, ties.method="first"),
                label = paste0("scRNAseq: ", gene))

# UMAP

umap_ggobj <- ggplot(expr_df, aes(x = UMAP_1, y = UMAP_2, order = expr_rank)) +
  geom_point_rast(show.legend = F, size = 0.5, col = "#d3d3d3") +
  geom_point_rast(data = expr_df %>%
                    dplyr::filter(expr>0), aes(col = expr), show.legend = T, size = 0.5) +
  labs(y = "",
       x = "") +
  guides(colour = guide_legend(override.aes = list(size = 3),
                               title = "nUMIs")) +
  facet_wrap(~label, nrow = 2, ncol = 1) +
  scale_color_viridis() +
  theme_minimal() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.pos = "bottom",
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text.x = element_text(face = "bold"))

pdf(width = 16, height = 8, file = umap_pdf)
print(umap_ggobj)
dev.off()

# Boxplot

boxplot_ggobj <- ggplot(expr_df, aes(x = manual_l3, y = expr)) +
  geom_jitter_rast(alpha = 0.5) +
  geom_boxplot(alpha = 0.75, outlier.shape = NA) +
  facet_wrap(~label, nrow = 2, ncol = 1) +
  labs(y = "nUMIs") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.pos = "bottom",
        axis.title.x = element_blank(),
        strip.text.x = element_text(face = "bold"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

pdf(width = 16, height = 8, file = boxplot_pdf)
print(boxplot_ggobj)
dev.off()

sessionInfo()