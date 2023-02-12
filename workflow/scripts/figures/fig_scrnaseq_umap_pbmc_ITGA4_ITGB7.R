#!/usr/bin/env Rscript

# The goal of this script is to create a scatterplot of the first two UMAP dimensions of all PBMCs colored by manual_l3

library(Seurat)
library(dplyr)
library(ggplot2)
library(ggrastr)
library(ggrepel)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop(paste0("Script needs 2 arguments. Current input is:", args))
}

seurat_rds <- args[1]
manual_l3_order_xlsx <- args[2]
umap_pdf <- args[3]

seuratObject <- readRDS(seurat_rds)
manual_l3_order <- readxl::read_excel(manual_l3_order_xlsx, col_names = F)

umap_df <- data.frame(CB = colnames(seuratObject),
                      Embeddings(seuratObject[["umap"]]),
                      expr = GetAssayData(seuratObject)["ITGA4",],
                      gene = "ITGA4",
                      seuratObject@meta.data) %>%
  dplyr::bind_rows(data.frame(CB = colnames(seuratObject),
                              Embeddings(seuratObject[["umap"]]),
                              expr = GetAssayData(seuratObject)["ITGB7",],
                              gene = "ITGB7",
                              seuratObject@meta.data)) %>%
  dplyr::mutate(manual_l3 = factor(manual_l3, levels = pull(manual_l3_order)),
                manual_l3_number = as.numeric(manual_l3),
                manual_l3_w_number = paste0(as.numeric(manual_l3), ". ", manual_l3))

manual_l3_w_number_order <- umap_df %>% 
  dplyr::select(manual_l3_number, manual_l3_w_number) %>% 
  unique() %>% arrange(manual_l3_number) %>% 
  pull(manual_l3_w_number)

umap_df$manual_l3_w_number <- factor(umap_df$manual_l3_w_number, levels = manual_l3_w_number_order)

plotobj <- ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point_rast(show.legend = F, size = 0.5, col = "#d3d3d3") +
  geom_point_rast(data = umap_df %>%
                    dplyr::filter(expr>0), aes(col = expr), show.legend = T, size = 0.5) +
  labs(y = "",
       x = "") +
  # geom_label_repel(data = umap_df %>%
  #                    dplyr::group_by(manual_l3_number, manual_l3_w_number) %>%
  #                    summarize(x = median(x = UMAP_1),
  #                              y = median(x = UMAP_2)),
  #                  mapping = aes(label = manual_l3_w_number, x = x, y = y),
  #                  alpha = 0.5, 
  #                  show.legend = F) +
  guides(colour = guide_legend(override.aes = list(size = 3))) +
  facet_wrap(~gene, nrow = 1, ncol = 2) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.pos = "bottom",
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text.x = element_text(face = "bold"))

pdf(width = 20, height = 10, file = umap_pdf)
print(plotobj)
dev.off()

png(width = 10, height = 5, units = "in", res = 120, file = umap_png)
print(plotobj)
dev.off()

sessionInfo()