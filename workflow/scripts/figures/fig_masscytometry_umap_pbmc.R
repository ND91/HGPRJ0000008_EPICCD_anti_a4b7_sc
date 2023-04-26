#!/usr/bin/env Rscript

# The goal of this script is to create a scatterplot of the first two UMAP dimensions of all PBMCs colored by manual_l3 as obtained through mass cytometry

library(dplyr)
library(ggplot2)
library(ggrastr)
library(ggrepel)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop(paste0("Script needs 3 arguments. Current input is:", args))
}

umap_csv <- args[1]
manual_l3_order_xlsx <- args[2]
umap_pdf <- args[3]

umap_df <- read.csv(umap_csv)
manual_l3_order <- readxl::read_excel(manual_l3_order_xlsx) %>%
  dplyr::mutate(number = as.numeric(factor(Celltype, levels = Celltype)),
                celltype_number = paste0(number, ". ", Celltype),
                celltype_number = factor(celltype_number, levels = celltype_number))

umap_df <- umap_df %>%
  dplyr::left_join(manual_l3_order, by = c("manual_l3" = "Celltype")) %>%
  dplyr::mutate(manual_l3 = factor(manual_l3, levels = manual_l3_order$Celltype))

manual_l3_colors <- manual_l3_order$Color
names(manual_l3_colors) <- manual_l3_order$celltype_number

ggobj <- ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, col = celltype_number)) +
  geom_point_rast(show.legend = T, size = 0.5, col = "black") +
  geom_point_rast(show.legend = T, size = 0.25) +
  labs(y = "",
       x = "",
       title = "Mass cytometry",
       subtitle = paste0("subsampled to ", nrow(umap_df), " cells")) +
  geom_label_repel(data = umap_df %>%
                     dplyr::group_by(number, celltype_number) %>%
                     summarize(x = median(x = UMAP_1),
                               y = median(x = UMAP_2)),
                   mapping = aes(label = number, x = x, y = y, fill = celltype_number),
                   alpha = 0.75, 
                   show.legend = F,
                   col = "black") +
  guides(colour = guide_legend(override.aes = list(size = 3))) +
  scale_color_manual(values = manual_l3_colors, drop = F) +
  scale_fill_manual(values = manual_l3_colors, drop = F) +
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

pdf(width = 10, height = 10, file = umap_pdf)
print(ggobj)
dev.off()

sessionInfo()