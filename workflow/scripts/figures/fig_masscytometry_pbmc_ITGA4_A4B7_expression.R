#!/usr/bin/env Rscript

# The goal of this script is to create a scatterplot of the first two UMAP dimensions of all PBMCs colored by ITGA4 and A4B7 protein expression (VDZ-binding) as obtained through mass cytometry.

library(SingleCellExperiment)
library(dplyr)
library(ggplot2)
library(ggrastr)
library(ggrepel)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) {
  stop(paste0("Script needs 5 arguments. Current input is:", args))
}

sce_ss_rds <- args[1]
sce_full_rds <- args[2]
manual_l3_order_xlsx <- args[3]
umap_pdf <- args[4]
boxplot_pdf <- args[5]

sce_ss <- readRDS(sce_ss_rds)
sce_full <- readRDS(sce_full_rds)

manual_l3_order <- readxl::read_excel(manual_l3_order_xlsx) %>%
  dplyr::mutate(number = as.numeric(factor(Celltype, levels = Celltype)),
                celltype_number = paste0(number, ". ", Celltype),
                celltype_number = factor(celltype_number, levels = celltype_number))

# UMAP

umap_df <- as.matrix(assay(sce_ss)[c("ITGA4", "A4B7"),]) %>%
  data.frame(., protein = rownames(.)) %>%
  tidyr::pivot_longer(-protein, names_to = "cellID", values_to = "expr") %>%
  dplyr::left_join(data.frame(colData(sce_ss), reducedDims(sce_ss)[[1]], cellID = colnames(sce_ss)),
                   by = "cellID") %>%
  dplyr::left_join(manual_l3_order, by = c("manual_l3" = "Celltype")) %>%
  dplyr::mutate(manual_l3 = factor(manual_l3, levels = manual_l3_order$Celltype),
                expr_rank = rank(expr, ties.method="first"),
                label = factor(paste0("Mass cytometry: ", protein), levels = c("Mass cytometry: ITGA4", "Mass cytometry: A4B7")),
                protein = factor(protein, levels = c("ITGA4", "A4B7")))

umap_ggobj <- ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, order = expr_rank)) +
  geom_point_rast(show.legend = F, size = 0.5, col = "#d3d3d3") +
  geom_point_rast(data = umap_df %>%
                    dplyr::filter(expr>0), aes(col = expr), show.legend = T, size = 0.5) +
  labs(y = "",
       x = "") +
  guides(colour = guide_legend(override.aes = list(size = 3),
                               title = "Expression")) +
  facet_wrap(~label, nrow = 2, ncol = 1) +
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

pdf(width = 20, height = 10, file = umap_pdf)
print(umap_ggobj)
dev.off()

# Boxplot

boxplot_df <- as.matrix(assay(sce_full)[c("ITGA4", "A4B7"),]) %>%
  data.frame(., protein = rownames(.)) %>%
  tidyr::pivot_longer(-protein, names_to = "cellID", values_to = "expr") %>%
  dplyr::left_join(data.frame(colData(sce_full), cellID = colnames(sce_full)),
                   by = "cellID") %>%
  dplyr::left_join(manual_l3_order, by = c("manual_l3" = "Celltype")) %>%
  dplyr::mutate(manual_l3 = factor(manual_l3, levels = manual_l3_order$Celltype),
                expr_rank = rank(expr, ties.method="first"),
                label = factor(paste0("Mass cytometry: ", protein), levels = c("Mass cytometry: ITGA4", "Mass cytometry: A4B7")),
                protein = factor(protein, levels = c("ITGA4", "A4B7")))

boxplot_ggobj <- ggplot(boxplot_df, aes(x = manual_l3, y = expr)) +
  geom_jitter_rast(alpha = 0.1) +
  geom_boxplot(alpha = 0.75, outlier.shape = NA) +
  facet_wrap(~label, nrow = 2, ncol = 1) +
  labs(y = "Expression") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.pos = "bottom",
        axis.title.x = element_blank(),
        strip.text.x = element_text(face = "bold"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

pdf(width = 20, height = 10, file = boxplot_pdf)
print(boxplot_ggobj)
dev.off()

sessionInfo()