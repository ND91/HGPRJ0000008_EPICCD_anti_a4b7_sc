#!/usr/bin/env Rscript

# The goal of this script is to create a scatterplot of the first two UMAP dimensions of all PBMCs colored by manual_l3

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
pca_pdf <- args[3]

seuratObject <- readRDS(seurat_rds)
manual_l3_order <- readxl::read_excel(manual_l3_order_xlsx) %>%
  dplyr::mutate(number = as.numeric(factor(Celltype, levels = Celltype)),
                celltype_number = paste0(number, ". ", Celltype),
                celltype_number = factor(celltype_number, levels = celltype_number))

counts_mat <- GetAssayData(seuratObject, slot = "counts")

seuratObject@meta.data$Sample_ID_l3 <- paste0(seuratObject@meta.data$SampleID, "_", seuratObject@meta.data$manual_l3)

sample_metadata <- seuratObject@meta.data %>%
  dplyr::select(SampleID, PBMC_scrnaseq, Sex, Age, Smoker, Drug, Response) %>%
  unique()

pb_sample_l3 <- counts_mat %*% model.matrix(~0+Sample_ID_l3, data = seuratObject@meta.data)
colnames(pb_sample_l3) <- gsub("Sample_ID_l3", "", colnames(pb_sample_l3))

pb_sample_l3_metadata <- data.frame(SampleID = gsub("(F[0-9]+)_.+$", "\\1", colnames(pb_sample_l3)),
                                    Celltype = gsub("F[0-9]+_(.+)$", "\\1", colnames(pb_sample_l3))) %>%
  dplyr::left_join(sample_metadata, by = "SampleID")

pb_sample_l3 <- pb_sample_l3[which(Matrix::rowSums(pb_sample_l3) != 0),]

vst_sample_l3 <- vst(as.matrix(pb_sample_l3))

svd_sample_l3 <- svd(t(vst_sample_l3))

pca_df <- data.frame(PC1 = svd_sample_l3$u[,1],
                     PC2 = svd_sample_l3$u[,2],
                     pb_sample_l3_metadata) %>%
  dplyr::left_join(manual_l3_order, by = c("Celltype")) %>%
  dplyr::mutate(Celltype = factor(Celltype, levels = manual_l3_order$Celltype))

manual_l3_colors <- manual_l3_order$Color
names(manual_l3_colors) <- manual_l3_order$celltype_number

ggobj <- ggplot(pca_df, aes(x = PC1, y = PC2, col = celltype_number)) +
  #geom_point_rast(show.legend = T, size = 0.5, col = "black") +
  #geom_point_rast(show.legend = T, size = 0.25) +
  geom_point_rast(aes(shape = Response), size = 3) +
  labs(title = "PCA") +
  geom_label_repel(data = . %>%
                     dplyr::group_by(number, Celltype, celltype_number) %>%
                     summarize(x = median(x = PC1),
                               y = median(x = PC2)),
                   mapping = aes(label = Celltype, x = x, y = y, fill = celltype_number),
                   alpha = 0.75,
                   show.legend = F,
                   col = "black", max.overlaps = 100) +
  guides(colour = guide_legend(override.aes = list(size = 3))) +
  scale_color_manual(values = manual_l3_colors, drop = F) +
  scale_fill_manual(values = manual_l3_colors, drop = F) +
  theme_bw() +
  #facet_wrap(~Celltype, scales = "free") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.pos = "bottom",
        legend.title = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text.x = element_text(face = "bold"))

pdf(width = 10, height = 12.5, file = pca_pdf)
print(ggobj)
dev.off()

sessionInfo()