#!/usr/bin/env Rscript

# The goal of this script is to create figure 3 for the manuscript.
# Legend: Figure 3 Differential abundance identifies. Uniform manifold approximation and projection (UMAP) visualization of the PBMCs as obtained through (A) single-cell RNA-sequencing (scRNAseq) and (B) mass cytometry by time of flight (CyTOF) from CD patients on VDZ that respond (R; N = 4) and that do not respond (NR; N = 4) colored by the cell identity. Expression of the markers used to annotate the PBMCs at the level of (C) gene expression through a dotplot where size and color intensity represent the percentage cells with measurable expression and the median expression, respectively, and a (D) protein expression through a heatmap with the color representing the median expression. (E) Scatterplot representing the percentage celltypes relative to all PBMCs for scRNAseq on the x-axis and CyTOF on the y-axis colored by lineage. (F) Visualization of ITGA4 and ITGB7 gene expression as visualized through UMAP (left) and a jitterplot superposed on a boxplot with celltype on the x-axis and normalized gene expression on the y-axis. (G) Visualization of integrin α4 and integrin α4β7 protein expression as visualized through UMAP (left) and a jitterplot superposed on a boxplot with celltype on the x-axis and normalized gene expression on the y-axis.

library(dplyr)
library(ggplot2)
library(ggrastr)
library(ggrepel)
library(ggpubr)
library(ggforce)
library(Seurat)
library(SingleCellExperiment)
library(viridis)
library(scales)
library(DESeq2)
library(png)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 9) {
  stop(paste0("Script needs 9 arguments. Current input is:", args))
}

seurat_rds <- args[1]
scrnaseq_dacs_csv <- args[2]
scrnaseq_degs_l3_rds <- args[3]
gse134809_seurat_rds <- args[4]
facs_flowexample_png <- args[5]
facs_proportions_csv <- args[6]
facs_dacs_csv <- args[7]
manual_l3_order_xlsx <- args[8]
fig3_pdf <- args[9]

scrnaseq_seuratObject <- readRDS(seurat_rds)
scrnaseq_dacs <- read.csv(scrnaseq_dacs_csv)
scrnaseq_degs_l3 <- readRDS(scrnaseq_degs_l3_rds)
gse134809_seuratObject <- readRDS(gse134809_seurat_rds)
facs_proportions <- read.csv(facs_proportions_csv)
facs_dacs <- read.csv(facs_dacs_csv)
manual_l3_order <- readxl::read_excel(manual_l3_order_xlsx, col_names = T)  %>%
  dplyr::mutate(number = as.numeric(factor(Celltype, levels = Celltype)),
                celltype_number = paste0(number, ". ", Celltype),
                celltype_number = factor(celltype_number, levels = celltype_number))

manual_l3_colors <- manual_l3_order$Color
names(manual_l3_colors) <- manual_l3_order$celltype_number

response_colors <- c(`Responder` = "#009E73",
                     `Non-responder` = "#D55E00")

# Fig3A

fig3A <- scrnaseq_seuratObject@meta.data %>%
  dplyr::mutate(manual_l3 = as.factor(manual_l3)) %>%
  dplyr::group_by(SampleID, Response, manual_l3, .drop = F) %>%
  dplyr::summarize(Nl3sample = n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(SampleID, Response) %>%
  dplyr::mutate(Ncells = sum(Nl3sample),
                Ncellprop = Nl3sample/Ncells,
                Ncellprop = ifelse(is.na(Ncellprop), 0, Ncellprop)) %>%
  ggplot(aes(x = forcats::fct_reorder(manual_l3, -Ncellprop), y = Ncellprop*100, col = Response)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  geom_point(position = position_dodge(width=0.75), alpha = 0.5) +
  labs(title = "Percentage relative to PBMCs",
       y = "%PBMCs") +
  scale_color_manual(values = response_colors) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        legend.pos = "bottom",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# Fig3B

# fig3B_left <- scrnaseq_seuratObject@meta.data %>%
#   dplyr::mutate(manual_l3 = as.factor(manual_l3)) %>%
#   dplyr::group_by(SampleID, Response, manual_l3, .drop = F) %>%
#   dplyr::summarize(Nl3sample = n()) %>%
#   dplyr::ungroup() %>%
#   dplyr::group_by(SampleID, Response) %>%
#   dplyr::mutate(Ncells = sum(Nl3sample),
#                 Ncellprop = Nl3sample/Ncells,
#                 Ncellprop = ifelse(is.na(Ncellprop), 0, Ncellprop)) %>%
#   dplyr::filter(manual_l3 == "MAIT") %>%
#   ggplot(aes(x = Response, y = Ncellprop, col = Response)) +
#   geom_boxplot(alpha = 0.5, outlier.shape = NA, show.legend = F) +
#   geom_jitter(alpha = 0.5, show.legend = F) +
#   labs(title = "MAIT",
#        subtitle = paste0("p-value = ", round(scrnaseq_dacs$P.Value[scrnaseq_dacs$BaselineProp.clusters == "MAIT"], 3)),
#        y = "Proportion") +
#   scale_color_manual(values = response_colors) +
#   theme_bw() +
#   theme(panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(),
#         axis.title.x = element_blank(),
#         legend.pos = "bottom",
#         axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

fig3B_right <- scrnaseq_seuratObject@meta.data %>%
  dplyr::mutate(manual_l3 = as.factor(manual_l3)) %>%
  dplyr::group_by(SampleID, Response, manual_l3, .drop = F) %>%
  dplyr::summarize(Nl3sample = n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(SampleID, Response) %>%
  dplyr::mutate(Ncells = sum(Nl3sample),
                Ncellprop = Nl3sample/Ncells,
                Ncellprop = ifelse(is.na(Ncellprop), 0, Ncellprop)) %>%
  dplyr::filter(manual_l3 == "PDC") %>%
  ggplot(aes(x = Response, y = Ncellprop*100, col = Response)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA, show.legend = F) +
  geom_jitter(alpha = 0.5, show.legend = F) +
  labs(title = "PDCs",
       subtitle = paste0("p-value = ", round(scrnaseq_dacs$P.Value[scrnaseq_dacs$BaselineProp.clusters == "PDC"], 3)),
       y = "%PBMCs") +
  scale_color_manual(values = response_colors) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        legend.pos = "bottom",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# Fig3C

facs_example <- readPNG(facs_flowexample_png)

fig3C <- ggplot() + 
  background_image(facs_example) +
  theme(plot.margin = margin(t=0.5, l=0.5, r=0.5, b=0.5, unit = "cm"))

# Fig3D

fig3D <- facs_proportions %>%
  dplyr::left_join(facs_dacs, by = "Celltype") %>%
  dplyr::filter(Celltype %in% c("PDC")) %>%
  dplyr::mutate(label = paste0(Celltype, "\np-value = ", round(pvalue_wilcox, 2))) %>%
  ggplot(aes(x = Response, y = Percentage_relative_PBMCs, col = Response)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  geom_point(position = position_dodge(width=0.75), alpha = 0.5) +
  scale_color_manual(values = response_colors) +
  labs(y = "%PBMCs") +
  facet_wrap(~label, ncol = 1, nrow = 2, scales = "free_y") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        legend.pos = "bottom",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# Fig3E

fig3E <- scrnaseq_degs_l3$PDC$degs %>%
  data.frame() %>%
  dplyr::mutate(Significance = ifelse(padj<0.05, "Significant", "NS"),
                label = ifelse(gene %in% c("ITGA4", "ITGB7"), gene, NA)) %>%
  ggplot(aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_vline(xintercept = 0, col = "#d3d3d3") +
  geom_hline(yintercept = 0, col = "#d3d3d3") +
  geom_point_rast(aes(col = Significance)) +
  geom_point_rast(data = . %>%
                    dplyr::filter(gene %in% c("ITGA4", "ITGB7")),
                  col = "red") +
  geom_label_repel(aes(label = label)) +
  labs(title = "PDC differential expression",
       subtitle = "Responders relative to non-responders",
       y = bquote('-'~log[10]~'(p-value)'),
       x = bquote('-'~log[2]~'(fold change)')) +
  scale_color_manual(values = c(`Significant` = "#000000", `NS` = "#d3d3d3"), drop = FALSE) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.pos = "bottom")

# Fig3F 

scrnaseq_degs_l3_PDC_rld <- rlog(scrnaseq_degs_l3$PDC$dds)

fig3F <- data.frame(SampleID = colData(scrnaseq_degs_l3_PDC_rld)$SampleID,
                    Response = colData(scrnaseq_degs_l3_PDC_rld)$Response,
                    t(assay(scrnaseq_degs_l3_PDC_rld)[c("ITGA4", "ITGB7"),])) %>%
  tidyr::pivot_longer(-c(SampleID, Response), names_to = "Gene", values_to = "Expression") %>%
  dplyr::left_join(scrnaseq_degs_l3$PDC$degs %>% 
                     data.frame() %>%
                     dplyr::filter(gene %in% c("ITGA4", "ITGB7")), 
                   by = c("Gene" = "gene")) %>%
  dplyr::mutate(label = paste0(Gene, "\np-value = ", round(pvalue, 2))) %>%
  ggplot(aes(x = Response, y = Expression, col = Response)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter_rast() +
  labs(y = bquote(log[2]~'(expression)')) +
  scale_color_manual(values = response_colors) +
  facet_wrap(~label, scales = "free_y", nrow = 2) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.pos = "bottom")

# Fig3G

gse134809_pdc_expr_df <- data.frame(CB = colnames(gse134809_seuratObject),
                                    Embeddings(gse134809_seuratObject[["umap"]]),
                                    expr = GetAssayData(gse134809_seuratObject)[c("HLA-DRA"),],
                                    gene = "HLA-DRA",
                                    gse134809_seuratObject@meta.data) %>%
  dplyr::bind_rows(data.frame(CB = colnames(gse134809_seuratObject),
                              Embeddings(gse134809_seuratObject[["umap"]]),
                              expr = GetAssayData(gse134809_seuratObject)["ITGAM",],
                              gene = "ITGAM",
                              gse134809_seuratObject@meta.data)) %>%
  dplyr::bind_rows(data.frame(CB = colnames(gse134809_seuratObject),
                              Embeddings(gse134809_seuratObject[["umap"]]),
                              expr = GetAssayData(gse134809_seuratObject)["IL3RA",],
                              gene = "IL3RA",
                              gse134809_seuratObject@meta.data)) %>%
  dplyr::bind_rows(data.frame(CB = colnames(gse134809_seuratObject),
                              Embeddings(gse134809_seuratObject[["umap"]]),
                              expr = GetAssayData(gse134809_seuratObject)["CLEC4C",],
                              gene = "CLEC4C",
                              gse134809_seuratObject@meta.data)) %>%
  dplyr::bind_rows(data.frame(CB = colnames(gse134809_seuratObject),
                              Embeddings(gse134809_seuratObject[["umap"]]),
                              expr = GetAssayData(gse134809_seuratObject)["PTPRC",],
                              gene = "PTPRC",
                              gse134809_seuratObject@meta.data)) %>%
  dplyr::bind_rows(data.frame(CB = colnames(gse134809_seuratObject),
                              Embeddings(gse134809_seuratObject[["umap"]]),
                              expr = GetAssayData(gse134809_seuratObject)["VIL1",],
                              gene = "EPCAM",
                              gse134809_seuratObject@meta.data)) %>%
  #dplyr::left_join(manual_l3_order, by = c("manual_l3" = "Celltype")) %>%
  dplyr::mutate(#manual_l3 = factor(manual_l3, levels = manual_l3_order$Celltype),
                expr_rank = rank(expr, ties.method="first"),
                gene = factor(gene, levels = c("PTPRC", "EPCAM", "HLA-DRA", "ITGAM", "IL3RA", "CLEC4C")))

fig3G_left <- ggplot(gse134809_pdc_expr_df, aes(x = UMAP_1, y = UMAP_2, order = expr_rank)) +
  geom_point_rast(show.legend = F, size = 0.5, col = "#d3d3d3") +
  geom_point_rast(data = gse134809_pdc_expr_df %>%
                    dplyr::filter(expr>0), aes(col = expr), show.legend = T, size = 0.5) +
  labs(y = "",
       x = "",
       title = "Ileal biopsy PDCs",
       subtitle = "Source: Martin et al. 2019 (GSE134809)") +
  guides(colour = guide_legend(override.aes = list(size = 3),
                               title = "Expression")) +
  facet_wrap(~gene, nrow = 2, ncol = 3) +
  scale_color_viridis() +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.pos = "left",
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text.x = element_text(face = "bold"))

fig3G_right <- data.frame(CB = colnames(gse134809_seuratObject),
                          Embeddings(gse134809_seuratObject[["umap"]]),
                          gse134809_seuratObject@meta.data) %>%
  dplyr::mutate(label = ifelse(manual_l3 == "PDC", "PDC", NA)) %>%
  ggplot(aes(x = UMAP_1, y = UMAP_2)) +
  geom_point_rast(show.legend = T, size = 0.5, col = "#d3d3d3") +
  geom_point_rast(data = . %>%
                    dplyr::filter(manual_l3 == "PDC"),
                  show.legend = T, size = 0.25, col = "#B15928") +
  # geom_mark_hull(data =
  labs(y = "",
       x = "",
       title = "",
       subtitle = "") +
  guides(colour = guide_legend(override.aes = list(size = 3))) +
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

# Fig3H

gse134809_pdcrimmune_proportions <- gse134809_seuratObject@meta.data %>%
  dplyr::filter(manual_l0 %in% c("Immune")) %>%
  dplyr::mutate(manual_l3 = as.factor(manual_l3)) %>%
  dplyr::group_by(sampleID, Phenotype, manual_l3, .drop = F) %>%
  dplyr::summarize(Nl3sample = n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(sampleID, Phenotype) %>%
  dplyr::mutate(Ncells = sum(Nl3sample),
                Ncellprop = Nl3sample/Ncells,
                Ncellprop = ifelse(is.na(Ncellprop), 0, Ncellprop)) %>%
  dplyr::filter(manual_l3 == "PDC") 

fig3H <- gse134809_pdcrimmune_proportions %>%
  ggplot(aes(x = Phenotype, y = Ncellprop*100, col = Phenotype)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA, show.legend = F) +
  geom_jitter(alpha = 0.5, show.legend = F) +
  labs(title = "Ileal PDC from GSE134809",
       subtitle = paste0("p-value = ", round(wilcox.test(gse134809_pdcrimmune_proportions$Ncellprop ~ gse134809_pdcrimmune_proportions$Phenotype, paired = T)$p.value, 3)),
       y = "%Immune cells") +
  scale_color_manual(values = c(`Involved` = "#CC79A7", `Uninvolved` = "#E69F00")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        legend.pos = "bottom",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# Compiled

fig3AB <- ggarrange(fig3A, fig3B_right, nrow = 1, ncol = 2, labels = c("A", ""), common.legend = T, legend = "bottom", align = "hv", widths = c(4, 0.5))
#width = 18, height = 8.1
#fig3CD <- ggarrange(fig3C, fig3D, nrow = 1, ncol = 2, labels = c("C", "D"), align = "hv", widths = c(1, 0.25))
#fig3DEF <- ggarrange(fig3D, fig3E, fig3F, nrow = 1, ncol = 3, labels = c("D", "E", "F"), align = "hv", widths = c(0.5, 1, 0.5))
fig3CDEF <- ggarrange(ggplot(), fig3D, fig3E, fig3F, nrow = 1, ncol = 4, labels = c("B", "C", "D", "E"), align = "hv", widths = c(1, 0.25, 0.5, 0.25))
#fig3EFGH <- ggarrange(fig3E, fig3F, fig3G_left, fig3G_right, fig3H, nrow = 1, ncol = 5, labels = c("E", "F", "G", "", "H"), align = "hv", widths = c(1, 0.5, 1.5, 1, 0.5))
fig3GH <- ggarrange(fig3G_left, fig3G_right, fig3H, nrow = 1, ncol = 3, labels = c("F", "", "G"), align = "hv", widths = c(1.5, 1, 0.5))
#width = 18, height = 5

#fig3 <- ggarrange(plotlist = list(fig3AB, fig3C, fig3DEF, fig3GH), nrow = 4, ncol = 1, heights = c(1, 0.75, 0.75, 0.75))
#fig3 <- ggarrange(plotlist = list(fig3AB, fig3CD, fig3EFGH), nrow = 3, ncol = 1, heights = c(1, 1, 1))
fig3 <- ggarrange(plotlist = list(fig3AB, fig3CDEF, fig3GH), nrow = 3, ncol = 1, heights = c(1, 1, 1))

pdf(width = 16.5, height = 15, file = fig3_pdf)
print(fig3)
dev.off()

sessionInfo()