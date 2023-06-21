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

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 8) {
  stop(paste0("Script needs 8 arguments. Current input is:", args))
}

seurat_rds <- args[1]
scrnaseq_dacs_csv <- args[2]
scrnaseq_degs_l3_rds <- args[3]
gse134809_seurat_rds <- args[4]
facs_proportions_csv <- args[5]
facs_dacs_csv <- args[6]
manual_l3_order_xlsx <- args[7]
fig4_pdf <- args[8]

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

manual_l3_number_colors <- manual_l3_order$Color
names(manual_l3_number_colors) <- manual_l3_order$celltype_number

response_colors <- c(`Responder` = "#009E73",
                     `Non-responder` = "#D55E00")

# fig4B

PDC_perc_PBMC <- scrnaseq_seuratObject@meta.data %>%
  dplyr::mutate(manual_l3 = as.factor(manual_l3)) %>%
  dplyr::group_by(SampleID, Response, manual_l3, .drop = F) %>%
  dplyr::summarize(Nl3sample = n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(SampleID, Response) %>%
  dplyr::mutate(Ncells = sum(Nl3sample),
                Ncellprop = Nl3sample/Ncells,
                Ncellprop = ifelse(is.na(Ncellprop), 0, Ncellprop),
                Ncellperc = Ncellprop*100,
                Modality = "scRNAseq") %>%
  dplyr::filter(manual_l3 == "PDC") %>%
  dplyr::left_join(scrnaseq_dacs, by = c("manual_l3" = "BaselineProp.clusters")) %>%
  dplyr::select(-c(X, Nl3sample, Ncellprop, Ncells, BaselineProp.Freq, PropMean.Non.responder, PropMean.Responder, PropRatio, Tstatistic, FDR)) %>%
  dplyr::rename(Celltype = manual_l3) %>%
  dplyr::ungroup() %>%
  dplyr::add_row(
    facs_proportions %>%
      dplyr::left_join(facs_dacs, by = "Celltype") %>%
      dplyr::filter(Celltype %in% c("PDC")) %>%
      dplyr::select(-c(X.x, Run_ID, X.y)) %>%
      dplyr::rename(SampleID = Sample_ID,
                    Ncellperc = Percentage_relative_PBMCs,
                    P.Value = pvalue_wilcox) %>%
      dplyr::mutate(Modality = "FACS") %>%
      dplyr::select(SampleID, Response, Celltype, Ncellperc, Modality, P.Value)
  ) %>%
  dplyr::mutate(label = factor(paste0(Modality, "\np-value = ", round(P.Value, 3)), levels = unique(paste0(Modality, "\np-value = ", round(P.Value, 3)))))
  
fig4B <- PDC_perc_PBMC %>%
  ggplot(aes(x = Response, y = Ncellperc, col = Response)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA, show.legend = F) +
  geom_jitter(alpha = 0.5, show.legend = F) +
  labs(title = "PDC",
       y = "%PBMCs") +
  scale_color_manual(values = response_colors) +
  facet_wrap(~label, nrow = 2, scales = "free_y") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        legend.pos = "bottom",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

# fig4C

fig4C <- scrnaseq_degs_l3$PDC$degs %>%
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

# fig4ED

scrnaseq_degs_l3_PDC_rld <- rlog(scrnaseq_degs_l3$PDC$dds)

fig4D <- data.frame(SampleID = colData(scrnaseq_degs_l3_PDC_rld)$SampleID,
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
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.pos = "none")

# fig4E

gse134809_pdc_expr_df <- data.frame(CB = colnames(gse134809_seuratObject),
                                    Embeddings(gse134809_seuratObject[["umap"]]),
                                    expr = GetAssayData(gse134809_seuratObject)[c("HLA-DRA"),],
                                    gene = "HLA-DRA",
                                    gse134809_seuratObject@meta.data) %>%
  dplyr::bind_rows(data.frame(CB = colnames(gse134809_seuratObject),
                              Embeddings(gse134809_seuratObject[["umap"]]),
                              expr = GetAssayData(gse134809_seuratObject)["ITGAX",],
                              gene = "ITGAX",
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
                gene = factor(gene, levels = c("PTPRC", "EPCAM", "HLA-DRA", "ITGAX", "IL3RA", "CLEC4C")))

fig4E_left <- ggplot(gse134809_pdc_expr_df, aes(x = UMAP_1, y = UMAP_2, order = expr_rank)) +
  geom_point_rast(show.legend = F, size = 0.25, col = "#d3d3d3") +
  geom_point_rast(data = gse134809_pdc_expr_df %>%
                    dplyr::filter(expr>1), aes(col = expr), show.legend = T, size = 0.25) +
  labs(y = "",
       x = "",
       title = "Ileal biopsies",
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

fig4E_right <- data.frame(CB = colnames(gse134809_seuratObject),
                          Embeddings(gse134809_seuratObject[["umap"]]),
                          gse134809_seuratObject@meta.data) %>%
  dplyr::mutate(label = ifelse(manual_l3 == "PDC", "PDC", NA)) %>%
  ggplot(aes(x = UMAP_1, y = UMAP_2)) +
  geom_point_rast(show.legend = T, size = 0.5, col = "#d3d3d3") +
  geom_point_rast(data = . %>%
                    dplyr::filter(manual_l3 == "PDC"),
                  show.legend = T, size = 0.25, col = "#B15928") +
  labs(y = "",
       x = "",
       title = "",
       subtitle = "PDCs") +
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

# fig4F

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

fig4F <- gse134809_pdcrimmune_proportions %>%
  ggplot(aes(x = Phenotype, y = Ncellprop*100, col = Phenotype)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA, show.legend = F) +
  geom_jitter(alpha = 0.5, show.legend = F, height = 0) +
  labs(title = "Ileal PDC abundance",
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

fig4ABCDE <- ggarrange(ggplot(), fig4B, fig4C, fig4D, nrow = 1, ncol = 4, labels = c("A", "B", "C", "D"), align = "hv", widths = c(1, 0.25, 0.5, 0.25))
fig4FG <- ggarrange(fig4E_left, fig4E_right, fig4F, nrow = 1, ncol = 3, labels = c("E", "", "F"), align = "hv", widths = c(1.5, 1, 0.5))
#width = 18, height = 5

fig4 <- ggarrange(plotlist = list(fig4ABCDE, fig4FG), nrow = 2, ncol = 1, heights = c(1, 1))

pdf(width = 17, height = 11.5, file = fig4_pdf)
print(fig4)
dev.off()

sessionInfo()