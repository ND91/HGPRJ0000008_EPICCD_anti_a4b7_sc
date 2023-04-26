#!/usr/bin/env Rscript

# The goal of this script is to create a boxplot of the response criteria between responders and non-responders.

library(dplyr)
library(ggplot2)
library(ggrastr)
library(ggrepel)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop(paste0("Script needs 3 arguments. Current input is:", args))
}

sample_metadata_xlsx <- args[1]
response_cols <- c(`Non-responder` = "#0C7BDC", `Responder` =  "#FFC20A")
boxplot_pdf <- args[3]

sample_metadata <- readxl::read_excel(sample_metadata_xlsx) %>%
  dplyr::mutate(perc_SES_T1 = (as.numeric(SES_T1)-as.numeric(SES_T2))/as.numeric(SES_T1)*100,
                perc_fCalpro_T1 = (as.numeric(fCalpro_T1)-as.numeric(fCalpro_T2))/as.numeric(fCalpro_T1)*100,
                perc_CRP_T1 = (as.numeric(CRP_T1)-as.numeric(CRP_T2))/as.numeric(CRP_T1)*100,
                drop_HBI = as.numeric(HBI_T1)-as.numeric(HBI_T2)) %>%
  dplyr::filter(PBMC_scrnaseq == "Yes")

SES_ggobj <- sample_metadata %>%
  dplyr::select(Donor_ID, Response, perc_SES_T1) %>%
  ggplot(aes(x = Response, y = perc_SES_T1, fill = Response)) +
  geom_hline(yintercept = 0) +
  geom_jitter(shape = 21) +
  geom_boxplot(outlier.shape = NA, aes(fill = Response), alpha = 0.5) + 
  scale_fill_manual(values = response_cols) +
  labs(title = "SES-CD",
       subtitle = "Responders < 50%",
       y = "% drop from T1") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        title = element_text(face = "bold"), 
        legend.pos = "bottom")

pdf(width = 10, height = 10, file = umap_pdf)
print(SES_ggobj)
dev.off()

fCalpro_ggobj <- sample_metadata %>%
  dplyr::select(Donor_ID, Response, perc_fCalpro_T1) %>%
  ggplot(aes(x = Response, y = perc_fCalpro_T1, fill = Response)) +
  geom_hline(yintercept = 0) +
  geom_jitter(shape = 21) +
  geom_boxplot(outlier.shape = NA, aes(fill = Response), alpha = 0.5) + 
  scale_fill_manual(values = response_cols) +
  labs(title = "Fecal calprotectin",
       subtitle = "Responders < 50%",
       y = "% drop from T1") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        title = element_text(face = "bold"), 
        legend.pos = "bottom")

pdf(width = 10, height = 10, file = umap_pdf)
print(fCalpro_ggobj)
dev.off()

CRP_ggobj <- sample_metadata %>%
  dplyr::select(Donor_ID, Response, perc_CRP_T1) %>%
  ggplot(aes(x = Response, y = perc_CRP_T1, fill = Response)) +
  geom_hline(yintercept = 0) +
  geom_jitter(shape = 21) +
  geom_boxplot(outlier.shape = NA, aes(fill = Response), alpha = 0.5) + 
  scale_fill_manual(values = response_cols) +
  labs(title = "CRP",
       subtitle = "Responders < 50%",
       y = "% drop from T1") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        title = element_text(face = "bold"), 
        legend.pos = "bottom")

pdf(width = 10, height = 10, file = umap_pdf)
print(CRP_ggobj)
dev.off()

HBI_ggobj <- sample_metadata %>%
  dplyr::select(Donor_ID, Response, drop_HBI) %>%
  ggplot(aes(x = Response, y = drop_HBI, fill = Response)) +
  geom_hline(yintercept = 0) +
  geom_jitter(shape = 21) +
  geom_boxplot(outlier.shape = NA, aes(fill = Response), alpha = 0.5) + 
  scale_fill_manual(values = response_cols) +
  labs(title = "HBI",
       subtitle = "Responders < 3",
       y = "Points drop from T1") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        title = element_text(face = "bold"), 
        legend.pos = "bottom")

sessionInfo()