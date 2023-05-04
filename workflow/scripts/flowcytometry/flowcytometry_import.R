#!/usr/bin/env Rscript

# The goal of this script is to import and preprocess data obtained from FlowJo.

library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop(paste0("Script needs 4 arguments. Current input is:", args))
}

flowjo_wsp_xlsx <- args[1]
sample_metadata_xlsx <- args[2]
files_metadata_xlsx <- args[3]
celltype_percentage_csv <- args[4]

sample_metadata <- readxl::read_excel(sample_metadata_xlsx)
files_metadata <- readxl::read_excel(files_metadata_xlsx)

sample_metadata_complete <- files_metadata %>%
  dplyr::left_join(sample_metadata, by = "Sample_ID")

flowjo_wsp <- readxl::read_excel(flowjo_wsp_xlsx) %>%
  dplyr::rename(Run_ID = ...1) %>%
  dplyr::filter(!Run_ID %in% c("Mean", "SD")) %>%
  dplyr::mutate(Run_ID = gsub(".fcs", "", Run_ID),
                `PBMCs/Single Cells/Single Cells/Lin-ve LD-ve | Freq. of Parent` = as.numeric(gsub(" %", "", `PBMCs/Single Cells/Single Cells/Lin-ve LD-ve | Freq. of Parent`)),
                `PBMCs/Single Cells/Single Cells/Lin-ve LD-ve/MHCII+ve | Freq. of Parent` = as.numeric(gsub(" %", "", `PBMCs/Single Cells/Single Cells/Lin-ve LD-ve/MHCII+ve | Freq. of Parent`)),
                `PBMCs/Single Cells/Single Cells/Lin-ve LD-ve/MHCII+ve/CD14++CD16- | Freq. of Parent` = as.numeric(gsub(" %", "", `PBMCs/Single Cells/Single Cells/Lin-ve LD-ve/MHCII+ve/CD14++CD16- | Freq. of Parent`)),
                `PBMCs/Single Cells/Single Cells/Lin-ve LD-ve/MHCII+ve/CD14+CD16+ | Freq. of Parent` = as.numeric(gsub(" %", "", `PBMCs/Single Cells/Single Cells/Lin-ve LD-ve/MHCII+ve/CD14+CD16+ | Freq. of Parent`)),
                `PBMCs/Single Cells/Single Cells/Lin-ve LD-ve/MHCII+ve/CD14+CD16++ | Freq. of Parent` = as.numeric(gsub(" %", "", `PBMCs/Single Cells/Single Cells/Lin-ve LD-ve/MHCII+ve/CD14+CD16++ | Freq. of Parent`)),
                `PBMCs/Single Cells/Single Cells/Lin-ve LD-ve/MHCII+ve/CD14-ve CD16-ve | Freq. of Parent` = as.numeric(gsub(" %", "", `PBMCs/Single Cells/Single Cells/Lin-ve LD-ve/MHCII+ve/CD14-ve CD16-ve | Freq. of Parent`)),
                `PBMCs/Single Cells/Single Cells/Lin-ve LD-ve/MHCII+ve/CD14-ve CD16-ve/CD11c- CD123+ | Freq. of Parent` = as.numeric(gsub(" %", "", `PBMCs/Single Cells/Single Cells/Lin-ve LD-ve/MHCII+ve/CD14-ve CD16-ve/CD11c- CD123+ | Freq. of Parent`)),
                `PBMCs/Single Cells/Single Cells/Lin-ve LD-ve/MHCII+ve/CD14-ve CD16-ve/CD123- | Freq. of Parent` = as.numeric(gsub(" %", "", `PBMCs/Single Cells/Single Cells/Lin-ve LD-ve/MHCII+ve/CD14-ve CD16-ve/CD123- | Freq. of Parent`)),
                `PBMCs/Single Cells/Single Cells/Lin-ve LD-ve/MHCII+ve/CD14-ve CD16-ve/CD123-/CD11c+ CD1c+ | Freq. of Parent` = as.numeric(gsub(" %", "", `PBMCs/Single Cells/Single Cells/Lin-ve LD-ve/MHCII+ve/CD14-ve CD16-ve/CD123-/CD11c+ CD1c+ | Freq. of Parent`)),
                `PBMCs/Single Cells/Single Cells/Lin-ve LD-ve/MHCII+ve/CD14-ve CD16-ve/CD11c- CD123+ | Freq. of root` = `PBMCs/Single Cells/Single Cells/Lin-ve LD-ve/MHCII+ve/CD14-ve CD16-ve/CD11c- CD123+ | Freq. of Parent`/100*
                  `PBMCs/Single Cells/Single Cells/Lin-ve LD-ve/MHCII+ve/CD14-ve CD16-ve | Freq. of Parent`/100*
                  `PBMCs/Single Cells/Single Cells/Lin-ve LD-ve/MHCII+ve | Freq. of Parent`/100*
                  `PBMCs/Single Cells/Single Cells/Lin-ve LD-ve | Freq. of Parent`,
                `PBMCs/Single Cells/Single Cells/Lin-ve LD-ve/MHCII+ve/CD14++CD16- | Freq. of root` = `PBMCs/Single Cells/Single Cells/Lin-ve LD-ve/MHCII+ve/CD14++CD16- | Freq. of Parent`/100*
                  `PBMCs/Single Cells/Single Cells/Lin-ve LD-ve/MHCII+ve | Freq. of Parent`/100*
                  `PBMCs/Single Cells/Single Cells/Lin-ve LD-ve | Freq. of Parent`,
                `PBMCs/Single Cells/Single Cells/Lin-ve LD-ve/MHCII+ve/CD14+CD16+ | Freq. of root` = `PBMCs/Single Cells/Single Cells/Lin-ve LD-ve/MHCII+ve/CD14+CD16+ | Freq. of Parent`/100*
                  `PBMCs/Single Cells/Single Cells/Lin-ve LD-ve/MHCII+ve | Freq. of Parent`/100*
                  `PBMCs/Single Cells/Single Cells/Lin-ve LD-ve | Freq. of Parent`,
                `PBMCs/Single Cells/Single Cells/Lin-ve LD-ve/MHCII+ve/CD14+CD16++ | Freq. of root` = `PBMCs/Single Cells/Single Cells/Lin-ve LD-ve/MHCII+ve/CD14+CD16++ | Freq. of Parent`/100*
                  `PBMCs/Single Cells/Single Cells/Lin-ve LD-ve/MHCII+ve | Freq. of Parent`/100*
                  `PBMCs/Single Cells/Single Cells/Lin-ve LD-ve | Freq. of Parent`,
                `PBMCs/Single Cells/Single Cells/Lin-ve LD-ve/MHCII+ve/CD14-ve CD16-ve/CD123-/CD11c+ CD1c+ | Freq. of root` = `PBMCs/Single Cells/Single Cells/Lin-ve LD-ve/MHCII+ve/CD14-ve CD16-ve/CD123-/CD11c+ CD1c+ | Freq. of Parent`/100*
                  `PBMCs/Single Cells/Single Cells/Lin-ve LD-ve/MHCII+ve/CD14-ve CD16-ve/CD123- | Freq. of Parent`/100*
                  `PBMCs/Single Cells/Single Cells/Lin-ve LD-ve/MHCII+ve/CD14-ve CD16-ve | Freq. of Parent`/100*
                  `PBMCs/Single Cells/Single Cells/Lin-ve LD-ve/MHCII+ve | Freq. of Parent`/100*
                  `PBMCs/Single Cells/Single Cells/Lin-ve LD-ve | Freq. of Parent`) %>%
  dplyr::left_join(sample_metadata_complete, by = "Run_ID")

celltype_percentage <- flowjo_wsp %>%
  dplyr::select(Run_ID,
                Sample_ID,
                Response,
                `PBMCs/Single Cells/Single Cells/Lin-ve LD-ve/MHCII+ve/CD14++CD16- | Freq. of root`, #Classical monocytes
                `PBMCs/Single Cells/Single Cells/Lin-ve LD-ve/MHCII+ve/CD14+CD16+ | Freq. of root`, #Intermediate monocytes
                `PBMCs/Single Cells/Single Cells/Lin-ve LD-ve/MHCII+ve/CD14+CD16++ | Freq. of root`, #Non-classical monocytes
                `PBMCs/Single Cells/Single Cells/Lin-ve LD-ve/MHCII+ve/CD14-ve CD16-ve/CD123-/CD11c+ CD1c+ | Freq. of root`, #CDCs
                `PBMCs/Single Cells/Single Cells/Lin-ve LD-ve/MHCII+ve/CD14-ve CD16-ve/CD11c- CD123+ | Freq. of root`) %>% #PDC
  dplyr::rename(`Classical monocyte` = `PBMCs/Single Cells/Single Cells/Lin-ve LD-ve/MHCII+ve/CD14++CD16- | Freq. of root`, 
                `Intermediate monocyte` = `PBMCs/Single Cells/Single Cells/Lin-ve LD-ve/MHCII+ve/CD14+CD16+ | Freq. of root`, 
                `Non-classical monocyte` = `PBMCs/Single Cells/Single Cells/Lin-ve LD-ve/MHCII+ve/CD14+CD16++ | Freq. of root`, 
                CDC2 = `PBMCs/Single Cells/Single Cells/Lin-ve LD-ve/MHCII+ve/CD14-ve CD16-ve/CD123-/CD11c+ CD1c+ | Freq. of root`,
                PDC = `PBMCs/Single Cells/Single Cells/Lin-ve LD-ve/MHCII+ve/CD14-ve CD16-ve/CD11c- CD123+ | Freq. of root`) %>%
  tidyr::pivot_longer(-c(Run_ID, Sample_ID, Response), names_to = "Celltype", values_to = "Percentage_relative_PBMCs")

write.csv(celltype_percentage, celltype_percentage_csv)

sessionInfo()