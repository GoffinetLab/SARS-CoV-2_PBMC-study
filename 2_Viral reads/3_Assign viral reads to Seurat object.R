
library(tidyverse)
library(Seurat)
library(cowplot)
library(patchwork)

#----------------------------------------------------------------------------

## Assign viral read counts to cells 

## Donor11
  #SARS2 reads
d11_sars2_sars2 <- read.table("./bc_viral_SARS2_CG_JK_11_sars2.txt")
colnames(d11_sars2_sars2) <- c("barcode", "umi", "tag")
d11_sars2_sars2$donor <- "d11" 
d11_sars2_sars2$virus <- "S2"
d11_s2_s2_bc <- data.frame(d11_sars2_sars2[, c(1,4,5)])
d11_s2_s2_bc <- unique(d11_s2_s2_bc)


  #SARS1 reads
d11_sars1_sars1 <- read.table("./bc_viral_SARS1_CG_JK_11_sars1.txt")
colnames(d11_sars1_sars1) <- c("barcode", "umi", "tag")
d11_sars1_sars1$donor <- "d11" 
d11_sars1_sars1$virus <- "S1"
d11_s1_s1_bc <- data.frame(d11_sars1_sars1[, c(1,4,5)])
d11_s1_s1_bc <- unique(d11_s1_s1_bc)



## D12
  #SARS2 reads
d12_sars2_sars2 <- read.table("./bc_viral_SARS2_CG_JK_12_sars2.txt")
colnames(d12_sars2_sars2) <- c("barcode", "umi", "tag")
d12_sars2_sars2$barcode <- gsub("1", "2", d12_sars2_sars2$barcode)
d12_sars2_sars2$donor <- "d12" 
d12_sars2_sars2$virus <- "S2"
d12_s2_s2_bc <- data.frame(d12_sars2_sars2[, c(1,4,5)])
d12_s2_s2_bc <- unique(d12_s2_s2_bc)

  #SARS1 reads
d12_sars1_sars1 <- read.table("./bc_viral_SARS1_CG_JK_12_sars1.txt")
colnames(d12_sars1_sars1) <- c("barcode", "umi", "tag")
d12_sars1_sars1$barcode <- gsub("1", "2", d12_sars1_sars1$barcode)
d12_sars1_sars1$donor <- "d12" 
d12_sars1_sars1$virus <- "S1"
d12_s1_s1_bc <- data.frame(d12_sars1_sars1[, c(1,4,5)])
d12_s1_s1_bc <- unique(d12_s1_s1_bc)


## get reads per cell
reads_cell <- rbind(d11_sars1_sars1, d11_sars2_sars2, d12_sars1_sars1, d12_sars2_sars2) %>%
  group_by(barcode) %>%
  summarise(viral_count = n_distinct(umi), donor=donor, virus=virus ) %>%
  unique() 
reads_cell$barcode <- gsub("-", "_", reads_cell$barcode)


## read in annotated Seurat object

sars2.all <- readRDS("path/to/annotated/Seurat_object.rds")


## assign status as viral RNA positive or negative to cells in Seurat object

sars2.all@meta.data$viral <- ifelse(sars2.all@meta.data$barcode %in% reads_cell$barcode, "pos", "neg")
#viral status and counts
viral_count <- plyr::join(sars2.all@meta.data, reads_cell[, c("barcode", "viral_count", "virus")])
viral_count[is.na(viral_count)] <- 0
sars2.all@meta.data$viral_count <- viral_count$viral_count

sars2.all@meta.data$virus <- ifelse(sars2.all$viral == "pos" & sars2.all$HTO_maxID == "SARS2", "S2",
                                    ifelse(sars2.all$viral == "pos" & sars2.all$HTO_maxID == "SARS1", "S1", "uninfected"))



