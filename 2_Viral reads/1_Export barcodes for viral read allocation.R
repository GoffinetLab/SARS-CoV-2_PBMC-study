library(Seurat)
library(dplyr)
library(cowplot)
library(patchwork)
library(ggplot2)
library(pdftools)
library(ggpubr)
library(monocle)
library(wesanderson)
library(tibble)
library(ggplotify)
#------------------------------------------------------------------------------------


sars2.all <- readRDS("path/to/annotated_object.rds")
head(sars2.all@meta.data)

unique(sars2.all@meta.data$orig.ident.HTO_maxID)
barcode <- sars2.all@meta.data

######################################################################################
# code to get barcodes of all cells within donor within a treatment


## Donor1 #####
#MOCK
bc_d11_mock <- subset(barcode, barcode$orig.ident.HTO_maxID == "JK11_Mock")
bc_d11_mock <- as.data.frame(rownames(bc_d11_mock))
names(bc_d11_mock) <- "bcode"

#SARS1
bc_d11_sars1 <- subset(barcode, barcode$orig.ident.HTO_maxID == "JK11_SARS1")
bc_d11_sars1 <- as.data.frame(rownames(bc_d11_sars1))
names(bc_d11_sars1) <- "bcode"

#SARS2  
bc_d11_sars2 <- subset(barcode, barcode$orig.ident.HTO_maxID == "JK11_SARS2")
bc_d11_sars2 <- as.data.frame(rownames(bc_d11_sars2))
names(bc_d11_sars2) <- "bcode"

## Donor2 ####
#MOCK
bc_d12_mock <- subset(barcode, barcode$orig.ident.HTO_maxID == "JS12_Mock")
bc_d12_mock <- as.data.frame(rownames(bc_d12_mock))
names(bc_d12_mock) <- "bcode"

#SARS1
bc_d12_sars1 <- subset(barcode, barcode$orig.ident.HTO_maxID == "JS12_SARS1")
bc_d12_sars1 <- as.data.frame(rownames(bc_d12_sars1))
names(bc_d12_sars1) <- "bcode"

#SARS2  
bc_d12_sars2 <- subset(barcode, barcode$orig.ident.HTO_maxID == "JS12_SARS2")
bc_d12_sars2 <- as.data.frame(rownames(bc_d12_sars2))
names(bc_d12_sars2) <- "bcode"

##create list
df_list <- list("bc_d11_mock"= bc_d11_mock,
                "bc_d11_sars1"=  bc_d11_sars1, 
                "bc_d11_sars2"=  bc_d11_sars2,
                "bc_d12_mock"=  bc_d12_mock,
                "bc_d12_sars1"=  bc_d12_sars1,
                "bc_d12_sars2"=  bc_d12_sars2)



## Generate TSV tables with barcodes of cells corresponding to treatment (based on HTO identifier)

##export lists of barcodes
for (i in 1:length(df_list)){
  b <- df_list[[i]]
  b$bcode <- gsub("_", "-", b$bcode)
  #b <- b[b$bcode %in% viral,]
  #b <- as.vector(viral$)
 # b <- dplyr::data_frame(bcode= b)
  b$bcode <- paste0("CB:Z:", b$bcode)
  b$bcode <- gsub("-2", "-1", b$bcode)
  output.file <- file(paste0("./", names(df_list)[i], ".txt"), "wb")
  write.table(b, 
              file = output.file, quote = F, row.names = F, col.names = F, sep="\t" )
  close(output.file)
  
}























##################################################################################################
# code to get bcs for virus + cells found by CF
# deprecated approach!


#viral + and -

#viralreads <- read.csv("C:/Users/Dylan2/Google Drive/UNI WORK/2020 Internship/SARS2/viralreads.csv")
#donor11
viralreads <- read.csv("C:/Users/Dylan2/Google Drive/UNI WORK/2020 Internship/SARS2/viral/CG_JK_11_10x_AY_NC_viral.csv")
head(viralreads)
viralreads$Barcode <- gsub("-1", "-1", viralreads$Barcode)

#donor12
viralreads2 <- read.csv("C:/Users/Dylan2/Google Drive/UNI WORK/2020 Internship/SARS2/viral/CG_JK_12_10x_AY_NC_viral_cloupe.csv")
head(viralreads2)
viralreads2$Barcode <- gsub("-1", "-2", viralreads2$Barcode)


viral <- rbind(viralreads, viralreads2)
head(viral)

length(viral$viral)
#viral$Barcode <- gsub("-", "_", viral$Barcode)
viral <- as.vector(viral$Barcode)
viral

sars2.all[["viral"]] <- ifelse(rownames(sars2.all@meta.data) %in% viral, "plus", "minus")
unique(sars2.all@meta.data$viral)
table(sars2.all@meta.data$viral)
head(sars2.all@meta.data)
table(sars2.all@meta.data$viral, sars2.all$HTO_maxID)
table(sars2.all@meta.data$viral, sars2.all$orig.ident.HTO_maxID)
table(sars2.all@meta.data$viral, sars2.all$CellTypes)


#original viral align
#viral <- read.csv("C:/Users/Dylan2/Google Drive/UNI WORK/2020 Internship/SARS2/viralreads.csv")

#071220 revised viral align, first alignment w/ added donor12 reads
viral <- read.csv("C:/Users/Dylan2/Google Drive/UNI WORK/2020 Internship/SARS2/viral/CG_JK_11_10x_AY_NC_viral.csv")
viral
names(viral)
unique(viral$viral)
viral <- viral[viral$viral != "",]
nrow(viral)
##
bc_d11_mock
##


## D11 #####
  #MOCK
bc_d11_mock <- subset(barcode, barcode$orig.ident.HTO_maxID == "JK11_Mock")
bc_d11_mock <- as.data.frame(rownames(bc_d11_mock))
names(bc_d11_mock) <- "bcode"

  #SARS1
bc_d11_sars1 <- subset(barcode, barcode$orig.ident.HTO_maxID == "JK11_SARS1")
bc_d11_sars1 <- as.data.frame(rownames(bc_d11_sars1))
names(bc_d11_sars1) <- "bcode"

  #SARS2  
bc_d11_sars2 <- subset(barcode, barcode$orig.ident.HTO_maxID == "JK11_SARS2")
bc_d11_sars2 <- as.data.frame(rownames(bc_d11_sars2))
names(bc_d11_sars2) <- "bcode"

## D12 ####
  #MOCK
bc_d12_mock <- subset(barcode, barcode$orig.ident.HTO_maxID == "JS12_Mock")
bc_d12_mock <- as.data.frame(rownames(bc_d12_mock))
names(bc_d12_mock) <- "bcode"

  #SARS1
bc_d12_sars1 <- subset(barcode, barcode$orig.ident.HTO_maxID == "JS12_SARS1")
bc_d12_sars1 <- as.data.frame(rownames(bc_d12_sars1))
names(bc_d12_sars1) <- "bcode"

  #SARS2  
bc_d12_sars2 <- subset(barcode, barcode$orig.ident.HTO_maxID == "JS12_SARS2")
bc_d12_sars2 <- as.data.frame(rownames(bc_d12_sars2))
names(bc_d12_sars2) <- "bcode"

##create list
df_list <- list("bc_d11_mock"= bc_d11_mock,
                "bc_d11_sars1"=  bc_d11_sars1, 
                "bc_d11_sars2"=  bc_d11_sars2,
                "bc_d12_mock"=  bc_d12_mock,
                "bc_d12_sars1"=  bc_d12_sars1,
                "bc_d12_sars2"=  bc_d12_sars2)


b <- bc_d11_sars2
head(b)
viral
##export lists of barcodes
for (i in 1:length(df_list)){
  b <- df_list[[i]]
  b$bcode <- gsub("_", "-", b$bcode)
  b <- b[b$bcode %in% viral,]
  #b <- as.vector(viral$)
  b <- dplyr::data_frame(bcode= b)
  b$bcode <- paste0("CB:Z:", b$bcode)
  b$bcode <- gsub("-2", "-1", b$bcode)
  output.file <- file(paste0("C:/Users/Dylan2/Google Drive/UNI WORK/2020 Internship/SARS2/viral/", names(df_list)[i], ".txt"), "wb")
  write.table(b, 
              file = output.file, quote = F, row.names = F, col.names = F, sep="\t" )
  close(output.file)
  
}



b <- bc_d12_mock
head(b)
viral





















lapply(df_list, function(x){
  print(objname(x))
})



for (bc in list("bc_d11_mock"= bc_d11_mock,
                "bc_d11_sars1"=  bc_d11_sars1, 
                "bc_d11_sars2"=  bc_d11_sars2,
                "bc_d12_mock"=  bc_d12_mock,
                "bc_d12_sars1"=  bc_d12_sars1,
                "bc_d12_sars2"=  bc_d12_sars2)){
  
  b <- bc
  print(objname(bc))
  
  b$bcode <- gsub("_", "-", b$bcode)
  b$bcode <- paste0("CB:Z:", b$bcode)
  
  
  
  
  
  output.file <- file(paste0("C:/Users/Dylan2/Google Drive/UNI WORK/2020 Internship/SARS2/viral/", objname(b), ".txt"), "wb")
  write.table(b, 
              file = output.file, quote = F, row.names = F, col.names = F, sep="\t" )
  
  close(output.file)
  
}




for (bc in list(bc_d11_mock, bc_d11_sars1)){
  
  b <- bc
  print(objname(bc))
  
  b$bcode <- gsub("_", "-", b$bcode)
  b$bcode <- paste0("CB:Z:", b$bcode)
  
  
 
 
  
  output.file <- file(paste0("C:/Users/Dylan2/Google Drive/UNI WORK/2020 Internship/SARS2/viral/", objname(b), ".txt"), "wb")
  write.table(b, 
              file = output.file, quote = F, row.names = F, col.names = F, sep="\t" )
  
  close(output.file)

}

d <- list(bc_d11_mock, bc_d11_sars1)

for (b in c(bc_d11_mock)){
  
  b <- gsub("_", "-", b)
  b <- paste0("CB:Z:", b)
  b <- as.data.frame(b)
  
  output.file <- file(paste0("C:/Users/Dylan2/Google Drive/UNI WORK/2020 Internship/SARS2/viral/", objname(b), ".txt"), "wb")
  write.table(b, 
              file = output.file, quote = F, row.names = F, col.names = F, sep="\t" )
  
  close(output.file)
}

  

for bc in c(bc_d11_mock){
  e <- toString(bc)
}
  
 
e <- toString(bc_d11_mock)
e
d

barcode <- gsub("_", "-", barcode)
barcode <- paste0("CB:Z:", barcode)
barcode <- as.data.frame(barcode)

output.file <- file("C:/Users/Dylan2/Google Drive/UNI WORK/2020 Internship/SARS2/viral/barcode_filter.txt", "wb")
write.table(barcode, 
            file = output.file, quote = F, row.names = F, col.names = F, sep="\t" )

close(output.file)
head(barcode)

list(bc_d11_mock, bc_d11_sars1)
