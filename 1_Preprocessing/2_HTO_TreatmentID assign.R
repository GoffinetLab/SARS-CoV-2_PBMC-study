library(Seurat)


#### PROCESS DONOR 1 ####


#Basic setup
##Read in data
pbmc.umis <- Read10X("./CG_JK_11_10x/umi_count/", gene.column = 1)
row.names(pbmc.umis)<-c("Mock", "SARS-1", "SARS-2", "Unmapped")
pbmc.htos <- Read10X(data.dir = "./CG_JK_11_10x/read_count/", gene.column = 1)
row.names(pbmc.htos)<-c("Mock", "SARS-1", "SARS-2", "Unmapped")

joint.bcs <- intersect(colnames(pbmc.umis), colnames(pbmc.htos))
pbmc.umis <- pbmc.umis[, joint.bcs]
pbmc.htos <- as.matrix(pbmc.htos[, joint.bcs])
rownames(pbmc.htos)

## export HTO count csv
fwrite(as.data.table(pbmc.htos, keep.rownames="feature"), "./CG_JK_11_10x/CG_JK_11_10x_HTO_counts.csv")

##Setup Seurat object and add in the HTO data
pbmc.hashtag <- CreateSeuratObject(counts = pbmc.umis)
pbmc.hashtag <- NormalizeData(pbmc.hashtag)
pbmc.hashtag <- FindVariableFeatures(pbmc.hashtag, selection.method = "mean.var.plot")
pbmc.hashtag <- ScaleData(pbmc.hashtag, features = VariableFeatures(pbmc.hashtag))
#Adding HTO data as an independent assay
pbmc.hashtag[["HTO"]] <- CreateAssayObject(counts = pbmc.htos)
pbmc.hashtag <- NormalizeData(pbmc.hashtag, assay = "HTO", normalization.method = "CLR")
#Demultiplex cells based on HTO enrichment
pbmc.hashtag <- HTODemux(pbmc.hashtag, assay = "HTO", positive.quantile = 0.99)
#Visualize demultiplexing results
table(pbmc.hashtag$HTO_classification.global)
table(pbmc.hashtag$hash.ID)
table(pbmc.hashtag$HTO_classification.global)

##Visualize enrichment for selected HTOs with ridge plots
Idents(pbmc.hashtag) <- "HTO_maxID"
RidgePlot(pbmc.hashtag, assay = "HTO", features = rownames(pbmc.hashtag[["HTO"]])[1:3], ncol = 3, log = T)
##Visualize pairs of HTO signals to confirm mutual exclusivity in singlets
c("Mock", "SARS-1", "SARS-2", "Unmapped")
FeatureScatter(pbmc.hashtag, feature1 = "Mock", feature2 = "SARS-1", pt.size = .1)
FeatureScatter(pbmc.hashtag, feature1 = "SARS-2", feature2 = "SARS-1", pt.size = .1)
##Compare number of UMIs for singlets, doublets and negative cells
Idents(pbmc.hashtag) <- "HTO_classification.global"
table(pbmc.hashtag$nCount_RNA)
VlnPlot(pbmc.hashtag, features = "nCount_RNA", pt.size = 0.1, log = TRUE)

## export for processing in 10X Loupe Browser


cloupecord<-pbmc.hashtag$HTO_maxID
#names(cloupecord)<-paste0(names(cloupecord),"-1_2")
names(cloupecord)<-paste0(names(cloupecord),"-1_1")

cloupecord<-data.frame(Barcode=as.character(gsub("[1-9]_", "", names(cloupecord))),
                       HTO=cloupecord, 
                       stringsAsFactors = F)

cloupecord_real<-data.frame(Barcode=as.character(gsub("[1-9]_", "", names(seurat$stim))),
                            Suerat_sclusters=seurat$integrated_snn_res.0.2, 
                            stringsAsFactors = F)

cloupecord<-subset(cloupecord, Barcode %in% intersect((cloupecord_real$Barcode), (cloupecord$Barcode)))
names(cloupecord)[2]<-"HTO_11"

write.table(x = cloupecord, file = paste0("SP056_CsC_Sinit_HTO_CG_JK_11_10x", "_categories_.csv"), row.names = F, quote = F, sep = ",")
  

#### PROCESS DONOR 2 ####


#Basic setup
##Read in data
pbmc.umis <- Read10X("./CG_JK_12_10x/umi_count/", gene.column = 1)
row.names(pbmc.umis)<-c("Mock", "SARS-1", "SARS-2", "Unmapped")
pbmc.htos <- Read10X(data.dir = "./CG_JK_12_10x/read_count/", gene.column = 1)
row.names(pbmc.htos)<-c("Mock", "SARS-1", "SARS-2", "Unmapped")

joint.bcs <- intersect(colnames(pbmc.umis), colnames(pbmc.htos))
pbmc.umis <- pbmc.umis[, joint.bcs]
pbmc.htos <- as.matrix(pbmc.htos[, joint.bcs])
rownames(pbmc.htos)

##Setup Seurat object and add in the HTO data
pbmc.hashtag <- CreateSeuratObject(counts = pbmc.umis)
pbmc.hashtag <- NormalizeData(pbmc.hashtag)
pbmc.hashtag <- FindVariableFeatures(pbmc.hashtag, selection.method = "mean.var.plot")
pbmc.hashtag <- ScaleData(pbmc.hashtag, features = VariableFeatures(pbmc.hashtag))
#Adding HTO data as an independent assay
pbmc.hashtag[["HTO"]] <- CreateAssayObject(counts = pbmc.htos)
pbmc.hashtag <- NormalizeData(pbmc.hashtag, assay = "HTO", normalization.method = "CLR")
#Demultiplex cells based on HTO enrichment
pbmc.hashtag <- HTODemux(pbmc.hashtag, assay = "HTO", positive.quantile = 0.99)
#Visualize demultiplexing results
table(pbmc.hashtag$HTO_classification.global)
table(pbmc.hashtag$hash.ID)
table(pbmc.hashtag$HTO_classification.global)

##Visualize enrichment for selected HTOs with ridge plots
Idents(pbmc.hashtag) <- "HTO_maxID"
RidgePlot(pbmc.hashtag, assay = "HTO", features = rownames(pbmc.hashtag[["HTO"]])[1:3], ncol = 3, log = T)
##Visualize pairs of HTO signals to confirm mutual exclusivity in singlets
c("Mock", "SARS-1", "SARS-2", "Unmapped")
FeatureScatter(pbmc.hashtag, feature1 = "Mock", feature2 = "SARS-1", pt.size = .1)
FeatureScatter(pbmc.hashtag, feature1 = "SARS-2", feature2 = "SARS-1", pt.size = .1)
##Compare number of UMIs for singlets, doublets and negative cells
Idents(pbmc.hashtag) <- "HTO_classification.global"
table(pbmc.hashtag$nCount_RNA)
VlnPlot(pbmc.hashtag, features = "nCount_RNA", pt.size = 0.1, log = TRUE)


cloupecord<-pbmc.hashtag$HTO_maxID
names(cloupecord)<-paste0(names(cloupecord),"-1_2")
#names(cloupecord)<-paste0(names(cloupecord),"-1_1")

cloupecord<-data.frame(Barcode=as.character(gsub("[1-9]_", "", names(cloupecord))),
                       HTO=cloupecord, 
                       stringsAsFactors = F)

cloupecord_real<-data.frame(Barcode=as.character(gsub("[1-9]_", "", names(seurat$stim))),
                            Suerat_sclusters=seurat$integrated_snn_res.0.2, 
                            stringsAsFactors = F)

cloupecord<-subset(cloupecord, Barcode %in% intersect((cloupecord_real$Barcode), (cloupecord$Barcode)))
names(cloupecord)[2]<-"HTO_12"

write.table(x = cloupecord, file = paste0( "SP056_CsC_Sinit_HTO_CG_JK_12_10x", "_categories_.csv"), row.names = F, quote = F, sep = ",")

## export HTO count csv
fwrite(as.data.table(pbmc.htos, keep.rownames="feature"), "./CG_JK_11_10x/CG_JK_12_10x_HTO_counts.csv")
