library(Seurat)
library(future)
library(openxlsx)
library(cowplot)


remove(namingvars); namingvars<-list()
setwd("/wd")

namingvars$group<-c("human_10X")
namingvars$indiv<-c("CG_JK_11_10x",
                    "CG_JK_12_10x")

cond1.data <- Read10X_h5("./CG_JK_11_10x/outs/filtered_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE)
cond2.data <- Read10X_h5("./CG_JK_12_10x/outs/filtered_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE)


## Standard Seurat Preprocessing

setupMultiobject<-function(x,y,z) {
  cond <- CreateSeuratObject(counts = z, project = y, min.cells = 3)
  cond$stim <- x
  cond <- subset(cond, subset = nFeature_RNA > 100)  
  cond <- NormalizeData(cond, verbose = FALSE)
  cond <- FindVariableFeatures(cond, selection.method = "vst", nfeatures = 2000)
  cond
}


cond1<-setupMultiobject(namingvars$indiv[1], namingvars$group, cond1.data)
cond2<-setupMultiobject(namingvars$indiv[2], namingvars$group, cond2.data)

conditions<-list(cond1, cond2)
remove(cond1.data, cond2.data, cond3.data)

options(future.globals.maxSize = 4 * 1000 * 1024^2) # x * 1 GB

## integrate data to remove batch effects

immune.anchors <- FindIntegrationAnchors(object.list = conditions, dims = 1:20)
remove(conditions)
immune.combined <- IntegrateData(anchorset = immune.anchors, dims = 1:20)
remove(immune.anchors)

## Perform an integrated analysis


# Run the standard workflow for visualization and clustering
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
immune.combined <- RunTSNE(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindClusters(immune.combined, resolution = .2)

plot_grid(
  DimPlot(immune.combined, reduction = "umap", label = T) + NoAxes() + coord_fixed() +ggtitle(""),
  DimPlot(immune.combined, reduction = "umap", group.by = "stim", label = F, pt.size = .1) + NoAxes() + coord_fixed()+ggtitle("")
)
plot_grid(DimPlot(immune.combined, reduction = "umap", split.by = "stim", label = T) + NoAxes() + coord_fixed(), ncol=1)

## perform alternative cluster similarity spectrum processing

library(simspec)
seurat<-immune.combined
unique(seurat$stim)
DefaultAssay(seurat)<-"integrated"
seurat <- cluster_sim_spectrum(seurat, label_tag = "stim", use_scale = T)
seurat <- RunUMAP(seurat, reduction = "css", dims = 1:ncol(Embeddings(seurat, "css")), reduction.name = "css", reduction.key = "CSS")

plot_grid(
  DimPlot(seurat, reduction = "css", label = T) + NoAxes() + coord_fixed() +ggtitle(""),
  DimPlot(seurat, reduction = "css", group.by = "stim", label = F, pt.size = .1) + NoAxes() + coord_fixed()+ggtitle("")
)
plot_grid(DimPlot(seurat, reduction = "css", split.by = "stim", label = T) + NoAxes() + coord_fixed(), ncol=1)

DefaultAssay(seurat) <- "RNA"

## define cell cycle scores

seurat<-CellCycleScoring(seurat,   g2m.features = cc.genes$g2m.genes, s.features = cc.genes$s.genes, set.ident = T)
DefaultAssay(seurat) <- "integrated"
DimPlot(seurat, reduction = "umap", group.by = "Phase")

plot_grid(
  DimPlot(seurat, reduction = "umap", group.by = "Phase", label = F) + NoAxes() + coord_fixed() +ggtitle(""),
  DimPlot(seurat, reduction = "css", group.by = "Phase", label = F) + NoAxes() + coord_fixed()+ggtitle("")
)

### cloupe generation for processing with 10X Loupe Browser

cloupecord<-data.frame(Barcode=as.character(gsub("[1-9]_", "", names(seurat$stim))),
                       UMAP_1=seurat@reductions$umap@cell.embeddings[,1], 
                       UMAP_2=seurat@reductions$umap@cell.embeddings[,2],
                       stringsAsFactors = F)
write.table(x = cloupecord, file = paste0("./", "SP056_CR3int_Sinit_umap", "_projection.csv"), row.names = F, quote = F, sep = ",")

cloupecord<-data.frame(Barcode=as.character(gsub("[1-9]_", "", names(seurat$stim))),
                       CSS_1=seurat@reductions$css@cell.embeddings[,1], 
                       CSS_2=seurat@reductions$css@cell.embeddings[,2],
                       stringsAsFactors = F)
write.table(x = cloupecord, file = paste0("./", "SP056_CR3int_Sinit_CSSumap", "_projection.csv"), row.names = F, quote = F, sep = ",")

cloupecord<-data.frame(Barcode=as.character(gsub("[1-9]_", "", names(seurat$stim))),
                       Suerat_sclusters=seurat$integrated_snn_res.0.2, 
                       stringsAsFactors = F)

write.table(x = cloupecord, file = paste0("./", "SP056_CR3int_Sinit", "_categories.csv"), row.names = F, quote = F, sep = ",")

