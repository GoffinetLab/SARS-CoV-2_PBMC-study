library(Seurat)
library(dplyr)
library(cowplot)
library(patchwork)
library(xlsx)
library(ggplot2)
library(pdftools)
library(ggpubr)
library(monocle)
library(DESeq2)
library(fgsea)
library(wesanderson)
library(msigdbr)
library(tibble)
library(MAST)
library(readtext)

#------------------------------------------------------------------------

## read in Reactome IFN signalling module genes

ifn_sig <- read.table("C./Reactome Interferon Signalling [R-HSA-913531].tsv",
                      sep = "\t", header = T)
ifn_sig$gene <- gsub("UniProt:.* (.*)$", "\\1", ifn_sig$MoleculeName)
ifn_sig <- ifn_sig$gene
ifn_sig <- list(ifn_sig)

## Define trajectory plotting functions

floor_dec <- function(x, level=1) round(x - 5*10^(-level-1), level)
ceiling_dec <- function(x, level=1) round(x + 5*10^(-level-1), level)


## functions to visualise pseudotime trajectories

# colour by Pseudotime
ptime.tjct <- function(cds,rev=F, max=20){
  
  print("Max pseudoval is:")
  print(max(cds$Pseudotime))
  print("Min pseudoval is:")
  print(min(cds$Pseudotime))
  
  pal <- wes_palette("Zissou1", 20, type = "continuous")
  if (rev==T){
    pal <- rev(wes_palette("Zissou1", 20, type = "continuous"))
  }
  
  name <- deparse(substitute(cds))
  name <- gsub("sars2.(.*)cds", "\\1", name)
  
  #without legend
  tjct <-plot_cell_trajectory(cds, color_by = "Pseudotime", cell_size=0.7)+
    scale_color_gradientn(colours=pal, limit =c(0, max))+
    labs(color= "Pseudotime")+
    theme(axis.title = element_blank())+
    theme(axis.text.x= element_text( face='bold'))+
    theme(axis.text.y= element_text( face='bold'))+
    theme(legend.position = "none")+
    NoAxes()
  
  ggsave(paste0(path,name,"_pseudo_pseudotime.pdf"), plot=tjct, width = 5, height =3.5, dpi=200)
  ggsave(paste0(path,name,"_pseudo_pseudotime.eps"), plot=tjct, width = 5, height =3.5, dpi=200)
  #extract legend separately
  tjct <- plot_cell_trajectory(cds, color_by = "Pseudotime", cell_size=0.7)+
    scale_color_gradientn(colours=pal, limit =c(0, max))+
    labs(color= "Pseudotime")+
    theme(axis.title = element_blank())+
    theme(axis.text.x= element_text( face='bold'))+
    theme(axis.text.y= element_text( face='bold'))+
    theme(legend.position = "bottom")+
    NoAxes()
  
  leg <- get_legend(tjct )
  leg <- as_ggplot(leg)
  leg
  ggsave(paste0(path,name,"_pseudo_leg.eps"), plot = leg, width=2.5, height=1, dpi=200)
  ggsave(paste0(path,name,"_pseudo_leg.pdf"), plot = leg, width=2.5, height=1, dpi=200)
  
}


## colour by subcluster
subc.tjct <- function(cds,rev=F){
  
  
  name <- deparse(substitute(cds))
  name <- gsub("sars2.(.*)cds", "\\1", name)
  
  
  #without legend
  tjct <- plot_cell_trajectory(cds, color_by = "seurat_clusters_int", cell_size=0.7)+
    scale_color_manual(labels=c('0', '1', '2', '3', 'undefined'), values=c("#E6CC80", "#E67F63", "#5B5F89", "#A9A693", "grey" ))+
    labs(color= "")+
    theme(axis.title = element_blank())+
    theme(axis.text.x= element_text( face='bold'))+
    theme(axis.text.y= element_text( face='bold'))+
    theme(legend.text=element_text(size=15, face="bold"))+
    theme(legend.position = "none")+
    NoAxes()
  
  ggsave(paste0(path,name,"_pseudo_subclust.pdf"), plot=tjct, width = 5, height =3.5, dpi=200)
  ggsave(paste0(path,name,"_pseudo_subclust.eps"), plot=tjct, width = 5, height =3.5, dpi=200)
  
  #extract legend separately
  tjct <- plot_cell_trajectory(cds, color_by = "seurat_clusters_int", cell_size=0.7)+
    scale_color_manual(labels=c('0', '1', '2', '3', 'undefined'), values=c("#E6CC80", "#E67F63", "#5B5F89", "#A9A693", "grey" ))+
    labs(color= "")+
    theme(axis.title = element_blank())+
    theme(axis.text.x= element_text( face='bold'))+
    theme(axis.text.y= element_text( face='bold'))+
    theme(legend.text=element_text(size=15, face="bold"))+
    theme(legend.position = "bottom")+
    NoAxes()
  
  leg <- get_legend(tjct )
  leg <- as_ggplot(leg)
  
  ggsave(paste0(path,name,"_subclust_leg.eps"), plot = leg, width=2.5, height=1, dpi=200)
  ggsave(paste0(path,name,"_subclust_leg.pdf"), plot = leg, width=2.5, height=1, dpi=200)
  
}

## colour by IFN module score

ifnmod.tjct <- function(cds, min=-0.1, max=0.5){
  
  
  high <- rgb(9,153,99, max=255)
  med <- "white"
  low <- rgb(159,4,77, max=255)
  
  print("Min ifn mod score is:")
  print(min(cds$ifn_sig1))
  print("Max ifn mod score is:")
  print(max(cds$ifn_sig1))
  
  name <- deparse(substitute(cds))
  name <- gsub("sars2.(.*)cds", "\\1", name)
  
  
  #without legend
  tjct <-plot_cell_trajectory(cds, color_by = "ifn_sig1", cell_size=0.7)+
    scale_color_gradient(high= high, low=low, limit=c(min, max)) +
    theme(axis.title = element_blank())+
    theme(axis.text.x= element_text( face='bold'))+
    theme(axis.text.y= element_text( face='bold'))+
    theme(legend.position = "none")+
    NoAxes()
  ggsave(paste0(path,name,"_pseudo_ifnmod.pdf"), plot=tjct, width = 5, height =3.5, dpi=200)
  ggsave(paste0(path,name,"_pseudo_ifnmod.eps"), plot=tjct, width = 5, height =3.5, dpi=200)
  #extract legend separately
  tjct <-plot_cell_trajectory(cds, color_by = "ifn_sig1", cell_size=0.7)+
    scale_color_gradient(high= high, low=low, limit=c(min, max)) +
    theme(axis.title = element_blank())+
    theme(axis.text.x= element_text( face='bold'))+
    theme(axis.text.y= element_text( face='bold'))+
    theme(legend.position = "bottom")+
    NoAxes()
  
  leg <- get_legend(tjct)
  leg <- as_ggplot(leg)
  leg
  ggsave(paste0(path,name,"_ifnmod_leg.eps"), plot = leg, width=2.5, height=1, dpi=200)
  ggsave(paste0(path,name,"_ifn_leg.pdf"), plot = leg, width=2.5, height=1, dpi=200)
  
  
}


## colour by treatment
treat.tjct <- function(cds){
  
  
  name <- deparse(substitute(cds))
  name <- gsub("sars2.(.*)cds", "\\1", name)
  
  
  #without legend
  tjct <-plot_cell_trajectory(cds, color_by = "HTO_maxID", cell_size=0.7)+
    scale_color_manual(breaks = c("Mock", "SARS1", "SARS2"), values=c("#626972", "#A0D2E1", "#649A7D" ))+
    labs(color= "")+
    theme(axis.title = element_blank())+
    theme(axis.text.x= element_text( face='bold'))+
    theme(axis.text.y= element_text( face='bold'))+
    theme(legend.text=element_text(size=15, face="bold"))+
    theme(legend.position = "none")+
    NoAxes()
  
  
  ggsave(paste0(path,name,"_pseudo_treat.pdf"), plot=tjct, width = 5, height =3.5, dpi=200)
  ggsave(paste0(path,name,"_pseudo_treat.eps"), plot=tjct, width = 5, height =3.5, dpi=200)
  
  #extract legend separately
  tjct <-plot_cell_trajectory(cds, color_by = "HTO_maxID", cell_size=0.7)+
    scale_color_manual(breaks = c("Mock", "SARS1", "SARS2"), values=c("#626972", "#A0D2E1", "#649A7D" ))+
    labs(color= "")+
    theme(axis.title = element_blank())+
    theme(axis.text.x= element_text( face='bold'))+
    theme(axis.text.y= element_text( face='bold'))+
    theme(legend.text=element_text(size=15, face="bold"))+
    theme(legend.position = "bottom")+
    NoAxes()
  
  leg <- get_legend(tjct )
  leg <- as_ggplot(leg)
  
  ggsave(paste0(path,name,"_treat_leg.eps"), plot = leg, width=3.5, height=1, dpi=200)
  ggsave(paste0(path,name,"_treat_leg.pdf"), plot = leg, width=3.5, height=1, dpi=200)
  
}



#------------------------------------------------------------------------

### read in fully annotated Seurat object (including viral reads)

sars2.all <- readRDS("path/to/annotated/Seurat_object_with_viral_read_counts.rds")



#### SUBCLUSTER by cell type #####

## Monocytes -------------------------------------------------------------


sars2.mono <-  subset(sars2.all, CellTypes == "Monocyte")

DefaultAssay(sars2.mono)<-"RNA"

mono.list <- SplitObject(sars2.mono, split.by = "orig.ident")

for (i in 1:length(mono.list)) {
  mono.list[[i]] <- SCTransform(mono.list[[i]], verbose = FALSE)
}
options(future.globals.maxSize =1 * 4000 * 1024^2)
pancreas.features <- SelectIntegrationFeatures(object.list = mono.list, nfeatures = 2000)


mono.list <- PrepSCTIntegration(object.list = mono.list, anchor.features = pancreas.features, 
                                verbose = FALSE)
pancreas.anchors <- FindIntegrationAnchors(object.list = mono.list, normalization.method = "SCT", 
                                           anchor.features = pancreas.features, verbose = FALSE)
pancreas.integrated <- IntegrateData(anchorset = pancreas.anchors, normalization.method = "SCT", verbose = FALSE)

sars2.mono <- RunPCA(pancreas.integrated, npcs = 10, verbose = FALSE)
sars2.mono <- RunUMAP(sars2.mono, dims = 1:10, reduction.name = "umapsub", reduction.key = "UMAPSUB")


sars2.mono <- FindNeighbors(sars2.mono, reduction = "pca", dims = 1:10)
sars2.mono <- FindClusters(sars2.mono, resolution = 0.3)
sars2.mono$seurat_clusters_int <- sars2.mono$seurat_clusters

rm(mono.list)
##UMAPS
clust_umap <- DimPlot(sars2.mono,label = T , reduction = "umapsub", group.by="seurat_clusters_int") + NoAxes() + coord_fixed() 
clust_umap
ggsave(paste0(path, "mono_umap_cluster.png"), plot= clust_umap)

infect_umap <- DimPlot(sars2.mono,label = T , reduction = "umapsub", group.by="virus") + NoAxes() + coord_fixed()
infect_umap

donor_umap <- DimPlot(sars2.mono,label = T , reduction = "umapsub", group.by="orig.ident") + NoAxes() + coord_fixed()
donor_umap

treat_umap<- DimPlot(sars2.mono,label = T , reduction = "umapsub", group.by = "HTO_maxID") + NoAxes() + coord_fixed() 
treat_umap
ggsave(paste0(path, "mono_umap_treat.png"), plot= treat_umap)

#add IFN signalling module score
sars2.mono <- AddModuleScore(object = sars2.mono, features = ifn_sig, ctrl = 100, name = 'ifn_sig', assay="RNA")
head(sars2.mono@meta.data)

#### MONOCYTE CDS ####


## IMPORT TO MONOCLE #

data <- as(as.matrix(sars2.mono@assays$RNA@counts), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = sars2.mono@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)

#Construct monocle cds
sars2.monocds <- newCellDataSet(data,
                                phenoData = pd,
                                featureData = fd,
                                lowerDetectionLimit = 0.5,
                                expressionFamily = negbinomial.size())



#estimate size factors and dispersion
sars2.monocds <- estimateSizeFactors(sars2.monocds)
sars2.monocds <- estimateDispersions(sars2.monocds)
#remove genes not expressed in at least 10 cells
sars2.monocds <- detectGenes(sars2.monocds, min_expr = 0.1)
print(head(fData(sars2.monocds)))
expressed_genes <- row.names(subset(fData(sars2.monocds),
                                    num_cells_expressed >= 10))




#determine ordering genes - diff expressed between infections


#Idents(sars2.mono) <- "HTO_maxID"
#ordering_genes <- FindAllMarkers(sars2.mono)

Idents(sars2.mono) <- "seurat_clusters_int"
unique(sars2.mono$seurat_clusters_int)
ordering_genes <- FindAllMarkers(sars2.mono)

ordering_genes <- ordering_genes[ordering_genes$p_val_adj < 0.05,]
ordering_genes$log2 <- log2(exp(ordering_genes$avg_logFC))

ordering_genes <- as.vector(ordering_genes$gene)

#set ordering genes
sars2.monocds <- setOrderingFilter(sars2.monocds, ordering_genes)
#reduce dimensions
sars2.monocds <- reduceDimension(sars2.monocds, max_components = 2,
                                 method = 'DDRTree')
#order cells along trajectory
sars2.monocds <- orderCells(sars2.monocds )
warnings()

#check root state
plot_cell_trajectory(sars2.monocds, color_by = "State") +
  facet_wrap(~State, nrow = 1)

plot_cell_trajectory(sars2.monocds, color_by = "HTO_maxID") +
  facet_wrap(~State, nrow = 1)


sars2.monocds <- orderCells(sars2.monocds, root_state = 6)

###MONOCYTE PSEUDOTIME ####


#plot pseudotime
ptime.tjct(sars2.monocds, max=18, rev=F)
#plot subcluster
subc.tjct(sars2.monocds)
#plot ifnmod
ifnmod.tjct(sars2.monocds)
#plot treatment
treat.tjct(sars2.monocds)
#plot virus pos cells
#### PSEUDOTIME VIRUS ####

unique(sars2.mono$virus)

mono.virus.tjct_pos <-plot_cell_trajectory(sars2.monocds, color_by ="virus", cell_size = 1, order=sars2.monocds$order)+
  scale_color_manual(labels=c('S1', 'S2', "uninfected"), values=c('S1' ="blue", 'S2' ="red", "uninfected" = NA))+
  labs(color= "")+
  theme(axis.title = element_blank())+
  theme(axis.text.x= element_text( face='bold'))+
  theme(axis.text.y= element_text( face='bold'))+
  theme(legend.text=element_text(size=15, face="bold"))+
  theme(legend.position = "none")+
  #scale_size_manual(labels=c('plus', 'minus'), values=c('plus'=10, 'minus'=0.2))+
  NoAxes()
mono.virus.tjct_pos

## B CELLS -----------------------------------------------------------------------------------------------------------------

#subcluster

sars2.bcell <-  subset(sars2.all, CellTypes == "B cell")

DefaultAssay(sars2.bcell)<-"RNA"

bcell.list <- SplitObject(sars2.bcell, split.by = "orig.ident")

for (i in 1:length(bcell.list)) {
  bcell.list[[i]] <- SCTransform(bcell.list[[i]], verbose = FALSE)
}
options(future.globals.maxSize =1 * 4000 * 1024^2)
pancreas.features <- SelectIntegrationFeatures(object.list = bcell.list, nfeatures = 2000)


bcell.list <- PrepSCTIntegration(object.list = bcell.list, anchor.features = pancreas.features, 
                                 verbose = FALSE)
pancreas.anchors <- FindIntegrationAnchors(object.list = bcell.list, normalization.method = "SCT", 
                                           anchor.features = pancreas.features, verbose = FALSE)
pancreas.integrated <- IntegrateData(anchorset = pancreas.anchors, normalization.method = "SCT", verbose = FALSE)

sars2.bcell <- RunPCA(pancreas.integrated, npcs = 10, verbose = FALSE)
sars2.bcell <- RunUMAP(sars2.bcell, dims = 1:10, reduction.name = "umapsub", reduction.key = "UMAPSUB")


sars2.bcell <- FindNeighbors(sars2.bcell, reduction = "pca", dims = 1:10)
sars2.bcell <- FindClusters(sars2.bcell, resolution = 0.3)
sars2.bcell$seurat_clusters_int <- sars2.bcell$seurat_clusters

rm(bcell.list)
##UMAPS
clust_umap <- DimPlot(sars2.bcell,label = T , reduction = "umapsub", group.by="seurat_clusters_int") + NoAxes() + coord_fixed() 
clust_umap
ggsave(paste0(path, "bcell_umap_cluster.png"), plot= clust_umap)

treat_umap<- DimPlot(sars2.bcell,label = T , reduction = "umapsub", group.by = "HTO_maxID") + NoAxes() + coord_fixed() 
ggsave(paste0(path, "bcell_umap_treat.png"), plot= treat_umap)


#add IFN signalling module score
sars2.bcell <- AddModuleScore(object = sars2.bcell, features = ifn_sig, ctrl = 100, name = 'ifn_sig', assay="RNA")

#### BCELL CDS ####


## IMPORT TO MONOCLE #

data <- as(as.matrix(sars2.bcell@assays$RNA@counts), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = sars2.bcell@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)

#Construct monocle cds
sars2.bcellcds <- newCellDataSet(data,
                                 phenoData = pd,
                                 featureData = fd,
                                 lowerDetectionLimit = 0.5,
                                 expressionFamily = negbinomial.size())


#estimate size factors and dispersion
sars2.bcellcds <- estimateSizeFactors(sars2.bcellcds)
sars2.bcellcds <- estimateDispersions(sars2.bcellcds)
#remove genes not expressed in at least 10 cells
sars2.bcellcds <- detectGenes(sars2.bcellcds, min_expr = 0.1)
print(head(fData(sars2.bcellcds)))
expressed_genes <- row.names(subset(fData(sars2.bcellcds),
                                    num_cells_expressed >= 10))




#determine ordering genes - diff expressed between infections
head(sars2.bcell@meta.data)



#Idents(sars2.bcell) <- "HTO_maxID"
#ordering_genes <- FindAllMarkers(sars2.bcell)

Idents(sars2.bcell) <- "seurat_clusters_int"
ordering_genes <- FindAllMarkers(sars2.bcell)

ordering_genes <- ordering_genes[ordering_genes$p_val_adj < 0.05,]
ordering_genes$log2 <- log2(exp(ordering_genes$avg_logFC))
ordering_genes <- ordering_genes[abs(ordering_genes$log2) > 2,]

ordering_genes <- as.vector(ordering_genes$gene)

#set ordering genes
sars2.bcellcds <- setOrderingFilter(sars2.bcellcds, ordering_genes)
#reduce dimensions
sars2.bcellcds <- reduceDimension(sars2.bcellcds, max_components = 2,
                                  method = 'DDRTree')
#order cells along trajectory
sars2.bcellcds <- orderCells(sars2.bcellcds )
warnings()

#check root state
plot_cell_trajectory(sars2.bcellcds, color_by = "HTO_maxID")

plot_cell_trajectory(sars2.bcellcds, color_by = "State") +
  facet_wrap(~State, nrow = 1)

cds <- sars2.bcellcds

##determine which state has most cells from mock
GM_state <- function(cds){
  if (length(unique(pData(cds)$State)) > 1){
    T0_counts <- table(pData(cds)$State, pData(cds)$HTO_maxID)[,"Mock"]
    return(as.numeric(names(T0_counts)[which
                                       (T0_counts == max(T0_counts))]))
  } else {
    return (1)
  }
}

root_state = GM_state(sars2.bcellcds)
root_state
sars2.bcellcds <- orderCells(sars2.bcellcds, root_state = GM_state(sars2.bcellcds))
plot_cell_trajectory(sars2.bcellcds, color_by = "Pseudotime")

###BCELL PSEUDOTIME ####

#plot pseudotime
ptime.tjct(sars2.bcellcds, max=18, rev=F)
#plot ifnmod
ifnmod.tjct(sars2.bcellcds)
#plot treatment
treat.tjct(sars2.bcellcds)



## NK CELLS -----------------------------------------------------------------------------------------------------------------

#subcluster

sars2.nk <-  subset(sars2.all, CellTypes == "NK + CD8+ T cell")

DefaultAssay(sars2.nk)<-"RNA"

nk.list <- SplitObject(sars2.nk, split.by = "orig.ident")

for (i in 1:length(nk.list)) {
  nk.list[[i]] <- SCTransform(nk.list[[i]], verbose = FALSE)
}
options(future.globals.maxSize =1 * 4000 * 1024^2)
pancreas.features <- SelectIntegrationFeatures(object.list = nk.list, nfeatures = 2000)


nk.list <- PrepSCTIntegration(object.list = nk.list, anchor.features = pancreas.features, 
                              verbose = FALSE)
pancreas.anchors <- FindIntegrationAnchors(object.list = nk.list, normalization.method = "SCT", 
                                           anchor.features = pancreas.features, verbose = FALSE)
pancreas.integrated <- IntegrateData(anchorset = pancreas.anchors, normalization.method = "SCT", verbose = FALSE)

sars2.nk <- RunPCA(pancreas.integrated, npcs = 10, verbose = FALSE)
sars2.nk <- RunUMAP(sars2.nk, dims = 1:10, reduction.name = "umapsub", reduction.key = "UMAPSUB")


sars2.nk <- FindNeighbors(sars2.nk, reduction = "pca", dims = 1:10)
sars2.nk <- FindClusters(sars2.nk, resolution = 0.3)
sars2.nk$seurat_clusters_int <- sars2.nk$seurat_clusters

##UMAPS
clust_umap <- DimPlot(sars2.nk,label = T , reduction = "umapsub", group.by="seurat_clusters_int") + NoAxes() + coord_fixed() 
clust_umap
ggsave(paste0(path, "nk_umap_cluster.png"), plot= clust_umap)

treat_umap<- DimPlot(sars2.nk,label = T , reduction = "umapsub", group.by = "HTO_maxID") + NoAxes() + coord_fixed() 
treat_umap
ggsave(paste0(path, "nk_umap_treat.png"), plot= treat_umap)


##add viral reads and ifn mod

#assign viral read status
sars2.nk[["viral"]] <- ifelse(rownames(sars2.nk@meta.data) %in% viralreads, "plus", "minus")
table(sars2.nk@meta.data$viral)

#add IFN signalling module score
sars2.nk <- AddModuleScore(object = sars2.nk, features = ifn_sig, ctrl = 100, name = 'ifn_sig', assay="RNA")


#### NK CDS ####


## IMPORT TO MONOCLE #

data <- as(as.matrix(sars2.nk@assays$RNA@counts), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = sars2.nk@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)

#Construct nkcle cds
sars2.nkcds <- newCellDataSet(data,
                              phenoData = pd,
                              featureData = fd,
                              lowerDetectionLimit = 0.5,
                              expressionFamily = negbinomial.size())


#estimate size factors and dispersion
sars2.nkcds <- estimateSizeFactors(sars2.nkcds)
sars2.nkcds <- estimateDispersions(sars2.nkcds)
#remove genes not expressed in at least 10 cells
sars2.nkcds <- detectGenes(sars2.nkcds, min_expr = 0.1)
print(head(fData(sars2.nkcds)))
expressed_genes <- row.names(subset(fData(sars2.nkcds),
                                    num_cells_expressed >= 10))

#determine ordering genes - diff expressed between infections

Idents(sars2.nk) <- "seurat_clusters_int"
ordering_genes <- FindAllMarkers(sars2.nk)

ordering_genes <- ordering_genes[ordering_genes$p_val_adj < 0.05,]
ordering_genes$log2 <- log2(exp(ordering_genes$avg_logFC))
ordering_genes <- ordering_genes[abs(ordering_genes$log2) > 2,]

ordering_genes <- as.vector(ordering_genes$gene)

#set ordering genes
sars2.nkcds <- setOrderingFilter(sars2.nkcds, ordering_genes)
#reduce dimensions
sars2.nkcds <- reduceDimension(sars2.nkcds, max_components = 2,
                               method = 'DDRTree')
#order cells along trajectory
sars2.nkcds <- orderCells(sars2.nkcds )

#check root state
plot_cell_trajectory(sars2.nkcds, color_by = "HTO_maxID")

plot_cell_trajectory(sars2.nkcds, color_by = "State") +
  facet_wrap(~State, nrow = 1)


##determine which state has most cells from mock
GM_state <- function(cds){
  if (length(unique(pData(cds)$State)) > 1){
    T0_counts <- table(pData(cds)$State, pData(cds)$HTO_maxID)[,"Mock"]
    return(as.numeric(names(T0_counts)[which
                                       (T0_counts == max(T0_counts))]))
  } else {
    return (1)
  }
}

sars2.nkcds <- orderCells(sars2.nkcds, root_state = GM_state(sars2.nkcds))
plot_cell_trajectory(sars2.nkcds, color_by = "Pseudotime")


###nk PSEUDOTIME ####

setwd("C:/Users/Dylan2/Google Drive/UNI WORK/2020 Internship/SARS2/revised 211120/pseudo_plots")

path <- "C:/Users/Dylan2/Google Drive/UNI WORK/2020 Internship/SARS2/revised 211120/pseudo_plots/other_celltypes/"

#plot pseudotime
ptime.tjct(sars2.nkcds, max=18, rev=F)
#plot ifnmod
ifnmod.tjct(sars2.nkcds)
#plot treatment
treat.tjct(sars2.nkcds)




## CD4 CELLS -----------------------------------------------------------------------------------------------------------------


#subcluster

sars2.cd4 <-  subset(sars2.all, CellTypes == "CD4+ T cell")

DefaultAssay(sars2.cd4)<-"RNA"

cd4.list <- SplitObject(sars2.cd4, split.by = "orig.ident")

for (i in 1:length(cd4.list)) {
  cd4.list[[i]] <- SCTransform(cd4.list[[i]], verbose = FALSE)
}
options(future.globals.maxSize =1 * 4000 * 1024^2)
pancreas.features <- SelectIntegrationFeatures(object.list = cd4.list, nfeatures = 2000)


cd4.list <- PrepSCTIntegration(object.list = cd4.list, anchor.features = pancreas.features, 
                               verbose = FALSE)
pancreas.anchors <- FindIntegrationAnchors(object.list = cd4.list, normalization.method = "SCT", 
                                           anchor.features = pancreas.features, verbose = FALSE)
pancreas.integrated <- IntegrateData(anchorset = pancreas.anchors, normalization.method = "SCT", verbose = FALSE)

sars2.cd4 <- RunPCA(pancreas.integrated, npcs = 10, verbose = FALSE)
sars2.cd4 <- RunUMAP(sars2.cd4, dims = 1:10, reduction.name = "umapsub", reduction.key = "UMAPSUB")


sars2.cd4 <- FindNeighbors(sars2.cd4, reduction = "pca", dims = 1:10)
sars2.cd4 <- FindClusters(sars2.cd4, resolution = 0.3)
sars2.cd4$seurat_clusters_int <- sars2.cd4$seurat_clusters

rm(cd4.list)
##UMAPS
clust_umap <- DimPlot(sars2.cd4,label = T , reduction = "umapsub", group.by="seurat_clusters_int") + NoAxes() + coord_fixed() 
clust_umap
ggsave(paste0(path, "cd4_umap_cluster.png"), plot= clust_umap)

treat_umap<- DimPlot(sars2.cd4,label = T , reduction = "umapsub", group.by = "HTO_maxID") + NoAxes() + coord_fixed() 
treat_umap
ggsave(paste0(path, "cd4_umap_treat.png"), plot= treat_umap)


##add viral reads and ifn mod


#add IFN signalling module score
sars2.cd4 <- AddModuleScore(object = sars2.cd4, features = ifn_sig, ctrl = 100, name = 'ifn_sig', assay="RNA")


#### CD4 CDS ####


## IMPORT TO MONOCLE #

data <- as(as.matrix(sars2.cd4@assays$RNA@counts), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = sars2.cd4@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)

#Construct cd4cle cds
sars2.cd4cds <- newCellDataSet(data,
                               phenoData = pd,
                               featureData = fd,
                               lowerDetectionLimit = 0.5,
                               expressionFamily = negbinomial.size())


#estimate size factors and dispersion
sars2.cd4cds <- estimateSizeFactors(sars2.cd4cds)
sars2.cd4cds <- estimateDispersions(sars2.cd4cds)
#remove genes not expressed in at least 10 cells
sars2.cd4cds <- detectGenes(sars2.cd4cds, min_expr = 0.1)
print(head(fData(sars2.cd4cds)))
expressed_genes <- row.names(subset(fData(sars2.cd4cds),
                                    num_cells_expressed >= 10))




#determine ordering genes - diff expressed between infections

#Idents(sars2.cd4) <- "HTO_maxID"
#ordering_genes <- FindAllMarkers(sars2.cd4)

Idents(sars2.cd4) <- "seurat_clusters_int"
ordering_genes <- FindAllMarkers(sars2.cd4)

ordering_genes <- ordering_genes[ordering_genes$p_val_adj < 0.05,]
ordering_genes$log2 <- log2(exp(ordering_genes$avg_logFC))
ordering_genes <- ordering_genes[abs(ordering_genes$log2) > 2,]

ordering_genes <- as.vector(ordering_genes$gene)

#set ordering genes
sars2.cd4cds <- setOrderingFilter(sars2.cd4cds, ordering_genes)
#reduce dimensions
sars2.cd4cds <- reduceDimension(sars2.cd4cds, max_components = 2,
                                method = 'DDRTree')
#order cells along trajectory
sars2.cd4cds <- orderCells(sars2.cd4cds )

#check root state
plot_cell_trajectory(sars2.cd4cds, color_by = "HTO_maxID")

plot_cell_trajectory(sars2.cd4cds, color_by = "State") +
  facet_wrap(~State, nrow = 1)

#cds <- sars2.cd4cds

##determine which state has most cells from mock
GM_state <- function(cds){
  if (length(unique(pData(cds)$State)) > 1){
    T0_counts <- table(pData(cds)$State, pData(cds)$HTO_maxID)[,"Mock"]
    return(as.numeric(names(T0_counts)[which
                                       (T0_counts == max(T0_counts))]))
  } else {
    return (1)
  }
}

sars2.cd4cds <- orderCells(sars2.cd4cds, root_state = GM_state(sars2.cd4cds))
plot_cell_trajectory(sars2.cd4cds, color_by = "Pseudotime")


###CD4  PSEUDOTIME ####

#plot pseudotime
ptime.tjct(sars2.cd4cds, max=18, rev=F)
#plot ifnmod
ifnmod.tjct(sars2.cd4cds)
#plot treatment
treat.tjct(sars2.cd4cds)





## CD8 CELLS -----------------------------------------------------------------------------------------------------------------

#subcluster

sars2.cd8 <-  subset(sars2.all, CellTypes == "CD8+ T cell")

DefaultAssay(sars2.cd8)<-"RNA"

cd8.list <- SplitObject(sars2.cd8, split.by = "orig.ident")

for (i in 1:length(cd8.list)) {
  cd8.list[[i]] <- SCTransform(cd8.list[[i]], verbose = FALSE)
}
options(future.globals.maxSize =1 * 4000 * 1024^2)
pancreas.features <- SelectIntegrationFeatures(object.list = cd8.list, nfeatures = 2000)


cd8.list <- PrepSCTIntegration(object.list = cd8.list, anchor.features = pancreas.features, 
                               verbose = FALSE)
pancreas.anchors <- FindIntegrationAnchors(object.list = cd8.list, normalization.method = "SCT", 
                                           anchor.features = pancreas.features, verbose = FALSE)
pancreas.integrated <- IntegrateData(anchorset = pancreas.anchors, normalization.method = "SCT", verbose = FALSE)

sars2.cd8 <- RunPCA(pancreas.integrated, npcs = 10, verbose = FALSE)
sars2.cd8 <- RunUMAP(sars2.cd8, dims = 1:10, reduction.name = "umapsub", reduction.key = "UMAPSUB")


sars2.cd8 <- FindNeighbors(sars2.cd8, reduction = "pca", dims = 1:10)
sars2.cd8 <- FindClusters(sars2.cd8, resolution = 0.3)
sars2.cd8$seurat_clusters_int <- sars2.cd8$seurat_clusters


##UMAPS
clust_umap <- DimPlot(sars2.cd8,label = T , reduction = "umapsub", group.by="seurat_clusters_int") + NoAxes() + coord_fixed() 
clust_umap
ggsave(paste0(path, "cd8_umap_cluster.png"), plot= clust_umap)

treat_umap<- DimPlot(sars2.cd8,label = T , reduction = "umapsub", group.by = "HTO_maxID") + NoAxes() + coord_fixed() 
treat_umap
ggsave(paste0(path, "cd8_umap_treat.png"), plot= treat_umap)


##add viral reads and ifn mod


#add IFN signalling module score
sars2.cd8 <- AddModuleScore(object = sars2.cd8, features = ifn_sig, ctrl = 100, name = 'ifn_sig', assay="RNA")

#### cd8 CDS ####


## IMPORT TO MONOCLE #

data <- as(as.matrix(sars2.cd8@assays$RNA@counts), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = sars2.cd8@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)

#Construct cd8cle cds
sars2.cd8cds <- newCellDataSet(data,
                               phenoData = pd,
                               featureData = fd,
                               lowerDetectionLimit = 0.5,
                               expressionFamily = negbinomial.size())


#estimate size factors and dispersion
sars2.cd8cds <- estimateSizeFactors(sars2.cd8cds)
sars2.cd8cds <- estimateDispersions(sars2.cd8cds)
#remove genes not expressed in at least 10 cells
sars2.cd8cds <- detectGenes(sars2.cd8cds, min_expr = 0.1)
print(head(fData(sars2.cd8cds)))
expressed_genes <- row.names(subset(fData(sars2.cd8cds),
                                    num_cells_expressed >= 10))



Idents(sars2.cd8) <- "seurat_clusters_int"
ordering_genes <- FindAllMarkers(sars2.cd8)

ordering_genes <- ordering_genes[ordering_genes$p_val_adj < 0.05,]
ordering_genes$log2 <- log2(exp(ordering_genes$avg_logFC))
ordering_genes <- ordering_genes[abs(ordering_genes$log2) > 2,]

ordering_genes <- as.vector(ordering_genes$gene)

#set ordering genes
sars2.cd8cds <- setOrderingFilter(sars2.cd8cds, ordering_genes)
#reduce dimensions
sars2.cd8cds <- reduceDimension(sars2.cd8cds, max_components = 2,
                                method = 'DDRTree')
#order cells along trajectory
sars2.cd8cds <- orderCells(sars2.cd8cds )


#check root state
plot_cell_trajectory(sars2.cd8cds, color_by = "HTO_maxID")

plot_cell_trajectory(sars2.cd8cds, color_by = "State") +
  facet_wrap(~State, nrow = 1)

#cds <- sars2.cd8cds

##determine which state has most cells from mock
GM_state <- function(cds){
  if (length(unique(pData(cds)$State)) > 1){
    T0_counts <- table(pData(cds)$State, pData(cds)$HTO_maxID)[,"Mock"]
    return(as.numeric(names(T0_counts)[which
                                       (T0_counts == max(T0_counts))]))
  } else {
    return (1)
  }
}

sars2.cd8cds <- orderCells(sars2.cd8cds, root_state = GM_state(sars2.cd8cds))
plot_cell_trajectory(sars2.cd8cds, color_by = "Pseudotime")
plot_cell_trajectory(sars2.cd8cds, color_by = "predicted.id")

###CD8  PSEUDOTIME ####

#plot pseudotime
ptime.tjct(sars2.cd8cds, max=18, rev=F)
#plot ifnmod
ifnmod.tjct(sars2.cd8cds)
#plot treatment
treat.tjct(sars2.cd8cds)


