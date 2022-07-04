library(Seurat)
library(tidyverse)

#--------------------------------------------------------------------------

## read in full, annotated PBMC Seurat object (with viral read allocation)
pbmc <- readRDS("path/to/annotated/Seurat_object_with_viral_read_counts.rds")
## read in subclustered Monocyte Seurat object
mono <- readRDS("subclustered_monocyte_object.rds")

#--------------------------------------------------------------------------


## creat dataframe to store annotation data for downstream plotting
annot_df <- pbmc@meta.data %>%
  subset(CellTypes = "Monocyte") %>%
  select(barcode, viral, virus, viral_count)

mono$barcode <- rownames(mono@meta.data)
mono_meta <- mono@meta.data

join <- plyr::join(mono_meta, annot_df)


mono$viral <- join$viral
mono$virus <- join$virus
mono$viral_count <- join$viral_count



## add IFN module score (Reactome IFN signalling pathway)

ifn_sig <- read.table("C:/Users/postmusd/OneDrive - Charité - Universitätsmedizin Berlin/Dokumente/Push/SARS2_PBMC/Reactome Interferon Signalling [R-HSA-913531].tsv",
                      sep = "\t", header = T)

head(ifn_sig)
ifn_sig$gene <- gsub("UniProt:.* (.*)$", "\\1", ifn_sig$MoleculeName)
ifn_sig <- ifn_sig$gene
ifn_sig <- ifn_sig[ifn_sig %in% rownames(mono) ]
ifn_sig <- list(ifn_sig)


## add module score
mono <- AddModuleScore(
  mono,
  ifn_sig,
  pool = NULL,
  nbin = 24,
  ctrl = 100,
  k = FALSE,
  assay = "RNA",
  name = "ifn_sig",
  seed = 1)

head(mono)

### VIOLINS

mono_df <- mono@meta.data %>%
  subset(HTO_maxID %in% c("SARS1", "SARS2"))

## visualise violins by treatment

vln.infect <- ggplot(mono_df, aes(HTO_maxID, ifn_sig1, fill=HTO_maxID)) +
  geom_violin() +
  theme_classic() 

vln.infect

ggsave("./ifnmod_vln/ifnmod_vln_byvirus.pdf", 
       plot =vln.infect)
ggsave("./ifnmod_vln_byvirus.eps",
       plot = vln.infect)

## with median

vln.infect <- ggplot(mono_df, aes(HTO_maxID, ifn_sig1, fill=HTO_maxID)) +
  geom_violin() +
  theme_classic() +
  stat_summary(fun = "mean",
               geom = "crossbar", 
               width = 0.2,
               colour = "black")

vln.infect

ggsave("./ifnmod_vln_byvirus_median.pdf", 
       plot =vln.infect)
ggsave("./ifnmod_vln_byvirus_median.eps",
       plot = vln.infect)

###############################################

## by viral positive or negative cells

vln.infect_viralpos <- ggplot(mono_df, aes(HTO_maxID, ifn_sig1, fill=viral)) +
  geom_violin(position = "dodge") +
  theme_classic() 

vln.infect_viralpos

ggsave("./ifnmod_vln_byviral.positive.pdf", 
       plot =vln.infect_viralpos)
ggsave("./ifnmod_vln_byviral.positive.eps",
       plot = vln.infect_viralpos)

## with median
vln.infect_viralpos <- ggplot(mono_df, aes(HTO_maxID, ifn_sig1, fill=viral)) +
  geom_violin(position = "dodge") +
  theme_classic() +
  stat_summary(fun = "median",
               geom = "crossbar", 
               width = 0.2,
               colour = "black",
               position = position_dodge(width = 0.9))
vln.infect_viralpos

ggsave("./ifnmod_vln_byviral.positive_median.pdf", 
       plot =vln.infect_viralpos)
ggsave("./ifnmod_vln_byviral.positive_median.eps",
       plot =vln.infect_viralpos)

#### STATISTICAL TEST ####

# test mono ifn score
# SARS2


# viral + vs viral -

mono_ifn <- mono@meta.data
head(mono_ifn)

wtest <- wilcox.test(mono_ifn[mono_ifn$viral == "pos" & mono_ifn$HTO_maxID == "SARS2",]$ifn_sig1, mono_ifn[mono_ifn$viral == "neg" & mono_ifn$HTO_maxID == "SARS2",]$ifn_sig1)
wtest
mono_ifn[mono_ifn$viral == "pos" & mono_ifn$HTO_maxID == "SARS2",]$ifn_sig1

# SARS1

wtest <- wilcox.test(mono_ifn[mono_ifn$viral == "pos" & mono_ifn$HTO_maxID == "SARS1",]$ifn_sig1, mono_ifn[mono_ifn$viral == "neg" & mono_ifn$HTO_maxID == "SARS1",]$ifn_sig1)
wtest

# SARS2 vs SARS1
wtest <- wilcox.test(mono_ifn[mono_ifn$HTO_maxID == "SARS2",]$ifn_sig1, mono_ifn[mono_ifn$HTO_maxID == "SARS1",]$ifn_sig1)
wtest

