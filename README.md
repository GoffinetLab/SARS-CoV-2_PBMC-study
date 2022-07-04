# SARS-CoV-2_PBMC-study

## ABSTRACT

Cell-intrinsic responses mounted in vivo in PBMCs during mild and severe COVID-19 differ quantitatively and qualitatively. Whether they are triggered by signals emitted by productively infected cells of the respiratory tract or are, at least partially, resulting from physical interaction with virus particles, remains unclear. Here, we analyzed susceptibility and expression profiles of PBMCs from healthy donors upon ex vivo exposure to SARS-CoV and SARS-CoV-2. In line with the absence of detectable ACE2 receptor expression, human PBMCs were refractory to productive infection. Bulk and single cell RNA-sequencing revealed JAK/STAT-dependent induction of interferon-stimulated genes (ISGs), but not pro-inflammatory cytokines. This SARS-CoV-2-specific response was most pronounced in monocytes. SARS-CoV-2-RNA-positive monocytes displayed a lower ISG signature as compared to bystander cells of the identical culture. This suggests a preferential invasion of cells with a low ISG base-line profile or delivery of a SARS-CoV-2-specific sensing antagonist upon efficient particle internalization. Together, non-productive physical interaction of PBMCs with SARS-CoV-2-, and to a much lesser extent, SARS-CoV particles stimulates JAK/STAT-dependent, monocyte-accentuated innate immune responses that resemble those detected in vivo in patients with mild COVID-19.


## OVERVIEW

This is the Github repository for the manuscript "Non-productive exposure of PBMCs to SARS-CoV-2 induces cell-intrinsic innate immune responses" by Kazmierski et al., 2022. It contains all the code utilised in the analysis of the single cell RNAseq data. The analysis is divided into "1_Preprocessing" and "2_Analysis". "1_Preprocessing" contains all the code code used for the initial preprocessing of the 10X V3.1 3' single cell RNAseq data, as well as the HTO classification CITE-seq data, the latter of which was used purely to assign cells to an original treatment after the cells for each donor were pooled for all three treatments ("Mock exposed", "SARS-CoV exposed" and SARS-CoV-2 exposed") for a multiplexed Chromium controller run. Each cell barcode was allocated to the different treatments based on the enrichment of the HTO oligo assigned to each treatment (Mock: GGTTGCCAGATGTCA; SARS-CoV: CAGTAGTCACGGTCA; SARS-CoV-2: CTCCTCTGCAATTAC) associated with that barcode. Included also are the results of the CITE-seq count pipeline in the form of two csv files (Donor1: CG_JK_11_10x_HTO_counts.csv; Donor2: CG_JK_11_10x_HTO_counts.csv) that can be used to assign treatments to cells based on the HTO IDs. "2_Analysis" contains code snippets that were used for further processing of the data, including subclustering and Pseudotime analysis.


## RAW AND PROCESSED DATA

The raw and processed data can be downloaded from GEO.
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE197665
