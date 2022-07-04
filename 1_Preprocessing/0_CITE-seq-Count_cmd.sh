#!/bin/bash

## Run CITE-seq count pipeline to assign barcodes of cells to original treatments based on HTOs associated with cells within 10X Chromium GEMs

## Donor1

CITE-seq-Count \
-R1 /.../CG_JK_11_HTO_S7_L005_R1_001.fastq.gz \
-R2 /.../CG_JK_11_HTO_S7_L005_R2_001.fastq.gz \
-t ./hto_mnc.csv \
-cbf 1 \
-cbl 16 \
-umif 17 \
-umil 28 \
-cells 12369 \
--bc_collapsing_dist 3 \
--umi_collapsing_dist 3 \
--max-error 3 \
--unmapped-tags unmapped.csv \
-T 10 \
-o CG_JK_11_10x_02 # orig id: CG_JK_11_10x_02, last number was optim. run id

## Donor2

CITE-seq-Count \
-R1 /.../CG_JK_12_HTO_S10_L006_R1_001.fastq.gz \
-R2 /.../CG_JK_12_HTO_S10_L006_R2_001.fastq.gz \
-t ./hto_mnc.csv \
-cbf 1 \
-cbl 16 \
-umif 17 \
-umil 28 \
-cells 15822 \
--bc_collapsing_dist 3 \
--umi_collapsing_dist 3 \
--sliding-window \
--max-error 3 \
--unmapped-tags unmapped.csv \
-T 10 \
-o CG_JK_12_10x # orig id: G_JK_12_10x_11, last number was optim. run id