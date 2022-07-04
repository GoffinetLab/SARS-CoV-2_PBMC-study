#!/bin/sh

#------------------------------------------------------------------------------------------------------------------------
## script to:
# align FASTQs to combined SARS-CoV / SARS-CoV-2 reference. To create the reference, use cellranger mkref command, using the files:
# 1. SARS1_2_AY_NC.fa
# 2. SARS1_2_AY_NC.gtf
# referred to here as SARS1_2_ref
#------------------------------------------------------------------------------------------------------------------------

## DONOR11

cellranger count \
--id=CG_JK_11_10x_AY_NC \
--transcriptome=/SARS1_2_ref \
--fastqs=/CG_JK_11_fastqs \
--sample=CG_JK_11_10x \
--chemistry=SC3Pv3

## DONOR12

cellranger count \
--id=CG_JK_12_10x_AY_NC \
--transcriptome=/SARS1_2_ref \
--fastqs=/CG_JK_12_fastqs \
--sample=CG_JK_12_10x \
--chemistry=SC3Pv3

