#!/bin/bash

#------------------------------------------------------------------------------------------------------------------------
## script to:
# 1. sort reads from cell ranger bams (alignment to SARS1 and 2 combo genomes) into mock, sars1 and sars2 for each donor based on a list of
# 	barcodes obtained from Seurat object and HTO classification and compile report of representation of reads mapped to each genome
# 2. filter subsetted bam files to include only reads mapping to reads from expected genomes, then get barcodes corresponding to reads as well as
#	UMIs for reads in cells for downstream analysis and identification of virus + cells
#------------------------------------------------------------------------------------------------------------------------

#DONOR 11 #######################################################################################

d=CG_JK_11_10x_AY_NC_viralreads.sorted.bam

##filter reads based on barcodes

for bc in `ls -1 | grep bc_d11.*.txt`;

do

# Save the header lines
samtools view -H $d > SAM_header

n=${d/_10x_AY_NC_viralreads.sorted.bam/} ; ###native bash function to rename variable
v=$(echo $bc | sed 's/bc_d[0-9]\{2\}_\(.*\).txt/\1/') ;

# Filter alignments using filter.txt. Use LC_ALL=C to set C locale instead of UTF-8
samtools view $d | LC_ALL=C grep -F -f $bc > filtered_SAM_body ;


# Combine header and body
cat SAM_header filtered_SAM_body > filtered.sam;

samtools view -b filtered.sam > "$n"_"$v"_bc_filter.bam ;
samtools index "$n"_"$v"_bc_filter.bam ; #index so that sorting works down the line

# Get list of barcodes

##SARS1 reads
samtools view "$n"_"$v"_bc_filter.bam AY310120.1 -b > samtools sort -o "$n"_"$v"_bc_filter_SARS1.bam ; 

## get barcode, umi and 10x umi tag + filter out secondary alignments (flag 256) and mapq of at least 30 
samtools view -F 256 -q 30 "$n"_"$v"_bc_filter_SARS1.bam | awk '{for(j=1;j<=NF;j++){if($j~/^CB:Z:/){print $j}}}' | sed 's/CB:Z://' > viral_allocation/test_bc;
samtools view -F 256 -q 30 "$n"_"$v"_bc_filter_SARS1.bam | awk '{for(j=1;j<=NF;j++){if($j~/^UB:Z:/){print $j}}}' | sed 's/UB:Z://' > viral_allocation/test_umi;
samtools view -F 256 -q 30 "$n"_"$v"_bc_filter_SARS1.bam | awk '{for(j=1;j<=NF;j++){if($j~/^xf:i:/){print $j}}}' | sed 's/UB:Z://' > viral_allocation/test_tag;
#filter out xf:i:19 - UMIs which mapped to features that most others did not
paste -d "\t" viral_allocation/test_bc viral_allocation/test_umi viral_allocation/test_tag | grep -v "xf:i:19" >  viral_allocation/bc_viral_SARS1_"$n"_"$v".txt;


##SARS2 reads
samtools view "$n"_"$v"_bc_filter.bam NC_045512.2 -b > samtools sort -o "$n"_"$v"_bc_filter_SARS2.bam ;

## get barcode, umi and 10x umi tag + filter out secondary alignments (flag 256)
samtools view -F 256 -q 30 "$n"_"$v"_bc_filter_SARS2.bam | awk '{for(j=1;j<=NF;j++){if($j~/^CB:Z:/){print $j}}}' | sed 's/CB:Z://' > viral_allocation/test_bc;
samtools view -F 256 -q 30 "$n"_"$v"_bc_filter_SARS2.bam | awk '{for(j=1;j<=NF;j++){if($j~/^UB:Z:/){print $j}}}' | sed 's/UB:Z://' > viral_allocation/test_umi;
samtools view -F 256 -q 30 "$n"_"$v"_bc_filter_SARS2.bam | awk '{for(j=1;j<=NF;j++){if($j~/^xf:i:/){print $j}}}' | sed 's/UB:Z://' > viral_allocation/test_tag;
#filter out xf:i:19 - UMIs which mapped to features that most others did not
paste -d "\t" viral_allocation/test_bc viral_allocation/test_umi viral_allocation/test_tag | grep -v "xf:i:19" >  viral_allocation/bc_viral_SARS2_"$n"_"$v".txt;

done



#DONOR 12 #######################################################################################

d=CG_JK_12_10x_AY_NC_viralreads.sorted.bam

##filter reads based on barcodes

for bc in `ls -1 | grep bc_d12.*.txt`;

do

# Save the header lines
samtools view -H $d > SAM_header

n=${d/_10x_AY_NC_viralreads.sorted.bam/} ; ###native bash function to rename variable
v=$(echo $bc | sed 's/bc_d[0-9]\{2\}_\(.*\).txt/\1/') ;

# Filter alignments using filter.txt. Use LC_ALL=C to set C locale instead of UTF-8
samtools view $d | LC_ALL=C grep -F -f $bc > filtered_SAM_body ;

# Combine header and body
cat SAM_header filtered_SAM_body > filtered.sam;

samtools view -b filtered.sam > "$n"_"$v"_bc_filter.bam ;
samtools index "$n"_"$v"_bc_filter.bam ; #index so that sorting works down the line

# Get list of barcodes

##SARS1 reads
samtools view "$n"_"$v"_bc_filter.bam AY310120.1 -b > samtools sort -o "$n"_"$v"_bc_filter_SARS1.bam ; 

## get barcode, umi and 10x umi tag + filter out secondary alignments (flag 256)
samtools view -F 256 -q 30 "$n"_"$v"_bc_filter_SARS1.bam | awk '{for(j=1;j<=NF;j++){if($j~/^CB:Z:/){print $j}}}' | sed 's/CB:Z://' > viral_allocation/test_bc;
samtools view -F 256 -q 30 "$n"_"$v"_bc_filter_SARS1.bam | awk '{for(j=1;j<=NF;j++){if($j~/^UB:Z:/){print $j}}}' | sed 's/UB:Z://' > viral_allocation/test_umi;
samtools view -F 256 -q 30 "$n"_"$v"_bc_filter_SARS1.bam | awk '{for(j=1;j<=NF;j++){if($j~/^xf:i:/){print $j}}}' | sed 's/UB:Z://' > viral_allocation/test_tag;
#filter out xf:i:19 - UMIs which mapped to features that most others did not
paste -d "\t" viral_allocation/test_bc viral_allocation/test_umi viral_allocation/test_tag | grep -v "xf:i:19" >  viral_allocation/bc_viral_SARS1_"$n"_"$v".txt;


##SARS2 reads
samtools view "$n"_"$v"_bc_filter.bam NC_045512.2 -b > samtools sort -o "$n"_"$v"_bc_filter_SARS2.bam ;

## get barcode, umi and 10x umi tag + filter out secondary alignments (flag 256)
samtools view -F 256 -q 30 "$n"_"$v"_bc_filter_SARS2.bam | awk '{for(j=1;j<=NF;j++){if($j~/^CB:Z:/){print $j}}}' | sed 's/CB:Z://' > viral_allocation/test_bc;
samtools view -F 256 -q 30 "$n"_"$v"_bc_filter_SARS2.bam | awk '{for(j=1;j<=NF;j++){if($j~/^UB:Z:/){print $j}}}' | sed 's/UB:Z://' > viral_allocation/test_umi;
samtools view -F 256 -q 30 "$n"_"$v"_bc_filter_SARS2.bam | awk '{for(j=1;j<=NF;j++){if($j~/^xf:i:/){print $j}}}' | sed 's/UB:Z://' > viral_allocation/test_tag;
#filter out xf:i:19 - UMIs which mapped to features that most others did not
paste -d "\t" viral_allocation/test_bc viral_allocation/test_umi viral_allocation/test_tag | grep -v "xf:i:19" >  viral_allocation/bc_viral_SARS2_"$n"_"$v".txt;

done

#=======================================================================================================================================









