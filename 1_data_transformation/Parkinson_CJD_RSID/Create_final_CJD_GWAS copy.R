library(dplyr)
library(tidyverse)

# Merge SNPs with GWAS summary statistics and creatae final dataset 

setwd("~/Desktop/Cambridge/Dissertation/Get Parkinson and CJD SNPs")

# Open all datasets

final_CJD_SNPs_chr1 <- read.delim("~/Desktop/Cambridge/Dissertation/Get Parkinson and CJD SNPs/completely_final_CJD/final_CJD_SNPs_chr1.tsv")
final_CJD_SNPs_chr2 <- read.delim("~/Desktop/Cambridge/Dissertation/Get Parkinson and CJD SNPs/completely_final_CJD/final_CJD_SNPs_chr2.tsv")
final_CJD_SNPs_chr3 <- read.delim("~/Desktop/Cambridge/Dissertation/Get Parkinson and CJD SNPs/completely_final_CJD/final_CJD_SNPs_chr3.tsv")
final_CJD_SNPs_chr4 <- read.delim("~/Desktop/Cambridge/Dissertation/Get Parkinson and CJD SNPs/completely_final_CJD/final_CJD_SNPs_chr4.tsv")
final_CJD_SNPs_chr5 <- read.delim("~/Desktop/Cambridge/Dissertation/Get Parkinson and CJD SNPs/completely_final_CJD/final_CJD_SNPs_chr5.tsv")
final_CJD_SNPs_chr6 <- read.delim("~/Desktop/Cambridge/Dissertation/Get Parkinson and CJD SNPs/completely_final_CJD/final_CJD_SNPs_chr6.tsv")
final_CJD_SNPs_chr7 <- read.delim("~/Desktop/Cambridge/Dissertation/Get Parkinson and CJD SNPs/completely_final_CJD/final_CJD_SNPs_chr7.tsv")
final_CJD_SNPs_chr8 <- read.delim("~/Desktop/Cambridge/Dissertation/Get Parkinson and CJD SNPs/completely_final_CJD/final_CJD_SNPs_chr8.tsv")
final_CJD_SNPs_chr9 <- read.delim("~/Desktop/Cambridge/Dissertation/Get Parkinson and CJD SNPs/completely_final_CJD/final_CJD_SNPs_chr9.tsv")
final_CJD_SNPs_chr10 <- read.delim("~/Desktop/Cambridge/Dissertation/Get Parkinson and CJD SNPs/completely_final_CJD/final_CJD_SNPs_chr10.tsv")
final_CJD_SNPs_chr11 <- read.delim("~/Desktop/Cambridge/Dissertation/Get Parkinson and CJD SNPs/completely_final_CJD/final_CJD_SNPs_chr11.tsv")
final_CJD_SNPs_chr12 <- read.delim("~/Desktop/Cambridge/Dissertation/Get Parkinson and CJD SNPs/completely_final_CJD/final_CJD_SNPs_chr12.tsv")
final_CJD_SNPs_chr13 <- read.delim("~/Desktop/Cambridge/Dissertation/Get Parkinson and CJD SNPs/completely_final_CJD/final_CJD_SNPs_chr13.tsv")
final_CJD_SNPs_chr14 <- read.delim("~/Desktop/Cambridge/Dissertation/Get Parkinson and CJD SNPs/completely_final_CJD/final_CJD_SNPs_chr14.tsv")
final_CJD_SNPs_chr15 <- read.delim("~/Desktop/Cambridge/Dissertation/Get Parkinson and CJD SNPs/completely_final_CJD/final_CJD_SNPs_chr15.tsv")
final_CJD_SNPs_chr16 <- read.delim("~/Desktop/Cambridge/Dissertation/Get Parkinson and CJD SNPs/completely_final_CJD/final_CJD_SNPs_chr16.tsv")
final_CJD_SNPs_chr17 <- read.delim("~/Desktop/Cambridge/Dissertation/Get Parkinson and CJD SNPs/completely_final_CJD/final_CJD_SNPs_chr17.tsv")
final_CJD_SNPs_chr18 <- read.delim("~/Desktop/Cambridge/Dissertation/Get Parkinson and CJD SNPs/completely_final_CJD/final_CJD_SNPs_chr18.tsv")
final_CJD_SNPs_chr19 <- read.delim("~/Desktop/Cambridge/Dissertation/Get Parkinson and CJD SNPs/completely_final_CJD/final_CJD_SNPs_chr19.tsv")
final_CJD_SNPs_chr20 <- read.delim("~/Desktop/Cambridge/Dissertation/Get Parkinson and CJD SNPs/completely_final_CJD/final_CJD_SNPs_chr20.tsv")
final_CJD_SNPs_chr21 <- read.delim("~/Desktop/Cambridge/Dissertation/Get Parkinson and CJD SNPs/completely_final_CJD/final_CJD_SNPs_chr21.tsv")
final_CJD_SNPs_chr22 <- read.delim("~/Desktop/Cambridge/Dissertation/Get Parkinson and CJD SNPs/completely_final_CJD/final_CJD_SNPs_chr22.tsv")

# Merge them

final_SNP_list <- rbind(final_CJD_SNPs_chr1,final_CJD_SNPs_chr2,final_CJD_SNPs_chr3,final_CJD_SNPs_chr4,final_CJD_SNPs_chr5,
                        final_CJD_SNPs_chr6,final_CJD_SNPs_chr7,final_CJD_SNPs_chr8,final_CJD_SNPs_chr9,final_CJD_SNPs_chr10,
                        final_CJD_SNPs_chr11,final_CJD_SNPs_chr12,final_CJD_SNPs_chr13,final_CJD_SNPs_chr14,final_CJD_SNPs_chr15,
                        final_CJD_SNPs_chr16,final_CJD_SNPs_chr17,final_CJD_SNPs_chr18,final_CJD_SNPs_chr19,final_CJD_SNPs_chr20,
                        final_CJD_SNPs_chr21,final_CJD_SNPs_chr22)

# Remove bp2 line

final_SNP_list <- final_SNP_list[,-3]

# Read CJD GWAS data

CJD <- read.delim("~/Desktop/Cambridge/Dissertation/Get Parkinson and CJD SNPs/CJD.tsv")

# Merge them by chromosome, base pair position, reference and alternative allele

CJD_GWAS_with_SNPs <- merge.data.frame(CJD, final_SNP_list, by = c("chromosome", "base_pair_location", "effect_allele", "other_allele"))

# Clean up data and save

CJD_GWAS_with_SNPs <- CJD_GWAS_with_SNPs[,-c(12,13,15,16)]
CJD_GWAS_with_SNPs <- CJD_GWAS_with_SNPs[,c(1,2,12,3,4,5,6,7,8,9,10,11)]

write.table(CJD_GWAS_with_SNPs, "final_CJD_with_SNPs.tsv", sep = "\t", col.names = TRUE)
