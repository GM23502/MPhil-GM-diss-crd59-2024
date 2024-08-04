####################### Two sample MR for eQTL and Parkinson ################

library(TwoSampleMR)
library(ggplot2)
library(tidyverse)
library(genetics)
library(gt)
library(patchwork)
library(grid)
library(annotate)
library(ggeasy)

# Set working directory

setwd("~/Desktop/MPhil-GM-diss-crd59-2024-main/CR1")

# Read exposure dataset

pQTL_exposure <- 
  read_exposure_data("pQTL_exposure.csv",
                     sep = ",",
                     clump = TRUE,
                     snp_col = "SNP",
                     beta_col = "B",
                     se_col = "SE",
                     eaf_col = "EAF",
                     effect_allele_col = "A1",
                     other_allele_col = "A2",
                     pval_col = "P",
                     ncase_col = "N",
                     gene_col = "gene",
                     chr_col = "CHR",
                     pos_col = "BP", 
                     min_pval = 0)

# Read outcome dataset

Alzheimer_outcome <-
  read_outcome_data("Alzheimer_outcome.csv",
                    sep = ",",
                    snp_col = "SNP",
                    beta_col = "B",
                    se_col = "SE",
                    eaf_col = "EAF",
                    effect_allele_col = "A1",
                    other_allele_col = "A2",
                    pval_col = "P",
                    ncase_col = "N",
                    chr_col = "CHR",
                    pos_col = "BP")

# Harmonize data

H_data <- harmonise_data(exposure_dat = pQTL_exposure, outcome_dat = Alzheimer_outcome)

# perform first-stage F-statistic

H_data$Fstat <- (2*H_data$eaf.exposure*(1 - H_data$eaf.exposure) * H_data$beta.exposure^2)*(31684 - 2)/(1 - (2*H_data$eaf.exposure*
                                                                                                               (1 - H_data$eaf.exposure) * H_data$beta.exposure^2))

# Save IVS

write.csv(H_data, "./IVs_pQTL_Alzheimer.csv")

#how much the MR-Egger intercept is non-zero
pleiotropy <- mr_pleiotropy_test(H_data)

write.csv(pleiotropy, "./CR1_pQTL_Alzheimer_Pleiotropy.csv")

#Q statistics for heterogeneity
heterogeneity <- mr_heterogeneity(H_data, method_list=c("mr_egger_regression", "mr_ivw"))

write.csv(heterogeneity, "./CR1_pQTL_Alzheimer_Heterogeneity.csv")

# Significant heterogeneity, will perform MR analyses
