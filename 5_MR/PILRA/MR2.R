#################### PILRA MR ANALYSES ########################


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

setwd("~/Desktop/MPhil-GM-diss-crd59-2024-main/PILRA")

##################### START MR PROCEDURE ######################################

# Read eQTL exposure

eQTL_exposure <- read_exposure_data("eQTL_exposure.csv",
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

# Read pQTL exposure

pQTL_exposure <- read_exposure_data("pQTL_exposure.csv",
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
                               chr_col = "CHR",
                               pos_col = "BP",
                               min_pval = 0)

# Read outcome

Alzheimer_outcome <- read_outcome_data("Alzheimer_outcome.csv",
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


############################ eQTL and AD #####################################

H_data <- harmonise_data(exposure_dat = eQTL_exposure, outcome_dat = Alzheimer_outcome)

#how much the MR-Egger intercept is non-zero
pleiotropy <- mr_pleiotropy_test(H_data)

write.csv(pleiotropy, "PIRLA_Pleiotropy_eQTL.csv")

#Q statistics for heterogeneity
heterogeneity <- mr_heterogeneity(H_data, method_list=c("mr_egger_regression", "mr_ivw"))

write.csv(heterogeneity, "heterogeneity_eQTL_PIRLA.csv")

H_data$Fstat <- (2*H_data$eaf.exposure*(1 - H_data$eaf.exposure) * H_data$beta.exposure^2)*(31684 - 2)/(1 - (2*H_data$eaf.exposure*
                                                                                                               (1 - H_data$eaf.exposure) * H_data$beta.exposure^2))

write.csv(H_data, "./IVs_PILRA_eQTLs_Alzheimer.csv")


# Due to significant heterogeneity we did not performt the MR analyses using the
# eQTL-derived IVs

############################ pQTL and AD #####################################

# Harmonize between pQTL and AD

H_data <- harmonise_data(exposure_dat = pQTL_exposure, outcome_dat = Alzheimer_outcome)

#how much the MR-Egger intercept is non-zero
pleiotropy <- mr_pleiotropy_test(H_data)

write.csv(pleiotropy, "PIRLA_pQTL_Pleiotropy.csv")

#Q statistics for heterogeneity
heterogeneity <- mr_heterogeneity(H_data, method_list=c("mr_egger_regression", "mr_ivw"))

write.csv(heterogeneity, "heterogeneity_pQTL_PIRLA.csv")

H_data$Fstat <- (2*H_data$eaf.exposure*(1 - H_data$eaf.exposure) * H_data$beta.exposure^2)*(31684 - 2)/(1 - (2*H_data$eaf.exposure*
                   (1 - H_data$eaf.exposure) * H_data$beta.exposure^2))

write.csv(H_data, "./IVs_PILRA_pQTLs_Alzheimer.csv")


############### SINGLE SNP ANALYSIS ########################################

# since there are only two IVs, we perfromed Single SNP MR from the get-go. 
# This is because with two SNPs the only estimate that can be calucalted is 
# IVW, which the mr_singlesnp() function includes in the analysis anyways

res_single <- mr_singlesnp(H_data)
res_single <- generate_odds_ratios(res_single)

write.csv(res_single, file="res_single_pQTL_Alzheimer.csv")

res_single <- res_single[-4,] # Remove non-existant MR Egger estimate

##### Create figure 21E

# Start with middle part

figure_middle_single <- 
  res_single |> 
  ggplot(aes(y = fct_rev(SNP))) + # have the variants
  theme_classic() +
  geom_point(aes(x=or), shape=15, size=3) + # Add OR values into the plot
  geom_linerange(aes(xmin=or_lci95, xmax=or_uci95)) + # Add CI values into the plot
  geom_vline(xintercept = 1, linetype="dashed") + # Add a vertical line at x = 1
  labs(x="", y="") +
  coord_cartesian(ylim=c(1,4), xlim=c(0.97, 1.05)) +
  annotate("text", x = -.32, y = 5, label = "") +
  theme(axis.line.y = element_blank(), # taking out the axis lines and text 
        axis.ticks.y= element_blank(),
        axis.text.y= element_blank(),
        axis.title.y= element_blank())
figure_middle_single

## Table for plot

table_for_plot_single <- res_single |>
  # round estimates and 95% CIs to 2 decimal places 
  mutate(across(
    c(or, or_lci95, or_uci95),
    ~ str_pad(
      round(.x, 2),
      width = 4,
      pad = "0",
      side = "right"
    )
  ),
  # add an "-" between HR estimate confidence intervals
  OR_CIs = paste0(or, " (", or_lci95, "-", or_uci95, ")")) |>
  # round p-values to two decimal places, except in cases where p < .001
  mutate(p = case_when(
    p < .001 ~ "<0.001",
    round(p, 2) == .05 ~ as.character(round(p,3)),
    p < .01 ~ str_pad( # if p < .01, add one more decimal place
      as.character(round(p, 3)),
      width = 4,
      pad = "0",
      side = "right"
    ),
    TRUE ~ str_pad( # otherwise just round to 2 decimal places and pad string
      as.character(round(p, 2)),
      width = 4,
      pad = "0",
      side = "right"
    )
  )) |>
  # Create a row with column names which (shown on the plot)
  bind_rows( 
    data.frame(
      SNP = "SNP",
      OR_CIs = "OR (95% CI)",
      Lower_CI = "",
      Upper_CI = "",
      p = "p"
    )) |>
  mutate(SNP = fct_rev(fct_relevel(SNP, "SNP")))

table_for_plot_single$p[1] <- 1

# LEFT SIDE OF PLOT

left_figure_single <-
  table_for_plot_single  |>
  ggplot(aes(y = SNP)) + #grab variants
  geom_text(aes(x = 0, label = c("rs12113113", "rs1859788", 
                                 "All- IVW", "SNP")), hjust = 0, 
            fontface = "bold", size = 3.5) +
  geom_text(
    aes(x = 1, label = OR_CIs), # Add OR and CIs
    hjust = 0,
    fontface = ifelse(table_for_plot_single$OR_CIs == "OR (95% CI)", "bold", "plain"), # Create column header
    size = 3.5) +
  theme_void() +
  coord_cartesian(xlim = c(0, 4))
left_figure_single


## RIGHT SIDE OF PLOT

right_figure_single <-
  table_for_plot_single  |>
  ggplot() +
  geom_text(
    aes(x = 0, y = SNP, label = p), # add p value
    hjust = 0,
    fontface = ifelse(table_for_plot_single$p == "p", "bold", "plain"),
    size = 3.5
  ) +
  theme_void() 

right_figure_single


## MERGE PLOTs

# Create layout

layout <- c(
  patchwork::area(t = 0, l = 0, b = 30, r = 10), 
  patchwork::area(t = 0, l = 6, b = 30, r = 14), 
  patchwork::area(t = 0, l = 14, b = 30, r = 15)
)

## CReate final plot, arrange, and save

plot_single <- left_figure_single + figure_middle_single + right_figure_single + plot_layout(design = layout) +
  ggtitle("PILRA IVs cis-pQTLs and Alzheimer MR results ") +
  theme(plot.title = element_text(size = 14,
                                  hjust = 1.35, 
                                  face = "bold",
                                  margin = margin(t= 0, l = 0, b = 15, r= 40)))

ggsave(plot_single, filename = "mr_singlesnp_plot.pdf", device = "pdf", width=9, height=4)


