####################### Two sample MR for pQTL and Parkinson ################

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

setwd("~/Desktop/MPhil-GM-diss-crd59-2024-main/HIP1R")

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

Parkinson_outcome <-
  read_outcome_data("Parkinson_outcome.csv",
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

# Perform harmonization

H_data <- harmonise_data(exposure_dat = pQTL_exposure, outcome_dat = Parkinson_outcome)

# Perform first-stage F statistic

H_data$Fstat <- (2*H_data$eaf.exposure*(1 - H_data$eaf.exposure) * H_data$beta.exposure^2)*(31684 - 2)/(1 - (2*H_data$eaf.exposure*
                                                                                                               (1 - H_data$eaf.exposure) * H_data$beta.exposure^2))
# Remove weak variants

H_data <- H_data[!H_data$SNP == "rs10744218",]
H_data <- H_data[!H_data$SNP == "rs7964876",]
H_data <- H_data[!H_data$SNP == "rs4759404",]

write.csv(H_data, "./IVs_pQTL_Parkinson.csv")

#how much the MR-Egger intercept is non-zero
pleiotropy <- mr_pleiotropy_test(H_data)

write.csv(pleiotropy, "./HIP1R_pQTL_Parkinson_Pleiotropy.csv")

#Q statistics for heterogeneity
heterogeneity <- mr_heterogeneity(H_data, method_list=c("mr_egger_regression", "mr_ivw"))

write.csv(heterogeneity, "./HIP1R_pQTL_Parkinson_Heterogeneity.csv")

################################ PRIMARY MR ###################################

# Run MR an

mr_results <- mr(H_data)
mr_results <- generate_odds_ratios(mr_results)

write.csv(mr_results, "./mr_results_pQTL_Parkinson.csv")


#### MR Results plot, Figure 20D

# Create the middle part

figure_middle <- 
  mr_results |> 
  ggplot(aes(y = fct_rev(method))) + # have the variants
  theme_classic() +
  geom_point(aes(x=or), shape=15, size=3) + # Add OR values into the plot
  geom_linerange(aes(xmin=or_lci95, xmax=or_uci95)) + # Add CI values into the plot
  geom_vline(xintercept = 1, linetype="dashed") + # Add a vertical line at x = 1
  labs(x="", y="") +
  coord_cartesian(ylim=c(1,2), xlim=c(0, 1.1)) +
  annotate("text", x = -.32, y = 5, label = "") +
  theme(axis.line.y = element_blank(), # taking out the axis lines and text 
        axis.ticks.y= element_blank(),
        axis.text.y= element_blank(),
        axis.title.y= element_blank())
figure_middle

## Table for plot

table_for_plot <- mr_results |>
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
  mutate(pval = case_when(
    pval < .001 ~ "<0.001",
    round(pval, 2) == .05 ~ as.character(round(pval,3)),
    pval < .01 ~ str_pad( # if p < .01, add one more decimal place
      as.character(round(pval, 3)),
      width = 4,
      pad = "0",
      side = "right"
    ),
    TRUE ~ str_pad( # otherwise just round to 2 decimal places and pad string
      as.character(round(pval, 2)),
      width = 4,
      pad = "0",
      side = "right"
    )
  )) |>
  # Create a row with column names which (shown on the plot)
  bind_rows( 
    data.frame(
      method = "MR method",
      OR_CIs = "OR (95% CI)",
      Lower_CI = "",
      Upper_CI = "",
      pval = "P"
    )) |>
  mutate(method = fct_rev(fct_relevel(method, "MR method")))

# LEFT SIDE OF PLOT

left_figure <-
  table_for_plot  |>
  ggplot(aes(y = method)) + #grab variants
  geom_text(aes(x = 0, label = c("rs12817488 \n Wald ratio", "SNP")), hjust = 0, 
            fontface = "bold", size = 3.5) +
  geom_text(
    aes(x = 1, label = OR_CIs), # Add OR and CIs
    hjust = 0,
    fontface = ifelse(table_for_plot$OR_CIs == "OR (95% CI)", "bold", "plain"), # Create column header
    size = 3.5) +
  theme_void() +
  coord_cartesian(xlim = c(0, 4))
left_figure


## RIGHT SIDE OF PLOT

right_figure <-
  table_for_plot  |>
  ggplot() +
  geom_text(
    aes(x = 0, y = method, label = pval), # add p value
    hjust = 0,
    fontface = ifelse(table_for_plot$pval == "Pval", "bold", "plain"),
    size = 3.5
  ) +
  theme_void() 

right_figure


## MERGE PLOTs

# Create layout

layout <- c(
  patchwork::area(t = 0, l = 0, b = 30, r = 10), 
  patchwork::area(t = 0, l = 6, b = 30, r = 14), 
  patchwork::area(t = 0, l = 14, b = 30, r = 15)
)

## CReate final plot, arrange, and save

plot <- left_figure + figure_middle + right_figure + plot_layout(design = layout) +
  ggtitle("HIP1R IV cis-pQTLs and Parkinson MR results") +
  theme(plot.title = element_text(size = 14,
                                  hjust = 1.35, 
                                  face = "bold",
                                  margin = margin(t= 0, l = 0, b = 15, r= 50)))

ggsave(plot, filename = "mr_plot.pdf", device = "pdf", width=9, height=4)
