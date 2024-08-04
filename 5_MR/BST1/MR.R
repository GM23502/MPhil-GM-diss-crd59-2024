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

### Define a function to create the scatter plots seen in Figure 15B. This code
### was written by the TwoSampleMR team for the mr_scatter_plot() function,
### all I do here is change the aesthetics because I think it looks better that
### way. No credit goes to me for the following code.

trything <- function(mr_results, dat)
{
  # dat <- subset(dat, paste(id.outcome, id.exposure) %in% paste(mr_results$id.outcome, mr_results$id.exposure))
  mrres <- plyr::dlply(dat, c("id.exposure", "id.outcome"), function(d)
  {
    d <- plyr::mutate(d)
    if(nrow(d) < 2 | sum(d$mr_keep) == 0)
    {
      return(blank_plot("Insufficient number of SNPs"))
    }
    d <- subset(d, mr_keep)
    index <- d$beta.exposure < 0
    d$beta.exposure[index] <- d$beta.exposure[index] * -1
    d$beta.outcome[index] <- d$beta.outcome[index] * -1
    mrres <- subset(mr_results, id.exposure == d$id.exposure[1] & id.outcome == d$id.outcome[1])
    mrres$a <- 0
    if("MR Egger" %in% mrres$method)
    {
      temp <- mr_egger_regression(d$beta.exposure, d$beta.outcome, d$se.exposure, d$se.outcome, default_parameters())
      mrres$a[mrres$method == "MR Egger"] <- temp$b_i
    }
    
    if("MR Egger (bootstrap)" %in% mrres$method)
    {
      temp <- mr_egger_regression_bootstrap(d$beta.exposure, d$beta.outcome, d$se.exposure, d$se.outcome, default_parameters())
      mrres$a[mrres$method == "MR Egger (bootstrap)"] <- temp$b_i
    }
    
    ggplot2::ggplot(data=d, ggplot2::aes(x=beta.exposure, y=beta.outcome)) +
      ggplot2::geom_errorbar(ggplot2::aes(ymin=beta.outcome-se.outcome, ymax=beta.outcome+se.outcome), colour="grey", width=0) +
      ggplot2::geom_errorbarh(ggplot2::aes(xmin=beta.exposure-se.exposure, xmax=beta.exposure+se.exposure), colour="grey", height=0) +
      ggplot2::geom_point() +
      ggplot2::geom_abline(data=mrres, ggplot2::aes(intercept=a, slope=b, colour=method), show.legend=TRUE) +
      ggplot2::scale_colour_manual(values=c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a", "#ffff99", "#b15928")) +
      ggplot2::labs(colour="MR Test", x=paste("SNP effect on", d$exposure[1]), y=paste("SNP effect on", d$outcome[1])) +
      ggplot2::guides(colour=ggplot2::guide_legend(ncol=2)) +
      ggplot2::theme_classic() +
      ggplot2::theme(legend.position="top", legend.direction="vertical")
  })
  mrres
}

# Set working directory

setwd("~/Desktop/MPhil-GM-diss-crd59-2024-main/BST1")

# Read exposure dataset

eQTL_exposure <- 
  read_exposure_data("eQTL_exposure.csv",
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

# Harmonize data

H_data <- harmonise_data(exposure_dat = eQTL_exposure, outcome_dat = Parkinson_outcome)

# Calculate first-stage F statistic

H_data$Fstat <- (2*H_data$eaf.exposure*(1 - H_data$eaf.exposure) * H_data$beta.exposure^2)*(31684 - 2)/(1 - (2*H_data$eaf.exposure*
                                                                                                               (1 - H_data$eaf.exposure) * H_data$beta.exposure^2))

write.csv(H_data, "./IVs_eQTL_Parkinson.csv")


#how much the MR-Egger intercept is non-zero

pleiotropy <- mr_pleiotropy_test(H_data)
write.csv(pleiotropy, "./BST1_eQTL_Parkinson_Pleiotropy.csv")


#Q statistics for heterogeneity

heterogeneity <- mr_heterogeneity(H_data, method_list=c("mr_egger_regression", "mr_ivw"))
write.csv(heterogeneity, "./BST1_eQTL_Parkinson_Heterogeneity.csv")


################################ PRIMARY MR ###################################


# Run MR and generate odds ratio

mr_results <- mr(H_data)
mr_results <- generate_odds_ratios(mr_results)

write.csv(mr_results, "./mr_results_eQTL_Parkinson.csv")


#### Create figure 15A

# Create middle part of the plot

figure_middle <- 
  mr_results |> 
  ggplot(aes(y = fct_rev(method))) + # have the MR methods
  theme_classic() +
  geom_point(aes(x=or), shape=15, size=3) + # Add OR values into the plot
  geom_linerange(aes(xmin=or_lci95, xmax=or_uci95)) + # Add CI values into the plot
  geom_vline(xintercept = 1, linetype="dashed") + # Add a vertical line at x = 1
  labs(x="", y="") +
  coord_cartesian(ylim=c(1,6), xlim=c(0.65, 1)) +
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
  # add an "-" between MR estimate confidence intervals
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
      pval = "pval"
    )) |>
  mutate(method = fct_rev(fct_relevel(method, "MR method")))



# LEFT SIDE OF PLOT

left_figure <-
  table_for_plot  |>
  ggplot(aes(y = method)) + #grab variants
  geom_text(aes(x = 0, label = c("MR Egger", "Weighted \n median", 
                                 "IVW", 
                                 "Simple mode", "Weighted \n mode", "MR method")), hjust = 0, 
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
    fontface = ifelse(table_for_plot$pval == "pval", "bold", "plain"),
    size = 3.5
  ) +
  theme_void() 

right_figure


## Merge plots

# Create layout

layout <- c(
  patchwork::area(t = 0, l = 0, b = 30, r = 10), 
  patchwork::area(t = 0, l = 6, b = 30, r = 14), 
  patchwork::area(t = 0, l = 14, b = 30, r = 15)
)

## Create final plot, arrange, and save

plot <- left_figure + figure_middle + right_figure + plot_layout(design = layout) +
  ggtitle("BST1 IVs cis-eQTLs and Parkinson MR results") +
  theme(plot.title = element_text(size = 14,
                                  hjust = 1.35, 
                                  face = "bold",
                                  margin = margin(t= 0, l = 0, b = 15, r= 50)))

ggsave(plot, filename = "mr_plot.pdf", device = "pdf", width=9, height=4)

scatter <- trything(mr_results, H_data)

ggsave(scatter[[1]], file="scatter_plot_eQTLs_Parkinson.pdf", device = "pdf", width=9, height=4)


############### SINGLE SNP ANALYSIS ########################################

# Run single SNP MR and generate odds ratio

res_single <- mr_singlesnp(H_data)
res_single <- generate_odds_ratios(res_single)

write.csv(res_single, file="res_single_eQTL_Parkinson.csv")

##### Create figure 15C

# Middle part of the plot

figure_middle_single <- 
  res_single |> 
  ggplot(aes(y = fct_rev(SNP))) + # have the variants
  theme_classic() +
  geom_point(aes(x=or), shape=15, size=3) + # Add OR values into the plot
  geom_linerange(aes(xmin=or_lci95, xmax=or_uci95)) + # Add CI values into the plot
  geom_vline(xintercept = 1, linetype="dashed") + # Add a vertical line at x = 1
  labs(x="", y="") +
  coord_cartesian(ylim=c(1,8), xlim=c(0.3, 2.3)) +
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


# LEFT SIDE OF PLOT

left_figure_single <-
  table_for_plot_single  |>
  ggplot(aes(y = SNP)) + #grab variants
  geom_text(aes(x = 0, label = c("rs1024757", "rs2302169", 
                                 "rs34559912", "rs35366218", "rs4257648",  
                                 "All- IVW", "All- MR Egger", "SNP")), hjust = 0, 
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
  ggtitle("BST1 IVs cis-eQTLs and Parkinson Single SNP MR results ") +
  theme(plot.title = element_text(size = 14,
                                  hjust = 1.35, 
                                  face = "bold",
                                  margin = margin(t= 0, l = 0, b = 15, r= 0)))

ggsave(plot_single, filename = "mr_singlesnp_plot.pdf", device = "pdf", width=9, height=4)



