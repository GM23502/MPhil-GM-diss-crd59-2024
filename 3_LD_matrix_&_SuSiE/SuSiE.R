## SUSIE chr 19 genes

if (!require("dplyr", quietly = TRUE))
  install.packages("dplyr")

if (!require("susieR", quietly = TRUE))
  install.packages("susieR")

if (!require("ggplot2", quietly = TRUE))
  install.packages("ggplot2")

library(dplyr)
library(susieR)
library(ggplot2)


# list of genes

genes <- c("ARID4B", "CD34", "CLSTN1", "CTSS", "DRAXIN", "ECE1", "EFNA1", "EPHA10", "F11R", "FCRL5",
           "FGR", "GBP4", "LBR", "LGALS8", "ENO1", "NMNAT1", "PADI4", "PARK7", "PMVK", "PRDX1", "RHOC",
           "SLC16A1", "TDRKH", "TNFRSF1B", "TNFRSF8", "TNFRSF9", "TXLNA", "AK2", "AMIGO1", "BOLA1",
           "CD164L2", "CGN", "DNM3", "EIF4G3", "FCER1A", "GBA", "GBP1", "KCNC4", "LDLRAP1", "LMOD1",
           "LYPLA2", "PLEKHO1", "PPT1", "RRP15", "SNAPIN", "SV2A", "TP73", "TTF2")

# Run it on every gene

for(gene in genes){
  
  # Read LD matrix 
  ld_matrix <- read.table(paste0("./", gene, "/", gene,"_ld.ld"), quote="\"", comment.char="")
  ld_matrix <- as.matrix(ld_matrix)
  
  # Read cis-eQTLs for that region
  processed_file <- read.csv(paste0("./", gene, "/", gene, "_ld.processed"), sep="")
  
  # Run SuSiE fine-mapping
  fitted <- susie_rss(processed_file$Z, R = ld_matrix, L = 10, n = 31684)
  
  # Save variants with the highest posterior inclusion probability
  best <- summary(fitted)$cs
  write.table(best, paste0("./", gene, "/", gene, "_fitted_cs.tsv"), sep = "\t", col.names = TRUE)

}