# Download necessary packages to run the MR analyses and plot 

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("annotate")

install.packages("TwoSampleMR", repos = c("https://mrcieu.r-universe.dev", "https://cloud.r-project.org"))
install.packages("ggplot2")
install.packages("tidyverse")
install.packages("genetics")
install.packages("gt")
install.packages("patchwork")
install.packages("ggeasy")
