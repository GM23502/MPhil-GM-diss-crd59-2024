## Install and load packages

 if (!require("devtools", quietly = TRUE))
  install.packages("devtools")

 library(devtools)

 if (!require("hyprcoloc", quietly = TRUE))
 install_github("cnfoley/hyprcoloc", build_opts = c("--resave-data", "--no-manual"), build_vignettes = TRUE)

library(hyprcoloc)

## Set working directory 

setwd("~/Desktop/MPhil-GM-diss-crd59-2024-main/4_HyPrColoc")

# Create list of genes on which we are going to perform hyprcoloc

genes <- c("BST1", "DDR1", "HIP1R", "PILRA")

# Create for loop in which HyPrColoc analyses are run for every gene

for(gene in genes){
  
  #### Read all datasets necessary
  
  # List of all beta regression coefficient values 
  beta_main <- read.table(paste0(gene, "/beta_main_", gene, ".tsv"), sep = "\t", header = TRUE)
  beta_main <- as.matrix(beta_main)
  
  # List of all standard error values
  se_main <- read.table(paste0(gene, "/se_main_", gene, ".tsv"), sep = "\t", header = TRUE)
  se_main <- as.matrix(se_main)
  
  # list of all snps that are being analysed
  load(paste0(gene, "/snp_list_main.rda"))
  
  # list of the names of all the traits being analysed
  load(paste0(gene, "/trait_list_main.rda"))
  
  # List of whether every trait is continuous (0) or binary (1)
  load(paste0(gene, "/bin_main.rda"))
  
  
  #### Run primary HyPrColoc analysis
  
  # snpscores = TRUE makes it so the posterior probability of a cluster 
  # explained by every SNP is reported
  main_first_hyprcoloc <- hyprcoloc(beta_main, se_main,
                                    trait.names=trait_list_main,
                                    snp.id=snp_list_main, 
                                    binary.outcomes= bin_main,
                                    snpscores = TRUE) 
  
  #Save
  write.table(as.data.frame(main_first_hyprcoloc$results),
              paste0("./", gene, "/main_first_hyprcoloc.tsv"), 
              sep= "\t", col.names = TRUE)
  
  # Create sensitivity analyses in which prior.c, the regional threshold, and 
  # Alignment threshold have different values to assess the stability of the cluster
  # Any cluster that colocalises with PP = 0.80 in over 75% of these configurations
  # is deemed a stable cluster.
  
  sens_plot <- sensitivity.plot(beta_main,
                                se_main, 
                                trait.names = trait_list_main, 
                                snp.id=snp_list_main, 
                   reg.thresh = c(0.6,0.7,0.8,0.9),
                   align.thresh = c(0.6,0.7,0.8,0.9),
                   prior.c = c(0.02, 0.01, 0.005), 
                   equal.thresholds = FALSE, 
                   uniform.priors = FALSE,
                   binary.outcomes= bin_main,
                   similarity.matrix = TRUE)
  
  # Save matrix
  write.table(as.data.frame(sens_plot[[2]]), 
              paste0("./", gene, "/matrix_main.tsv"),
              sep = "\t", col.names = TRUE)
}