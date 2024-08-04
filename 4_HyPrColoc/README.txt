This folder contains everything necessary to run our HyPrColoc analyses on BST1, DDR1, HIP1R, and PILRA. We did not include STX6 or ADAMTS4 because they required the PSP data to be included, and the PSP data is not publicly available.

The "hyprcoloc.R" file contains the code to run all HyPrColoc analyses at once, using a for loop. Importantly, we wrote this using R 4.3.3 and hyprcoloc 1.0

Inside each folder are 5 files:
- beta_main_GENENAME.tsv: This file contains the beta regression coefficient values for every SNP in the 1Mb cis block.
- se_main_GENENAME.tsv: This file contains the standard erorr values for all the SNPs in the 1Mb cis block
- bin_main.rda: This is a file comprised of only 1s and 0s. A 1 means that the trait is defines is binary (all NDs), and a 0 means that the trait it defines is continuous (cis-pQTLs and cis-eQTLS). It is not necessary to use for the HyPrColoc analysis, but we did it because it changes how the algorithm computes colocalisation.
- snp_list_main.rda: A file comprise of a list of SNP RSIDs in the 1Mb cis block
- trait_list_main.rda: A file with the name of every trait. 
