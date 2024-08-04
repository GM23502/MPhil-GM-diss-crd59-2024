## Beak them into chromosomes

library(dplyr)

diseases <- c("ALS", "Alzheimer", "CJD", "LBD", "MS", "Parkinson", "PSP")

for(disease in diseases){
  
  ND <- read.table(paste0(disease, ".tsv"), header = TRUE, sep = "\t")

  for(i in 1:22){

    ND_per_chrom <- ND[ND$CHR == i,]
    write.table(ND_per_chrom, file = paste0(disease, "/", disease, "_chr_", i, ".tsv"), sep = "\t", col.names = TRUE)
}
}