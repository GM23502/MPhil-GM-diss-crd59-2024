### ALIGN reads to cis-eQTLs

library(tidyverse)

# Create a list with all the diseases

diseases <- c("ALS", "Alzheimer", "CJD", "LBD", "MS", "Parkinson", "PSP")

# Create a for loop in which this is done in every chromosome individually
# This way it puts less strain on the computer

for(chrom in 1:22){

  
 # Read eQTL data for that chromosome
  
chr_eQTL <- read.table(paste0("./../eQTLs/Raw_eQTLs_per_chrom/chr", chrom, "_whole_blood_eQTLs.txt"), header = TRUE, sep = "\t")

# in each disease
for(disease in diseases) {
  
  # Read the GWAS for that disease in that chromosome
chr_disease <- read.table(paste0("./", disease, "/", disease, "_chr_", chrom, ".tsv"), header = TRUE, sep = "\t")

  # Merge the cis-eQTL data with the disease GWAS
merge_eQTL_disease <- merge(chr_eQTL, chr_disease, by = "SNP")

  # Remove duplicates
merge_eQTL_disease <- merge_eQTL_disease[!duplicated(merge_eQTL_disease[,c(1,2,3,4,5)]),]

  # Remove snp if neither the bases reported in a variant in the disease GWAS is
  # are different from the ones reported in cis-eQTL dataset (e.g., if cis-eQTL
  # had A/G at avariant and CJD had T/C)
merge_eQTL_disease <- merge_eQTL_disease[!(merge_eQTL_disease$A1.x != merge_eQTL_disease$A1.y &
                                   merge_eQTL_disease$A1.x != merge_eQTL_disease$A2.y),]
merge_eQTL_disease <- merge_eQTL_disease[!(merge_eQTL_disease$A2.x != merge_eQTL_disease$A2.y & 
                                  merge_eQTL_disease$A2.x != merge_eQTL_disease$A1.y),]


  # Grab one variant at a time in which the ref and alt alleles are switched (e.g., cis-eQTL
  # has A/T and CJD has T/A) and Switch the order of ref and alt alleles in the disease GWAS data
for(i in 1:nrow(merge_eQTL_disease)){
  if(merge_eQTL_disease$A1.x[i] != merge_eQTL_disease$A1.y[i] & 
    merge_eQTL_disease$A1.x[i] == merge_eQTL_disease$A2.y[i] &
    merge_eQTL_disease$A2.x[i] == merge_eQTL_disease$A1.y[i]) {
        merge_eQTL_disease$A2.y[i] <- merge_eQTL_disease$A2.x[i]
        merge_eQTL_disease$A1.y[i] <- merge_eQTL_disease$A1.x[i]
        merge_eQTL_disease$B[i] <- -merge_eQTL_disease$B[i]
  }
}

  # Remove unnecessary columns
merge_eQTL_disease <- merge_eQTL_disease[,c(1,15:ncol(merge_eQTL_disease))]


  # Format
if(ncol(merge_eQTL_disease) == 10){
  colnames(merge_eQTL_disease) <- c("SNP", "CHR", "BP", "A1", "A2", "EAF", "B", "SE", "P", "N")
}else if(ncol(merge_eQTL_disease) == 9){
  colnames(merge_eQTL_disease) <- c("SNP", "CHR", "BP", "A1", "A2", "EAF", "B", "SE", "P")
}

# Save
write.table(merge_eQTL_disease, file = paste0("Final_aligned_per_chrom/", disease, "/", 
                   disease, "_chr", chrom, ".tsv"), sep = "\t", col.names = TRUE) 

}
}