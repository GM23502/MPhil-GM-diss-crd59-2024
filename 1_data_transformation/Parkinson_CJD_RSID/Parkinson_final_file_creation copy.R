library(dplyr)

#Actual thing

Parkinson_chr <- read.delim("/gpfs3/well/combat/users/vzd660/diss/Parkinson per chr/Parkinson_chr22.tsv")

final_SNP_list_Parkinson_chr <- 
  read.delim("/gpfs3/well/combat/users/vzd660/diss/Final_SNPs_Parkinson/final_SNP_list_Parkinson_chr22.tsv")

#tests that there is only one entry per base pair position

length(unique(Parkinson_chr$base_pair_location)) == length(Parkinson_chr$chromosome)

dif_vals <- setdiff(Parkinson_chr$base_pair_location, final_SNP_list_Parkinson_chr$base_pair_location)

newrow <- Parkinson_chr[Parkinson_chr$base_pair_location %in% dif_vals,]
newrow$chr_imp <- newrow$chromosome
newrow$pos_imp <- newrow$base_pair_location
newrow$snp_imp <- "."
newrow$a1_imp <- newrow$effect_allele
newrow$a2_imp <- newrow$other_allele

final_SNP_list_Parkinson_chr <- rbind(final_SNP_list_Parkinson_chr, newrow)


for(i in dif_vals){
newrow <- c(unlist(Parkinson_chr[Parkinson_chr$base_pair_location == i,], use.names = FALSE), 
  unlist(Parkinson_chr[Parkinson_chr$base_pair_location == i,1:2], use.names = FALSE),
  ".",
  unlist(Parkinson_chr[Parkinson_chr$base_pair_location == i, 4:5], use.names = FALSE))

final_SNP_list_Parkinson_chr <- rbind(final_SNP_list_Parkinson_chr, newrow)
}

write.table(final_SNP_list_Parkinson_chr,
            "/gpfs3/well/combat/users/vzd660/diss/completely_final_Parkinson/final_Parkinson_SNPs_chr22.tsv",
            row.names = FALSE, sep = "\t")