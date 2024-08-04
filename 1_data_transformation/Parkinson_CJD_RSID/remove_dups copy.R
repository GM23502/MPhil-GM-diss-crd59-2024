library(stringr)
library(tidyverse)
library(dplyr)

AGCT_combinations <- c("A,A,A", "A,A,G", "A,A,T", "A,A,C", "A,G,A", "A,G,G", 
                       "A,G,T", "A,G,C", "A,T,A", "A,T,G", "A,T,T", "A,T,C", 
                       "A,C,A", "A,C,G", "A,C,T", "A,C,C", "G,A,A", "G,A,G", 
                       "G,A,T", "G,A,C", "G,G,A", "G,G,G", "G,G,T", "G,G,C",
                       "G,T,A", "G,T,G", "G,T,T", "G,T,C", "G,C,A", "G,C,G", 
                       "G,C,T", "G,C,C", "T,A,A", "T,A,G", "T,A,T", "T,A,C", 
                       "T,G,A", "T,G,G", "T,G,T", "T,G,C", "T,T,A", "T,T,G", 
                       "T,T,T", "T,T,C", "T,C,A", "T,C,G", "T,C,T", "T,C,C",
                       "C,A,A", "C,A,G", "C,A,T", "C,A,C", "C,G,A", "C,G,G", 
                       "C,G,T", "C,G,C", "C,T,A", "C,T,G", "C,T,T", "C,T,C",
                       "C,C,A", "C,C,G", "C,C,T", "C,C,C")

for(i in 1:22){

# Read inferred snps file 

final_CJD_1 <- read.delim(paste0("/gpfs3/well/combat/users/vzd660/diss/Parkinson_final_files/intersect_chr",i,"_Parkinson.tsv"), header=FALSE)

# only keep those with matching base pairs

final_CJD_1 <- final_CJD_1[final_CJD_1$V2 == final_CJD_1$V7,]

# Keep only those where at least one of the alleles is the reference

final_CJD_1 <- final_CJD_1[final_CJD_1$V4 == final_CJD_1$V9 | final_CJD_1$V5 == final_CJD_1$V9,]

# Remove irrelevant columns 

final_CJD_1 <- final_CJD_1[,-c(11,12,13)]

#Colnames

colnames(final_CJD_1) <- c("chromosome", "base_pair_location", "BP2", "effect_allele", "other_allele",
                                "chr_imp", "pos_imp", "snp_imp","a1_imp","a2_imp")

# Grab only base pairs that repeat

final_CJD_1 <- final_CJD_1[nchar(final_CJD_1$a2_imp) == 1 | 
                             (grepl(",", final_CJD_1$a2_imp) & nchar(final_CJD_1$a2_imp) == 3) |
                             final_CJD_1$a2_imp %in% AGCT_combinations,]

# Treat the rest individually

n_occur <- data.frame(table(final_CJD_1$base_pair_location))
only_dups <- final_CJD_1[final_CJD_1$base_pair_location %in% n_occur$Var1[n_occur$Freq > 1],]

write.table(final_CJD_1, file = paste0("./intdir_Parkinson/intfile_chr", i, ".tsv"), row.names = FALSE, sep = "\t")
write.table(only_dups, file = paste0("./intdir_Parkinson/only_dups_chr", i, ".tsv"), row.names = FALSE, sep = "\t")

}

intfile_chr <- read.delim2("/gpfs3/well/combat/users/vzd660/diss/intdir_Parkinson/intfile_chr22.tsv")
only_dups_chr <- read.delim("/gpfs3/well/combat/users/vzd660/diss/intdir_Parkinson/only_dups_chr22.tsv")

length(unique(only_dups_chr$snp_imp))
length(unique(only_dups_chr$base_pair_location))


imp_seq <- seq(from = 2, to = 90, by = 2)

for(i in imp_seq){
intfile_chr <- intfile_chr %>% filter(!snp_imp == only_dups_chr$snp_imp[i])
}

length(unique(intfile_chr$base_pair_location))

intfile_chr <- intfile_chr %>% filter(!snp_imp == "rs1235670236")
intfile_chr <- intfile_chr %>% filter(!snp_imp == "rs9554314")
intfile_chr <- intfile_chr %>% filter(!snp_imp == "rs1347058791")
intfile_chr <- intfile_chr %>% filter(!snp_imp == "rs372225628")
intfile_chr <- intfile_chr %>% filter(!snp_imp == "rs372225628")

#save final file

write.table(intfile_chr, "./Final_SNPs_Parkinson//final_SNP_list_Parkinson_chr22.tsv", row.names = FALSE, sep = "\t")
