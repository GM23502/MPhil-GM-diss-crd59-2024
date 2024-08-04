library(data.table)
library(dplyr)

# Format all GWAS and run quality control

########################### ALS ##############################################

ALS <- read.table("ALS.tsv", sep = "\t", header = TRUE)

colnames(ALS) <- c("SNP", "CHR", "BP", "A1", "A2", "EAF", "B", "SE", "P", "N")

# Remove SNPs with minor allele frequency of <0.05

ALS <- ALS %>% filter(!EAF <= 0.05)
ALS <- ALS %>% filter(!EAF >= 0.95)

length(unique(ALS$SNP)) # Check all SNPs unique
"/" %in% ALS$A1 # Check for non-biallelic variants
"/" %in% ALS$A2 # Check for non-biallelic variants
"/" %in% ALS$A1 # Check for non-biallelic variants
"/" %in% ALS$A2 # Check for non-biallelic variants

ALS <- ALS %>% filter(!N > (mean(ALS$N) + (5 * sd(ALS$N)))) # Remove row with sample number above 5 SD of the mean
ALS <- ALS %>% filter(!N < (mean(ALS$N) - (5 * sd(ALS$N)))) # Remove row with sample number below 5 SD of the mean

# Save

write.table(ALS, "ALS.tsv", sep = "\t", col.names = TRUE)


############################ CJD ##############################################

CJD <- read.table("final_CJD_with_SNPs.tsv", sep = "\t", header = TRUE)

CJD <- CJD[,c(3,1,2,4,5,6,11,12,10,7,8,9)]

colnames(CJD) <- c("SNP", "CHR", "BP", "A1", "A2", "EAF", "B", "SE", "P", "OR", "LCI", "UCI")

# Remove SNPs with minor allele frequency of <0.05

CJD <- CJD %>% filter(!EAF <= 0.05)
CJD <- CJD %>% filter(!EAF >= 0.95)

# Remove duplicated SNPs and SNPs without RSID

CJD <- CJD[!duplicated(CJD$SNP),]
CJD <- CJD[!CJD$SNP == ".",]

"/" %in% CJD$A1 # Check for non-biallelic variants
"/" %in% CJD$A2 # Check for non-biallelic variants
"," %in% CJD$A1 # Check for non-biallelic variants
"," %in% CJD$A2 # Check for non-biallelic variants

# Save

write.table(CJD, "CJD.tsv", sep = "\t", col.names = TRUE)

############################### LBD ############################################

LBD <- read.table("LBD_GWAS_GRCh37.tsv", sep = "\t", header = TRUE)

LBD <- LBD[,c(3,1,2,4,5,7,9,10,6,8)]

colnames(LBD) <- c("SNP", "CHR", "BP", "A1", "A2", "EAF", "B", "SE", "P", "OR")

# Remove SNPs with minor allele frequency of <0.05

LBD <- LBD %>% filter(!EAF <= 0.05)
LBD <- LBD %>% filter(!EAF >= 0.95)

# Remove duplicated SNPs and SNPs without RSID

LBD <- LBD[!duplicated(LBD$SNP),]
LBD <- LBD[!LBD$SNP == ".",]

"/" %in% LBD$A1 # Check for non-biallelic variants
"/" %in% LBD$A2 # Check for non-biallelic variants
"," %in% LBD$A1 # Check for non-biallelic variants
"," %in% LBD$A2 # Check for non-biallelic variants

#Save

write.table(LBD, "LBD.tsv", sep = "\t", col.names = TRUE)

################################ MS ###########################################

MS <- read.table("MS.tsv", sep = "\t", header = TRUE)

MS <- MS[,c(2,1,3,4,5,7,11,12,13,6)]

colnames(MS) <- c("SNP", "CHR", "BP", "A1", "A2", "EAF", "B", "SE", "P", "N")

# Remove SNPs with minor allele frequency of <0.05

MS <- MS %>% filter(!EAF <= 0.05)
MS <- MS %>% filter(!EAF >= 0.95)

# Remove duplicated SNPs and SNPs without RSID

MS <- MS[!duplicated(MS$SNP),]
MS <- MS[!MS$SNP == ".",]

"/" %in% MS$A1 # Check for non-biallelic variants
"/" %in% MS$A2 # Check for non-biallelic variants
"," %in% MS$A1 # Check for non-biallelic variants
"," %in% MS$A2 # Check for non-biallelic variants

MS <- MS %>% filter(!N > (mean(MS$N) + (5 * sd(MS$N)))) # Remove row with sample number above 5 SD of the mean
MS <- MS %>% filter(!N < (mean(MS$N) - (5 * sd(MS$N)))) # Remove row with sample number below 5 SD of the mean

# Save

write.table(MS, "MS.tsv", sep = "\t", col.names = TRUE)

################################### Parkinson #################################

Parkinson <- read.table("final_Parkinson_with_SNPs.tsv", sep ="\t", header = TRUE)

Parkinson <- Parkinson[,c(2,3,1,4,5,7,8,6,9)]

colnames(Parkinson) <- c("SNP", "CHR", "BP", "A1", "A2", "EAF", "B", "SE", "P")

# Remove SNPs with minor allele frequency of <0.05

Parkinson <- Parkinson %>% filter(!EAF <= 0.05)
Parkinson <- Parkinson %>% filter(!EAF >= 0.95)

# Remove duplicated SNPs and SNPs without RSID

Parkinson <- Parkinson[!duplicated(Parkinson$SNP),]
Parkinson <- Parkinson[!Parkinson$SNP == ".",]

"/" %in% Parkinson$A1 # Check for non-biallelic variants
"/" %in% Parkinson$A2 # Check for non-biallelic variants
"," %in% Parkinson$A1 # Check for non-biallelic variants
"," %in% Parkinson$A2 # Check for non-biallelic variants

# Save

write.table(Parkinson, "Parkinson.tsv", sep = "\t", col.names = TRUE)

########################### Alzheimer #########################################

Alzheimer <- read.table("Alzheimer_final.tsv", sep ="\t", header = TRUE)

Alzheimer <- Alzheimer[,1:10]

colnames(Alzheimer) <- c("SNP", "CHR", "BP", "A1", "A2", "EAF", "B", "SE", "P", "N")

# Remove SNPs with minor allele frequency of <0.05

Alzheimer <- Alzheimer %>% filter(!EAF <= 0.05)
Alzheimer <- Alzheimer %>% filter(!EAF >= 0.95)

# Remove duplicated SNPs and SNPs without RSID

Alzheimer <- Alzheimer[!duplicated(Alzheimer$SNP),]
Alzheimer <- Alzheimer[!Alzheimer$SNP == ".",]

"/" %in% Alzheimer$A1 # Check for non-biallelic variants
"/" %in% Alzheimer$A2 # Check for non-biallelic variants
"," %in% Alzheimer$A1 # Check for non-biallelic variants
"," %in% Alzheimer$A2 # Check for non-biallelic variants

MS <- MS %>% filter(!N > (mean(MS$N) + (5 * sd(MS$N)))) # Remove row with sample number above 5 SD of the mean
MS <- MS %>% filter(!N < (mean(MS$N) - (5 * sd(MS$N)))) # Remove row with sample number below 5 SD of the mean


# Save

write.table(Alzheimer, "Alzheimer.tsv", sep = "\t", col.names = TRUE)


