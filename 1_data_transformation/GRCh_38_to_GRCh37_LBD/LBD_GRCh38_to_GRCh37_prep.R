#GRCh38 to GRCh37 LBD data preparation: make into a .vcf

setwd("~/Desktop/Cambridge/Dissertation/LBD GRCh38 to GRCh37")

# load LBD data

LBD <- read.delim("~/Desktop/Cambridge/Dissertation/Data/Primary study/LBD.tsv")

# Grab Chromosome, base pair, SNP rsID, A1, A2 in that order (as it is in a vcf file)

LBD_vcf <- LBD[,c(3,4,1,5,6)]

# Add VCF column names

colnames(LBD_vcf) <- c("CHROM", "POS", "ID", "REF", "ALT")

# Create fake QUAL FILTER and INFO columns (not used when changing builds anyways)
# in INFO column add the rownames so as to make the switch easier later

LBD_vcf$QUAL <- "."
LBD_vcf$FILTER <- "."
LBD_vcf$INFO <- rownames(LBD)

# Change NAs in rsids to a dot

LBD_vcf$ID[is.na(LBD_vcf$ID)] <- "."

#write vcf, add a hash at the beginning of the header line manually afterwards

write.table(LBD_vcf, file = "LBD_GRCh38.vcf", quote = FALSE, row.names = FALSE, sep = "\t")

# Add the "#" manually at the left of "CHROM" afterwards