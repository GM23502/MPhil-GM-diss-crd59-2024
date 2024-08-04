# Convert the LBD GWAS data to hg19

setwd("~/Desktop/Cambridge/Dissertation/LBD GRCh38 to GRCh37")

# Read GWAS dataset

LBD <- read.delim("~/Desktop/Cambridge/Dissertation/LBD GRCh38 to GRCh37/LBD.tsv")

# Read output dataset with original ref/other allele order

normal_alleles_hg19 <- read.table("~/Desktop/Cambridge/Dissertation/LBD GRCh38 to GRCh37/out_LBD_hg19.vcf", quote="\"")

# Read output dataset with switched ref/other allele order and switch the order of the alleles

switched_alleles_hg19 <- 
  read.table("~/Desktop/Cambridge/Dissertation/LBD GRCh38 to GRCh37/out_ LBD_switch_alleles_hg19.vcf", quote="\"")

switched_alleles_hg19 <- switched_alleles_hg19[,c(1,2,3,5,4,6,7,8)]

# Read output dataset with reads that went unmaped and will be discarded from final dataset
# and remove the final column which only says Fail(unmapped)

unmaped_reads <- 
  read.delim("~/Desktop/Cambridge/Dissertation/LBD GRCh38 to GRCh37/out_LBD_switch_alleles_hg19.vcf.unmap",
             header=FALSE, comment.char="#")

unmaped_reads <- unmaped_reads[,-9]

# Step 1: Create extra column in GWAS dataset with rownames, name the equivalent 
# column in all datasets the same thing and provide column names for all datasets

LBD$row <- row.names(LBD)
colnames(normal_alleles_hg19) <- c("CHROM_hg19", "POS_hg19", "RSID_hg19", "A1_hg19", "A2_hg19", "QUAL", "FILTER", "row")
colnames(switched_alleles_hg19) <- c("CHROM_hg19", "POS_hg19", "RSID_hg19", "A1_hg19", "A2_hg19", "QUAL", "FILTER", "row")
colnames(unmaped_reads) <- c("CHROM_hg19", "POS_hg19", "RSID_hg19", "A1_hg19", "A2_hg19", "QUAL", "FILTER", "row")

# Step 2: Add an identifier to unmaped reads, this will be a "UNMAPPED" on the QUAL column

unmaped_reads$QUAL <- "UNMAPPED"

# Step 3: merge "normal_alleles_hg19", "switched_alleles_hg19", and "unmaped_reads"

Complete_GRCh37_vcf <- rbind(normal_alleles_hg19, switched_alleles_hg19, unmaped_reads)

# Step 4: Merge both data frames using the "row" column as a meeting point

merged_GWAS_and_GRCh37_vcf <- merge.data.frame(LBD, Complete_GRCh37_vcf, by = "row")

# Step 5: Removed unmapped reads

merged_GWAS_and_GRCh37_vcf <- merged_GWAS_and_GRCh37_vcf[!merged_GWAS_and_GRCh37_vcf$QUAL == "UNMAPPED",]

# Step 6: Make final GRCh37 GWAS summary statistics table

final_GRCh37_LBD_GWAS <- merged_GWAS_and_GRCh37_vcf[,c(12, 13, 14, 15, 16, 3, 8, 9, 10, 11)]

# Step 7: Save 

write.table(final_GRCh37_LBD_GWAS, file = "~/Desktop/Cambridge/Dissertation/LBD GRCh38 to GRCh37/LBD_GWAS_GRCh37", 
            col.names = TRUE, row.names = FALSE, sep = "\t")
