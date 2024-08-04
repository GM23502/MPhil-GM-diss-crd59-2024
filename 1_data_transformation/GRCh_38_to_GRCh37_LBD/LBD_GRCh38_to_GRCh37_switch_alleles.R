
## From "out_LBD_hg19.vcf.unmap" create "LBD_switch_alleles_GRCh38.vcf" by switching allele order

setwd("~/Desktop/Cambridge/Dissertation/LBD GRCh38 to GRCh37")

# Read "out_LBD_hg19.vcf.unmap"

out_LBD_hg19.vcf <- read.delim("~/Desktop/Cambridge/Dissertation/LBD GRCh38 to GRCh37/out_LBD_hg19.vcf.unmap",
                               header=FALSE, comment.char="#")

#Switch allele order, add headers, remove last column (explanation as to why it didnt work)

out_LBD_hg19.vcf <- out_LBD_hg19.vcf[,-9]
out_LBD_hg19.vcf <- out_LBD_hg19.vcf[,c(1,2,3,5,4,6,7,8)]

colnames(out_LBD_hg19.vcf) <- c("CHROM", "POS", "ID", "REF", "ALT")

# Create "LBD_switch_alleles_GRCh38.vcf"

write.table(out_LBD_hg19.vcf, file = "LBD_switch_alleles_GRCh38.vcf", quote = FALSE, row.names = FALSE, sep = "\t")

# Add the "#" manually at the left of "CHROM" afterwards