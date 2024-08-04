#### Getting the SNPs for GWAS without SNPs: Parkinson, CJD. 

## Extract CJD and Parkinson GWAS

CJD <- read.delim("~/Desktop/Cambridge/Dissertation/Data/Primary study/CJD.tsv")
Parkinson <- read.delim("~/Desktop/Cambridge/Dissertation/Data/Primary study/Parkinson.tsv")

## See which chromosomes are studied

unique(CJD$chromosome)
unique(Parkinson$chromosome)

## Create table neeeded: Chr, pos, pos, a1, a2

CJD_structured <- CJD[,c(1,2,4,3)]
CJD_structured$BP2 <- CJD_structured$base_pair_location
CJD_structured <- CJD_structured[,c(1,2,5,3,4)]

Parkinson_structured <- Parkinson[,c(1,2,3,4)]
Parkinson_structured$BP2 <- Parkinson_structured$base_pair_location
Parkinson_structured <- Parkinson_structured[,c(1,2,5,3,4)]

## Create per chromosome tsv files for CJD

for(i in unique(CJD_structured$chromosome)){
  chr <- CJD_structured[CJD_structured$chromosome == i,]
  write.table(chr, file = paste0("./CJD per chrom/CJD_chr", i, ".tsv"), quote = FALSE, row.names = FALSE, sep = "\t")
}

## Create per chromosome tsv files for Parkinson

for(i in unique(Parkinson_structured$chromosome)){
  chr <- Parkinson_structured[Parkinson_structured$chromosome == i,]
  write.table(chr, file = paste0("./Parkinson per chr/Parkinson_chr", i, ".tsv"), quote = FALSE, row.names = FALSE, sep = "\t")
}