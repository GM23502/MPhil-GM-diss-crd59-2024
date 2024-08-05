For this step it is necessary to download the 1000 genomes phase 3 reference panel, available at ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/
It is also necessary to download PAINTOR 3.0, available at https://github.com/gkichaev/PAINTOR_V3.0

There are two sample code: 
- LD_computation.sh: The code necessary to compute LD matrices. The sample code shows the code I ran to compute LD matrices for all genes in chromosome 1. For this one it is necessary to have the cis-eQTL variants for one region, the "integrated_call_samples_v3.20130502.ALL.panel" file from 1000 Genomes Phase 3 panel, and the ALL.chr1.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz" file for the chromosome the region is on, also from the 1000 Genomes Phase 3 panel. Computing the LD matrices is done with the PAINTOR 3 "CalcLD_1KG_VCF.py" utility command.
- SuSiE.R: Sample code that we ran to perform SuSiE fine-mapping, also for all the genes in chromosome 1. Inside this file is also instructions to download all the necessary R packages
