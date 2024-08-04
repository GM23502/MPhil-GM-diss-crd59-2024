The files in this folder have the code necessary to perform quality control on the datasets and align the ref and alt alleles of every variant with the one in the cis-eQTL dataset.

There are three files which shoud be run in the following order:
- Format&QC.R: An initial QC and formatting step for every dataset. 
- break_into_chrom.R: Breaks all the datasets into columns to make the following step less straining for the computer
- align_reads.R: Align the variants in every disease dataset with the ones present in the cis-eQTL dataset and reformat them. This is done one chromosome at a time, and then the last step is to bind all the chromosomes together.

