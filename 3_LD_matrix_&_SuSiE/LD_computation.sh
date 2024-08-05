#!/bin/bash

for i in ARID4B CD34 CLSTN1 CTSS DRAXIN ECE1 EFNA1 EPHA10 F11R FCRL5 FGR GBP4 LBR LGALS8 ENO1 NMNAT1 PADI4 PARK7 PMVK PRDX1 RHOC SLC16A1 TDRKH TNFRSF1B TNFRSF8 TNFRSF9 TXLNA AK2 AMIGO1 BOLA1 CD164L2 CGN DNM3 EIF4G3 FCER1A GBA GBP1 KCNC4 LDLRAP1 LMOD1 LYPLA2 PLEKHO1 PPT1 RRP15 SNAPIN SV2A TP73 TTF2
do
mkdir $i
python ~/Desktop/PAINTOR/utilities/CalcLD_1KG_VCF.py \  # The command
# The eQTL data for a given gene
--locus ~/Desktop/Cambridge/Dissertation/eQTLs/eQTLs_per_chrom/chr1_eQTLs/"$i"_chr1_eQTLs.hml \
# The reference from the 1000 Genomes project phase 3 panel
--reference ~/Desktop/1000_genomes/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz \
# The map form the 1000 Genomes project phase 3 panel
--map ~/Desktop/1000_genomes/integrated_call_samples_v3.20130502.ALL.panel \
# Effect allele
--effect_allele A1 \
# other allele
--alt_allele A2 \
# population
--population EUR \
# Z score
--Zhead Z \
# Output name
--out_name ./"$i"/"$i"_ld \
# Base pair position
--position BP
done
