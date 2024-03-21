## Create vcf of samples to check nucleotide diversity

bcftools query -l final_filteredvcf_bu1003_SRR567658_F6_removed_renamedchr_melas2019plusglobal.vcf.gz | wc -l

## generate nucleotide diversity metric

# for multiple vcfs
#for i in $(ls *SRR.vcf.gz);do vcftools --gzvcf $i --window-pi 10000 --out ${i}_nuc_div_window_10kb;done

vcftools --gzvcf bijagos_only_melas_phased.vcf.gz --window-pi 10000 --out bijagos_only_nuc_div_window_10kb

# Can look at chromosome 3L only because it is not impact by large inversion found in An. gambiae complex in 3R, 2R, 2L, or X.

bcftools view -r 3L bijagos_only_melas_phased.vcf.gz -Oz -o 3L_bijagos_only_melas_phased.vcf.gz

# calculate nucleotide diversity per site

vcftools --gzvcf 3L_bijagos_only_melas_phased.vcf.gz --site-pi --out bijagos_3L_only_nuc_div

# calculate mean 

awk '{ total += $3; count++ } END { print total/count }' bijagos_3L_only_nuc_div.sites.pi
# 0.115401 (11.5% of nucleotides are different between individuals in this population on average)

# calculate nucleotide diversity per window

vcftools --gzvcf 3L_bijagos_only_melas_phased.vcf.gz --window-pi 20000 --out bijagos_3L_only_nuc_div_window20k

awk '{if(NR>1) sum+=$5; count++} END {print sum/(count-1)}' bijagos_3L_only_nuc_div_window20k.windowed.pi
# 0.00302611 (0.3% of nucleotides are different between individuals in this population on average)
`
## Tajima's D

vcftools --gzvcf 3L_bijagos_only_melas_phased.vcf.gz --TajimaD 20000 --out 3L_bijagos_only_melas_phased