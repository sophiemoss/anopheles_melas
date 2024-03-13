## Create vcf of samples to check nucleotide diversity

bcftools query -l final_filteredvcf_bu1003_SRR567658_F6_removed_renamedchr_melas2019plusglobal.vcf.gz | wc -l

## generate nucleotide diversity metric

# for multiple vcfs
#for i in $(ls *SRR.vcf.gz);do vcftools --gzvcf $i --window-pi 10000 --out ${i}_nuc_div_window_10kb;done

vcftools --gzvcf bijagos_only_melas_phased.vcf.gz --window-pi 10000 --out bijagos_only_nuc_div_window_10kb

# or just to get mean nucleotide diversity for population

vcftools --gzvcf bijagos_only_melas_phased.vcf.gz --site-pi --out bijagos_only_nuc_div

# calcualte mean of per site pi

awk '{ total += $3; count++ } END { print total/count }' bijagos_only_nuc_div.sites.pi

0.120334

# calculate mean value of windowed pi

awk '{ total += $5; count++ } END { print total/count }' bijagos_only_nuc_div_window_10kb.windowed.pi

0.0029906

# see which is more commonly reported and what they mean
