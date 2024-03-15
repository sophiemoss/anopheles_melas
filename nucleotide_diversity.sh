## Create vcf of samples to check nucleotide diversity

bcftools query -l final_filteredvcf_bu1003_SRR567658_F6_removed_renamedchr_melas2019plusglobal.vcf.gz | wc -l

## generate nucleotide diversity metric

# for multiple vcfs
#for i in $(ls *SRR.vcf.gz);do vcftools --gzvcf $i --window-pi 10000 --out ${i}_nuc_div_window_10kb;done

vcftools --gzvcf bijagos_only_melas_phased.vcf.gz --window-pi 10000 --out bijagos_only_nuc_div_window_10kb

bcftools view -r 3L,3R bijagos_only_melas_phased.vcf.gz -Oz -o chrom_3_bijagos_only_melas_phased.vcf.gz

# or just to get mean nucleotide diversity for population

vcftools --gzvcf bijagos_only_melas_phased.vcf.gz --site-pi --out bijagos_only_nuc_div

vcftools --gzvcf chrom_3_bijagos_only_melas_phased.vcf.gz --window-pi 20000 --out chrom_3_bijagos_only_nuc_div_window_20kb
awk '{if(NR>1) sum+=$5; count++} END {print sum/(count-1)}' chrom_3_bijagos_only_nuc_div_window_20kb.windowed.pi (0.00309722)

## median
awk 'NR>1 {print $5}' chrom_3_bijagos_only_nuc_div_window_20kb.windowed.pi | sort -n | awk '{
    count[NR] = $1
}
END {
    if (NR % 2 == 1) {
        print count[(NR + 1) / 2]
    } else {
        median = (count[(NR / 2)] + count[(NR / 2) + 1]) / 2
        print median
    }
}'

(median = 0.003239)

# calculate mean value of per site pi

awk '{ total += $3; count++ } END { print total/count }' bijagos_only_nuc_div.sites.pi

0.120334

# calculate mean value of windowed pi

awk '{ total += $5; count++ } END { print total/count }' bijagos_only_nuc_div_window_10kb.windowed.pi

0.0029906

# see which is more commonly reported and what they mean
