# FILTERING OF COMBINED GENOTYPED VCF
# same filters used for both gambiae pipeline and melas pipeline

########### FILTER 1: remove INDELS with bcftools, here -M2 indicates that bcftools should split multi-allelic sites into multiple biallelic sites, 
# keeping this information
# -m2 is used in conjunction with -M2 to apply the minor-allele-based decomposition. 
# This means that bcftools will decompose multi-allelic sites using the minor allele as the reference allele in the biallelic split.
# here -v tells bcftools to only view SNPS, so indels are excluded
bcftools query -l 2019_merged_melas.vcf.gz

bcftools view -M2 -m2 -v snps 2019_merged_melas.vcf.gz -Oz -o bi_snps_2019_merged_melas.vcf.gz

# tabix index the compressed VCF file, creates .vcf.gz.tbi
tabix -p vcf bi_snps_2019_merged_melas.vcf.gz

# filter out the contigs from the VCF file, note that to produce a properly bgzipped vcf file you need the -Oz flag
# bcftools view bi_snps_2019_merged_melas.vcf.gz --regions 2L,2R,3L,3R,Mt,X,Y_unplaced | bcftools sort -Oz -o bi_snps_chr_2019_merged_melas.vcf.gzz

# view the unique chromosome names
#bcftools query -f '%CHROM\n' bi_snps_chr_gambiae_nov2022.2023_07_05.genotyped.vcf.gz | sort | uniq > unique_chromosomes_filtered.txt
#tabix -p vcf bi_snps_chr_gambiae_nov2022.2023_07_05.genotyped.vcf.gz

########### FILTER 2: Filter samples to keep those with 40% of genome with > 10x coverage, and min-ac=1 so that all variants that remain are still variants after sample removal
# I have identified samples that are above this threshold using my basic WGS stats
# create file with samples to keep: sample_40_10.txt

bcftools view -S samples_40_10.txt bi_snps_2019_merged_melas.vcf.gz --min-ac=1 -Oz -o miss40_mac_bi_snps_2019_merged_melas.vcf.gz

tabix -p vcf miss40_mac_bi_snps_2019_merged_melas.vcf.gz

## count variants, gambiae: 21679588
preVARIANTS=$(bcftools view -H miss40_mac_bi_snps_2019_merged_melas.vcf.gz | wc -l)
echo "The number of variants in the vcf before gatk hard filtering is $preVARIANTS"

########### FILTER 3: for GATK standard filtering
# GATK filter recommendations: https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variants
# excellent explanation of this hard filtering and additional soft filtering
# https://evomics.org/learning/population-and-speciation-genomics/2020-population-and-speciation-genomics/first-steps-in-genomic-data-analysis/#ex3.4

gatk VariantFiltration \
-R /mnt/storage11/sophie/reference_genomes/Amelas_CM1001059_vectorbase/VectorBase-62_AmelasCM1001059_A_Genome.fasta \
-V miss40_mac_bi_snps_2019_merged_melas.vcf.gz \
-filter "QD < 5.0" --filter-name "QD5" \
-filter "QUAL < 30.0" --filter-name "QUAL30" \
-filter "SOR > 3.0" --filter-name "SOR3" \
-filter "FS > 60.0" --filter-name "FS60" \
-filter "MQ < 40.0" --filter-name "MQ40" \
-filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
-filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
-O gatk_tagged_miss40_mac_bi_snps_2019_merged_melas.vcf.gz

## note there is a repeated warning here JEXL engine undefined variable, not a problem
## because this is coming from where there are positions with no coverage
## https://gatk.broadinstitute.org/hc/en-us/community/posts/4408733963803-GATK-Variant-Filtration-undefined-variable

## remove the SNPs that have been tagged with the above filters
bcftools view -f 'PASS' gatk_tagged_miss40_mac_bi_snps_2019_merged_melas.vcf.gz -Oz -o gatk_filtered_miss40_mac_bi_snps_2019_merged_melas.vcf.gz

## count variants, gambiae: 20037931, melas: 8951325
VARIANTS=$(bcftools view -H gatk_filtered_miss40_mac_bi_snps_2019_merged_melas.vcf.gz | wc -l)
echo "The number of variants in the vcf after gatk hard filtering is $VARIANTS" 

## ADDITIONAL SOFT FILTERING

## After hard quality filtering, we have a data set containing only variant sites that we trust with high confidence. 
## However, so far we have only removed variants that had skewed values across all samples – 
## we haven’t yet looked at the individual genotypes at each variant site. 
## For example, we may have kept a variant that has a high quality across all of the individual samples, but there 
## may be some individuals where the read depth (FMT/DP) or the genotype quality (GQ) for the individual genotype is very low, 
## so we don’t trust the calls for these individuals

########### FILTER 4: Filter reads that have a read depth below 5 OR a genotype quality below 20
## The difference is, that instead of the whole variant site, we are now considering single genotypes 
## (the information in the “FORMAT” fields”) using -e 'FMT/DP<3 | FMT/GQ<20' and we will not remove the whole variant site, 
## but only set the respective genotypes to missing (./.) by using the bcftools filter -S . command.
## Here we use the logical operator “|” instead of “||”, since using “||” would mean that every genotype at a variant site 
## is set to missing, even if only one genotype doesn’t fulfill the threshold
bcftools filter -S . -e 'FMT/DP<5 | FMT/GQ<20' -O z -o DP5_GQ20_gatk_filtered_miss40_mac_bi_snps_2019_merged_melas.vcf.gz gatk_filtered_miss40_mac_bi_snps_2019_merged_melas.vcf.gz

## check this, the genotype here should have been set to missing ./. as we are filtering for depth below 5
# bcftools query -i 'FMT/DP<5' -f '[GT=%GT\tDP=%DP\tGQ=%GQ\t]\n' DP5_GQ20_gatk_filtered_miss40_mac_bi_snps_2019_merged_melas.vcf.gz | less -S

## happy? index vcf
tabix -p vcf DP5_GQ20_gatk_filtered_miss40_mac_bi_snps_2019_merged_melas.vcf.gz

## count variants, gambiae: 20037931 melas: 8951325. Same number as all are filters above just set to missing, so genotypes removed but site still in VCF.
countVARIANTS=$(bcftools view -H DP5_GQ20_gatk_filtered_miss40_mac_bi_snps_2019_merged_melas.vcf.gz | wc -l)
echo "The number of variants after filtering FMT/DP<5 | FMT/GQ<20 is $countVARIANTS"

########### FILTER 5: exclude -e all sites at which no alternative alleles are called for any of the samples
bcftools filter -e 'AC==0' DP5_GQ20_gatk_filtered_miss40_mac_bi_snps_2019_merged_melas.vcf.gz -O z -o AC0_DP5_GQ20_gatk_filtered_miss40_mac_bi_snps_2019_merged_melas.vcf.gz

tabix -p vcf AC0_DP5_GQ20_gatk_filtered_miss40_mac_bi_snps_2019_merged_melas.vcf.gz

## count variants, gambiae: 19656287 melas: 8901793

countVARIANTSagain=$(bcftools view -H AC0_DP5_GQ20_gatk_filtered_miss40_mac_bi_snps_2019_merged_melas.vcf.gz | wc -l)
echo "The number of variants after filtering to exclude any non-alternate sites is $countVARIANTSagain"

########### FILTER 6
## Remove variants with a high amount of missing genotypes and filter on minor allele frequency
## The data set overall will now have lots of missing data, because we have replaced calls with depth <5 or quality below 20 with ./.  
## Therefore we will now remove all variants that have more than 20% missing genotypes or MAF < 0.01
## Filtering for MAF 0.01 means that remaining sites in the VCF need to have a minimum minor allele frequency of 1%
## This dataset includes 44 individuals (thus maximally 88 alleles per site). If we filter for MAF 0.01, 
## a variant is excluded if its minor allele is found in less than 0.88 genotypes (1% of 88). 
## This would remove alleles where the minor allele occurs 0.88 times or less, which for this particular dataset means a 
## variant would be removed if there are no minor alleles? So if there is no alternate allele the variant is removed.

bcftools filter -e 'F_MISSING > 0.2 || MAF <= 0.01' -O z -o F_MISSING_MAF_AC0_DP5_GQ20_gatk_filtered_miss40_mac_bi_snps_2019_merged_melas.vcf.gz AC0_DP5_GQ20_gatk_filtered_miss40_mac_bi_snps_2019_merged_melas.vcf.gz

# check this has worked. The minor allele count should be 1 or above in this dataset.
bcftools query F_MISSING_MAF_AC0_DP5_GQ20_gatk_filtered_miss40_mac_bi_snps_2019_merged_melas.vcf.gz -f'%AC\n' | sort -g | head

# count variants, gambiae: 16452859 melas: 8489243
variantcount=$(bcftools view -H F_MISSING_MAF_AC0_DP5_GQ20_gatk_filtered_miss40_mac_bi_snps_2019_merged_melas.vcf.gz | wc -l)
echo "the number of variants after filtering for F_MISSING and MAF is $variantcount"

tabix -p vcf F_MISSING_MAF_AC0_DP5_GQ20_gatk_filtered_miss40_mac_bi_snps_2019_merged_melas.vcf.gz

# Filtering complete!

# removed two samples and then remade vcf, and re-did the above filter (F6)

tabix -p vcf final_filteredvcf_bu1003_SRR567658_F6_removed_renamedchr_melas2019plusglobal.vcf.gz

bcftools view final_filteredvcf_bu1003_SRR567658_F6_removed_renamedchr_melas2019plusglobal.vcf.gz  | grep -v "^#" | wc -l

########## PHASE VCF FILE ###########

# Phase the filtered vcf using beagle, https://faculty.washington.edu/browning/beagle/beagle_5.2_13Oct21.pdf
# necessary for some selection statistical analyses
# could also use phasing pipeline from malariagen. Have a look at this, and understand gatk parameters above.

# use beagle conda environment. mamba install beagle.
conda activate beagle
beagle -Xmx500g gt=final_filteredvcf_bu1003_SRR567658_F6_removed_renamedchr_melas2019plusglobal.vcf.gz out=2019melasglobal_finalfiltered_gambiaealigned_phased

tabix -p vcf 2019melasglobal_finalfiltered_gambiaealigned_phased.vcf.gz

bcftools query -l 2019melasglobal_finalfiltered_gambiaealigned_phased.vcf.gz
bcftools index --stats 2019melasglobal_finalfiltered_gambiaealigned_phased.vcf.gz | cut -f1 | uniq

## Check the number of SNPs in the phased and unphased VCF files for one chromosome, these should be the same
# unphased melas SNPs in chr 2L 1725502
bcftools query -f '%CHROM\t%POS\n' final_filteredvcf_bu1003_SRR567658_F6_removed_renamedchr_melas2019plusglobal.vcf.gz | awk '$1=="2L"' | wc -l
# phased melas SNPs in chr 2L 1725502
bcftools query -f '%CHROM\t%POS\n' 2019melasglobal_finalfiltered_gambiaealigned_phased.vcf.gz | awk '$1=="2L"' | wc -l

# This final vcf contains these chromosome names, with X and Y been renamed.
2L
2R
3L
3R
anop_mito
anop_X
Y_unplaced

# A Mile thesis filtering: criteria: QD < 5; FS > 60; ReadPosRankSum < −8; HRun > 4. I have adjusted QD above to be QD <5. GATK Standard is QD <2. 
# Another phasing pipeline Delaneau et al 2023 - make use of read-backed and statistical phasing methods - seee Ag1000G 2017 paper.

###### snpEFF annotation of filtered VCF #####

## SnpEff annotation of filtered VCF (this one was unphased, could do either)
## Can now use the extract snp info from vcf python script to extract all of the SNP information into readable format

### Aligned bijagos and global melas samples to anopheles gambiae and doing snpEff annotation of those using Anopheles gambiae database ###
snpEff Anopheles_gambiae final_filteredvcf_bu1003_SRR567658_F6_removed_renamedchr_melas2019plusglobal.vcf.gz > final_filteredvcf_bu1003_SRR567658_F6_removed_renamedchr_melas2019plusglobal.ann.vcf.gz
tabix -p vcf final_filteredvcf_bu1003_SRR567658_F6_removed_renamedchr_melas2019plusglobal.ann.vcf.gz



## Making Bijagos only final filtered VCF (for nucleotide diversity and tajimad only)
## create bijagos only vcf

bcftools view -Oz -o bijagos_only_final_filteredvcf_bu1003_SRR567658_F6_removed_renamedchr_melas2019plusglobal.vcf.gz -S bijagos_filtered_samples.txt final_filteredvcf_bu1003_SRR567658_F6_removed_renamedchr_melas2019plusglobal.vcf.gz
tabix -p vcf bijagos_only_final_filteredvcf_bu1003_SRR567658_F6_removed_renamedchr_melas2019plusglobal.vcf.gz

# re-run filter 6
bcftools filter -e 'F_MISSING > 0.2 || MAF <= 0.01' bijagos_only_final_filteredvcf_bu1003_SRR567658_F6_removed_renamedchr_melas2019plusglobal.vcf.gz -Oz -o F_MISSING_MAF_bijagos_only_final_filteredvcf_bu1003_SRR567658_F6_removed_renamedchr_melas2019plusglobal.vcf.gz
tabix -p vcf F_MISSING_MAF_bijagos_only_final_filteredvcf_bu1003_SRR567658_F6_removed_renamedchr_melas2019plusglobal.vcf.gz

#conda activate beagle
beagle -Xmx500g gt=F_MISSING_MAF_bijagos_only_final_filteredvcf_bu1003_SRR567658_F6_removed_renamedchr_melas2019plusglobal.vcf.gz out=bijagos_only_melas_phased

##### 22nd March 2024 #####

## make a vcf for the compkaryo of just 2R and one of just 2L
# make from phased vcf with all samples

bcftools view 2019melasglobal_finalfiltered_gambiaealigned_phased.vcf.gz --regions 2L > 2L_only_2019melasglobal_finalfiltered_gambiaealigned_phased.vcf

bcftools view 2019melasglobal_finalfiltered_gambiaealigned_phased.vcf.gz --regions 2R > 2R_only_2019melasglobal_finalfiltered_gambiaealigned_phased.vcf.gz

# bgzip
# tabix

# to check the chromosome names
# zgrep -v "^#" bijagos_only_melas_phased.vcf.gz | cut -f1 | sort | uniq