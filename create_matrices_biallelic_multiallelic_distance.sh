## Convert the filtered vcf to matrix
## to make biallelic matrix use .bi.GT.

# python /mnt/storage11/sophie/fastq2matrix/scripts/vcf2matrix.py --vcf mt_melas_plus_global_subset_filtered.vcf.gz --no-iupacgt --threads 6

## to make multiallelic matrix use .GT.
# not making multiallelic matrix for anopheles, have already split multiallelic sites to be biallelic

## rename mitochondira within the vcf from Mt to anop_mito, otherwise plink removes the variants
vim chrom_map.txt 
Mt  anop_mito
X   anop_X

bcftools annotate --rename-chrs chrom_map.txt FMISSING_MAF_AC0_DP5_GQ20_gatk_filtered_miss_40_mac_bi_snps_melas_2019_plusglobal.2023_07_25.genotyped.ann.vcf.gz -Oz -o pca_filteredvcf_renamedchr_melas2019plusglobal.vcf.gz
tabix -p vcf pca_filteredvcf_renamedchr_melas2019plusglobal.vcf.gz

zgrep -v '^##' pca_filteredvcf_renamedchr_melas2019plusglobal.vcf.gz | cut -f1 | sort | uniq > unique_chromosomes.txt

# remove outlier sample from vcf
bcftools view -s ^bu1003_Combined pca_filteredvcf_renamedchr_melas2019plusglobal.vcf.gz -Oz -o pca_filteredvcf_bu1003removed_renamedchr_melas2019plusglobal.vcf.gz
tabix -p vcf pca_filteredvcf_bu1003removed_renamedchr_melas2019plusglobal.vcf.gz
# run filter 6 again to keep only variants
bcftools filter -e 'F_MISSING > 0.2 || MAF <= 0.01' pca_filteredvcf_bu1003removed_renamedchr_melas2019plusglobal.vcf.gz -Oz -o pca_filteredvcf_bu1003removed_F6_renamedchr_melas2019plusglobal.vcf.gz
tabix -p vcf pca_filteredvcf_bu1003removed_F6_renamedchr_melas2019plusglobal.vcf.gz

## Create a distance matrix for PCA with melas plus global samples
plink --vcf pca_filteredvcf_bu1003removed_F6_renamedchr_melas2019plusglobal.vcf.gz --distance square --double-id --out wholgenome_melas_plusglobal --threads 6 --allow-extra-chr

## now create for every chromosome
## using just the mitochondria
## create mitochondria only vcf using bcftools, then use plink to make distance matrix
bcftools view -r anop_mito pca_filteredvcf_bu1003removed_F6_renamedchr_melas2019plusglobal.vcf.gz -Oz -o mito_only_pca_filteredvcf_bu1003removed_F6_renamedchr_melas2019plusglobal.vcf.gz
tabix -p vcf mito_only_pca_filteredvcf_bu1003removed_F6_renamedchr_melas2019plusglobal.vcf.gz
plink --vcf mito_only_pca_filteredvcf_bu1003removed_F6_renamedchr_melas2019plusglobal.vcf.gz --distance square --double-id --out mito_only_melas_plusglobal --threads 6 --allow-extra-chr

## just the X chromosome
bcftools view -r anop_X pca_filteredvcf_bu1003removed_F6_renamedchr_melas2019plusglobal.vcf.gz -Oz -o X_only_pca_filteredvcf_bu1003removed_F6_renamedchr_melas2019plusglobal.vcf.gz
tabix -p vcf X_only_pca_filteredvcf_bu1003removed_F6_renamedchr_melas2019plusglobal.vcf.gz
plink --vcf X_only_pca_filteredvcf_bu1003removed_F6_renamedchr_melas2019plusglobal.vcf.gz --distance square --double-id --out X_only_melas_plusglobal --threads 6 --allow-extra-chr

## just the 2L chromosome
bcftools view -r 2L pca_filteredvcf_bu1003removed_F6_renamedchr_melas2019plusglobal.vcf.gz -Oz -o 2L_only_pca_filteredvcf_bu1003removed_F6_renamedchr_melas2019plusglobal.vcf.gz
tabix -p vcf 2L_only_pca_filteredvcf_bu1003removed_F6_renamedchr_melas2019plusglobal.vcf.gz
plink --vcf 2L_only_pca_filteredvcf_bu1003removed_F6_renamedchr_melas2019plusglobal.vcf.gz --distance square --double-id --out 2L_only_melas_plusglobal --threads 6 --allow-extra-chr

## just the 2R chromosome
bcftools view -r 2R pca_filteredvcf_bu1003removed_F6_renamedchr_melas2019plusglobal.vcf.gz -Oz -o 2R_only_pca_filteredvcf_bu1003removed_F6_renamedchr_melas2019plusglobal.vcf.gz
tabix -p vcf 2R_only_pca_filteredvcf_bu1003removed_F6_renamedchr_melas2019plusglobal.vcf.gz
plink --vcf 2R_only_pca_filteredvcf_bu1003removed_F6_renamedchr_melas2019plusglobal.vcf.gz --distance square --double-id --out 2R_only_melas_plusglobal --threads 6 --allow-extra-chr

## just the 3L chromosome
bcftools view -r 3L pca_filteredvcf_bu1003removed_F6_renamedchr_melas2019plusglobal.vcf.gz -Oz -o 3L_only_pca_filteredvcf_bu1003removed_F6_renamedchr_melas2019plusglobal.vcf.gz
tabix -p vcf 3L_only_pca_filteredvcf_bu1003removed_F6_renamedchr_melas2019plusglobal.vcf.gz
plink --vcf 3L_only_pca_filteredvcf_bu1003removed_F6_renamedchr_melas2019plusglobal.vcf.gz --distance square --double-id --out 3L_only_melas_plusglobal --threads 6 --allow-extra-chr

## just the 3R chromosome
bcftools view -r 3R pca_filteredvcf_bu1003removed_F6_renamedchr_melas2019plusglobal.vcf.gz -Oz -o 3R_only_pca_filteredvcf_bu1003removed_F6_renamedchr_melas2019plusglobal.vcf.gz
tabix -p vcf 3R_only_pca_filteredvcf_bu1003removed_F6_renamedchr_melas2019plusglobal.vcf.gz
plink --vcf 3R_only_pca_filteredvcf_bu1003removed_F6_renamedchr_melas2019plusglobal.vcf.gz --distance square --double-id --out 3R_only_melas_plusglobal --threads 6 --allow-extra-chr


### Specific gene ###

# creating matrix for VGSC only VCF

python /mnt/storage11/sophie/fastq2matrix/scripts/vcf2matrix.py --vcf VGSC_only_F_MISSING_MAF_AC0_DP5_GQ20_gatk_miss40_mac_bi_snps_gambiae_nov2022.2023_07_05.genotyped.vcf.gz --no-iupacgt --threads 6

# Creating a distance matrix for VGSC PCA 
plink --vcf VGSC_only_F_MISSING_MAF_AC0_DP5_GQ20_gatk_miss40_mac_bi_snps_gambiae_nov2022.2023_07_05.genotyped.vcf.gz --distance square --double-id --allow-extra-chr --out 2019_melas_dist_m --threads 6


### finding unique chromosomes ###

zgrep -v '^##' your_vcf_file.vcf | cut -f1 | sort | uniq > unique_chromosomes.txt
zgrep -v '^##' pca_filteredvcf_renamedchr_melas2019plusglobal.vcf.gz | cut -f1 | sort | uniq > unique_chromosomes.txt

## count variants in vcf
zgrep -v '^##' your_vcf_file.vcf | wc -l
zgrep -v '^##' pca_filteredvcf_bu1003removed_F6_renamedchr_melas2019plusglobal.vcf.gz | wc -l
zgrep -v '^##' X_only_pca_filteredvcf_bu1003removed_F6_renamedchr_melas2019plusglobal.vcf.gz | wc -l

6785670 variants whole genome
474330 variants X only


