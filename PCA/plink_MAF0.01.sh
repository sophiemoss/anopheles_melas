## Create a distance matrix for PCA with melas plus global samples
plink --vcf xxx --distance square --double-id --out darlingi --threads 6 --allow-extra-chr

## now create for every chromosome ## no variants remain???
## using just the mitochondria
## create mitochondria only vcf using bcftools, then use plink to make distance matrix
## have to rename Mt to be mito so that plink will work

bcftools view -r mito Complete_F_MISSING_MAF_AC0_DP5_GQ20_gatk_filtered_minac_filtrenamedchr_bi_snps_darlingi.genotyped.vcf.gz -Oz -o mt_only_F_MISSING_MAF_AC0_DP5_GQ20_gatk_filtered_minac_filtrenamedchr_bi_snps_darlingi.genotyped.vcf.gz
tabix -p vcf mt_only_F_MISSING_MAF_AC0_DP5_GQ20_gatk_filtered_minac_filtrenamedchr_bi_snps_darlingi.genotyped.vcf.gz
plink --vcf mt_only_F_MISSING_MAF_AC0_DP5_GQ20_gatk_filtered_minac_filtrenamedchr_bi_snps_darlingi.genotyped.vcf.gz --distance square --double-id --out mito_only_darlingi --threads 6 --allow-extra-chr

## just the X chromosome
bcftools view -r X darlingi_filtered_phased.vcf.gz -Oz -o X_only_darlingi_filtered_phased.vcf.gz
tabix -p vcf X_only_darlingi_filtered_phased.vcf.gz
plink --vcf X_only_darlingi_filtered_phased.vcf.gz --distance square --double-id --out X_only_darlingi --threads 6 --allow-extra-chr

## just the 2L chromosome
bcftools view -r 2L MAF0.02_final_filteredvcf_bu1003_SRR567658_F6_removed_renamedchr_melas2019plusglobal.ann.vcf.gz -Oz -o 2L_only_MAF0.02_final_filteredvcf_bu1003_SRR567658_F6_removed_renamedchr_melas2019plusglobal.ann.vcf.gz
tabix -p vcf 2L_only_MAF0.02_final_filteredvcf_bu1003_SRR567658_F6_removed_renamedchr_melas2019plusglobal.ann.vcf.gz
plink --vcf 2L_only_MAF0.02_final_filteredvcf_bu1003_SRR567658_F6_removed_renamedchr_melas2019plusglobal.ann.vcf.gz --distance square --double-id --out 2L_only_melas_plusglobal --threads 6 --allow-extra-chr

## just the 2R chromosome
bcftools view -r 2R MAF0.02_final_filteredvcf_bu1003_SRR567658_F6_removed_renamedchr_melas2019plusglobal.ann.vcf.gz -Oz -o 2R_only_MAF0.02_final_filteredvcf_bu1003_SRR567658_F6_removed_renamedchr_melas2019plusglobal.ann.vcf.gz
tabix -p vcf 2R_only_MAF0.02_final_filteredvcf_bu1003_SRR567658_F6_removed_renamedchr_melas2019plusglobal.ann.vcf.gz
plink --vcf 2R_only_MAF0.02_final_filteredvcf_bu1003_SRR567658_F6_removed_renamedchr_melas2019plusglobal.ann.vcf.gz --distance square --double-id --out 2R_only_melas_plusglobal --threads 6 --allow-extra-chr

## just the 3L chromosome
bcftools view -r 3L MAF0.02_final_filteredvcf_bu1003_SRR567658_F6_removed_renamedchr_melas2019plusglobal.ann.vcf.gz -Oz -o 3L_only_MAF0.02_final_filteredvcf_bu1003_SRR567658_F6_removed_renamedchr_melas2019plusglobal.ann.vcf.gz
tabix -p vcf 3L_only_MAF0.02_final_filteredvcf_bu1003_SRR567658_F6_removed_renamedchr_melas2019plusglobal.ann.vcf.gz
plink --vcf 3L_only_MAF0.02_final_filteredvcf_bu1003_SRR567658_F6_removed_renamedchr_melas2019plusglobal.ann.vcf.gz --distance square --double-id --out 3L_only_melas_plusglobal --threads 6 --allow-extra-chr

## just the 3R chromosome
bcftools view -r 3R MAF0.02_final_filteredvcf_bu1003_SRR567658_F6_removed_renamedchr_melas2019plusglobal.ann.vcf.gz -Oz -o 3R_only_MAF0.02_final_filteredvcf_bu1003_SRR567658_F6_removed_renamedchr_melas2019plusglobal.ann.vcf.gz
tabix -p vcf 3R_only_MAF0.02_final_filteredvcf_bu1003_SRR567658_F6_removed_renamedchr_melas2019plusglobal.ann.vcf.gz
plink --vcf 3R_only_MAF0.02_final_filteredvcf_bu1003_SRR567658_F6_removed_renamedchr_melas2019plusglobal.ann.vcf.gz --distance square --double-id --out 3R_only_melas_plusglobal --threads 6 --allow-extra-chr
