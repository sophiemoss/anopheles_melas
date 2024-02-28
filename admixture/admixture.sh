## Admixture

# Identify which samples are in your vcf:
bcftools query -l yourfile.vcf.gz

# STEP 1 EXTRA FILTERING
# MAF > 0.01 filter has already been applied

# STEP 2 CONVERT CHR NAMES IN VCF TO INTEGERS

zgrep -v "^#" final_filteredvcf_bu1003removed_F6_renamedchr_melas2019plusglobal.vcf.gz | cut -f1 | sort | uniq

zcat final_filteredvcf_bu1003removed_F6_renamedchr_melas2019plusglobal.vcf.gz | awk 'BEGIN{OFS=FS="\t"} /^#/ {print; next} {gsub(/^2L$/, "1", $1); gsub(/^2R$/, "2", $1); gsub(/^3L$/, "3", $1); gsub(/^3R$/, "4", $1); gsub(/^anop_mito$/, "5", $1); gsub(/^anop_X$/, "6", $1); gsub(/^Y_unplaced$/, "7", $1); print}' | bgzip > admixture_modified.vcf.gz

tabix -p vcf admixture_modified.vcf.gz
zgrep -v "^#" admixture_modified.vcf.gz | cut -f1 | sort | uniq
bcftools query -l admixture_modified.vcf.gz

# STEP 3 MAKE BED AND BIM FILES

plink --vcf admixture_modified.vcf.gz --set-missing-var-ids @:# --keep-allele-order --const-fid --allow-extra-chr --make-bed --out melas_global_gambiaealigned

# STEP 4 RUN ADMIXTURE

cat K_runs.txt | xargs -I {} sh -c 'admixture --cv=10 -j20 -s 14062 melas_global_gambiaealigned.bed {} | tee log{}.cv10.seed14062.out'

# Plot admixture
conda create -n radmix r-essentials r-base
install.packages(c("unikn", "countrycode", "optparse"))
setwd("/mnt/storage11/sophie/bijagos_mosq_wgs/2019_melas_fq2vcf_gambiae_aligned/genomics_database_melas2019plusglobal/genomics_database_melas2019plusglobal_vcf/melas_2019_plusglobal_filtering/admixture")

Rscript /mnt/storage11/sophie/gitrepos/anopheles_melas/admixture/generate_admix_barplot.R \
-d /mnt/storage11/sophie/bijagos_mosq_wgs/2019_melas_fq2vcf_gambiae_aligned/genomics_database_melas2019plusglobal/genomics_database_melas2019plusglobal_vcf/melas_2019_plusglobal_filtering/admixture \
--prefix melas_global_gambiaealigned \
--kval 4 \
-m metadata_melasplusglobal_admixture.tsv \
--filter_N 1 \
--label_id sample \
--label_region region \
--label_country country \
--label_site site \
--country_code TRUE