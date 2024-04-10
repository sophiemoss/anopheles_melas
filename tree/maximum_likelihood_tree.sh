# working directory
/mnt/storage11/sophie/bijagos_mosq_wgs/anopheles_tree

# For mitochondrial tree
# Created fasta files from the filtered melas and filtered gambiae vcfs, which contain bijagos samples (and some extra melas from ncbi)

## use vcf2fasta to make just fasta sequences. must input multisample vcf
python /mnt/storage11/sophie/fastq2matrix/scripts/vcf2fasta_noiupac.py \
        --vcf "mito_only_FMISSING_MAF_AC0_DP5_GQ20_gatk_filtered_miss_40_mac_bi_snps_melas_2019_plusglobal.2023_07_25.genotyped.ann.vcf.gz" \
        --ref "/mnt/storage11/sophie/reference_genomes/A_gam_P4_ensembl/Mt_only_Anopheles_gambiae.AgamP4.dna.toplevel.fa" \
        --threads 10 \
        --whole-genome

# melas fasta:
/mnt/storage11/sophie/bijagos_mosq_wgs/2019_melas_fq2vcf_gambiae_aligned/genomics_database_melas2019plusglobal/genomics_database_melas2019plusglobal_vcf/melas_2019_plusglobal_filtering/mito_sequences/mito_only_FMISSING_MAF_AC0_DP5_GQ20_gatk_filtered_miss_40_mac_bi_snps_melas_2019_plusglobal.2023_07_25.genotyped.ann.fa
# gambiae fasta:
/mnt/storage11/sophie/bijagos_mosq_wgs/2022_gambiae_fq2vcf_agamP4/gambiae_nov2022_genomicdb/gambiae_nov2022_genotypedvcf/gambiae_nov2022_filtering/mitochondrial_sequences/mito_only_F_MISSING_MAF_AC0_DP5_GQ20_gatk_miss40_mac_bi_snps_gambiae_nov2022.2023_07_05.genotyped.ann.fa

# Align mitochondrial fasta files to each other
# muscle -align seqs.fa -output seqs.afa

muscle -align mito_only_F_MISSING_MAF_AC0_DP5_GQ20_gatk_miss40_mac_bi_snps_gambiae_nov2022.2023_07_05.genotyped.ann.fa -output 2022_gambiae_bijagos_mito.afa
muscle -align mito_only_FMISSING_MAF_AC0_DP5_GQ20_gatk_filtered_miss_40_mac_bi_snps_melas_2019_plusglobal.2023_07_25.genotyped.ann.fa -output 2019_melas_global_bijagos_mito.afa

# download alignments to local disk
# use aliview to view, and trim the alignment 
# reupload trimmed alignments to server

# Bijagos only, trim and align in aliview, then run raxml-ng and plot
raxml-ng --all --msa Anopheles_mito_18022024_aligned_bijagos_with_global_melasv6.fa --model GTR --prefix Anopheles_mito_0823 --seed 706046 --bs-metric tbe --tree rand{1} --bs-trees 1000

# More species including outgroup
raxml-ng --all --msa Anopheles_mito_16022024_with_outgroup_aligned.fa --model GTR --prefix Anopheles_mito_0823 --seed 706046 --bs-metric tbe --tree rand{1} --bs-trees 1000

# Adding anopheles arabiensis from Holly
Aarab_Agam.2021_04_02.genotyped.vcf.gz
# Change chromosome names so that they match my reference for GATK filtering later: 2L, 2R, 3L, 3R, X, Mt
# Chromosomes
zgrep -v "^#" Aarab_Agam.2021_04_02.genotyped.vcf.gz | cut -f1 | sort | uniq

# Rename chromosomes
vim chrom_map.txt 
AgamP4_Mt       Mt
AgamP4_X        X
AgamP4_2L       2L
AgamP4_2R       2R
AgamP4_3L       3L
AgamP4_3R       3R

bcftools annotate --rename-chrs chrom_map.txt Aarab_Agam.2021_04_02.genotyped.vcf.gz -Oz -o Aarab_Agam.2021_04_02_renamedchroms.genotyped.vcf.gz
tabix -p vcf Aarab_Agam.2021_04_02_renamedchroms.genotyped.vcf.gz
zgrep -v "^#" Aarab_Agam.2021_04_02_renamedchroms.genotyped.vcf.gz | cut -f1 | sort | uniq

# filter with my filters and make tree
bcftools view -M2 -m2 -v snps Aarab_Agam.2021_04_02_renamedchroms.genotyped.vcf.gz -Oz -o bi_snps_Aarab_Agam.2021_04_02_renamedchroms.genotyped.vcf.gz
tabix -p vcf bi_snps_Aarab_Agam.2021_04_02_renamedchroms.genotyped.vcf.gz

bcftools view --min-ac=1 -Oz -o mac_bi_snps_Aarab_Agam.2021_04_02_renamedchroms.genotyped.vcf.gz bi_snps_Aarab_Agam.2021_04_02_renamedchroms.genotyped.vcf.gz
tabix -p vcf mac_bi_snps_Aarab_Agam.2021_04_02_renamedchroms.genotyped.vcf.gz

gatk VariantFiltration \
-R /mnt/storage11/sophie/reference_genomes/A_gam_P4_ensembl/Anopheles_gambiae.AgamP4.dna.toplevel.fa \
-V mac_bi_snps_Aarab_Agam.2021_04_02_renamedchroms.genotyped.vcf.gz \
-filter "QD < 5.0" --filter-name "QD5" \
-filter "QUAL < 30.0" --filter-name "QUAL30" \
-filter "SOR > 3.0" --filter-name "SOR3" \
-filter "FS > 60.0" --filter-name "FS60" \
-filter "MQ < 40.0" --filter-name "MQ40" \
-filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
-filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
-O gatk_tagged_mac_bi_snps_Aarab_Agam.2021_04_02_renamedchroms.genotyped.vcf.gz

bcftools view -f 'PASS' gatk_tagged_mac_bi_snps_Aarab_Agam.2021_04_02_renamedchroms.genotyped.vcf.gz -Oz -o gatk_filtered_mac_bi_snps_Aarab_Agam.2021_04_02_renamedchroms.genotyped.vcf.gz
tabix -p vcf gatk_filtered_mac_bi_snps_Aarab_Agam.2021_04_02_renamedchroms.genotyped.vcf.gz

bcftools filter -S . -e 'FMT/DP<5 | FMT/GQ<20' -O z -o DP5_GQ20_gatk_filtered_mac_bi_snps_Aarab_Agam.2021_04_02_renamedchroms.genotyped.vcf.gz gatk_filtered_mac_bi_snps_Aarab_Agam.2021_04_02_renamedchroms.genotyped.vcf.gz
tabix -p vcf DP5_GQ20_gatk_filtered_mac_bi_snps_Aarab_Agam.2021_04_02_renamedchroms.genotyped.vcf.gz

bcftools filter -e 'AC==0' DP5_GQ20_gatk_filtered_mac_bi_snps_Aarab_Agam.2021_04_02_renamedchroms.genotyped.vcf.gz -O z -o AC0_DP5_GQ20_gatk_filtered_mac_bi_snps_Aarab_Agam.2021_04_02_renamedchroms.genotyped.vcf.gz
tabix -p vcf AC0_DP5_GQ20_gatk_filtered_mac_bi_snps_Aarab_Agam.2021_04_02_renamedchroms.genotyped.vcf.gz

bcftools filter -e 'F_MISSING > 0.2 || MAF <= 0.01' -O z -o F_MISSING_MAF_AC0_DP5_GQ20_gatk_filtered_mac_bi_snps_Aarab_Agam.2021_04_02_renamedchroms.genotyped.vcf.gz DP5_GQ20_gatk_filtered_mac_bi_snps_Aarab_Agam.2021_04_02_renamedchroms.genotyped.vcf.gz
tabix -p vcf F_MISSING_MAF_AC0_DP5_GQ20_gatk_filtered_mac_bi_snps_Aarab_Agam.2021_04_02_renamedchroms.genotyped.vcf.gz

# subset vcf to be mitochondria only

bcftools view -r Mt F_MISSING_MAF_AC0_DP5_GQ20_gatk_filtered_mac_bi_snps_Aarab_Agam.2021_04_02_renamedchroms.genotyped.vcf.gz -Oz -o mito_only_F_MISSING_MAF_AC0_DP5_GQ20_gatk_filtered_mac_bi_snps_Aarab_Agam.2021_04_02_renamedchroms.genotyped.vcf.gz

# make fasta files
python /mnt/storage11/sophie/fastq2matrix/scripts/vcf2fasta_noiupac.py \
        --vcf "mito_only_F_MISSING_MAF_AC0_DP5_GQ20_gatk_filtered_mac_bi_snps_Aarab_Agam.2021_04_02_renamedchroms.genotyped.vcf.gz" \
        --ref "/mnt/storage11/sophie/reference_genomes/A_gam_P4_ensembl/Mt_only_Anopheles_gambiae.AgamP4.dna.toplevel.fa" \
        --threads 10 \
        --whole-genome

# Both types of tree
raxml-ng --all --msa Anopheles_mito_18022024_aligned_bijagos_with_global_melas_added_arabiensisv6.fa --model GTR --prefix Anopheles_mito_with_arabiensis --seed 706046 --bs-metric tbe --tree rand{1} --bs-trees 1000

raxml-ng --all --msa Anopheles_mito_16022024_with_outgroup_aligned_withdarlingiarabiensis.fa --model GTR --prefix Anopheles_mito_with_darlingi_arabiensis --seed 706046 --bs-metric tbe --tree rand{1} --bs-trees 1000

raxml-ng --all --msa Anopheles_gambiae_sensu_lato.afa --model GTR --prefix Anopheles_gambaie_sensu_lato_mito --seed 725364 --bs-metric tbe --tree rand{1} --bs-trees 1000

# Use best_tree for iTOL