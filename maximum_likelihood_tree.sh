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
raxml-ng --all --msa filename.fa --model GTR --prefix Anopheles_mito_0823 --seed 706046 --bs-metric tbe --tree rand{1} --bs-trees 1000

# More species including outgroup
raxml-ng --all --msa filename.fa --model GTR --prefix Anopheles_mito_0823 --seed 706046 --bs-metric tbe --tree rand{1} --bs-trees 1000

