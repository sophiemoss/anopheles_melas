# working directory
/mnt/storage11/sophie/bijagos_mosq_wgs/anopheles_tree

# For mitochondrial tree
# Align mitochondrial fasta files to each other, vcf2fasta script
# muscle -in seqs.fa -out seqs.afa
# use aliview to view, and trim the alignment 
# raxml-ng --all --msa Anopheles_mito_2023_aligned.fa --model GTR --prefix Anopheles_mito_0823 --seed 706046 --bs-metric tbe --tree rand{1} --bs-trees 1000
#  

## First need to generate fasta file of mitochondrial sequences
## Take unfilered vcf for gambiae-gambiae and melas-gambiae, combine, then filter, then use bcftools to subset for gene of interest.
## source unfilered VCFs, both of these vcfs aligned to A_gam_P4_ensembl/Anopheles_gambiae.AgamP4.dna.toplevel.fa

/mnt/storage11/sophie/bijagos_mosq_wgs/malariagen_wgs/gambiae_malariagen_GB_GM-ABC_Bijagos_merged.vcf.gz
/mnt/storage11/sophie/bijagos_mosq_wgs/2019_melas_fq2vcf_gambiae_aligned/genomics_database_melas2019plusglobal/genomics_database_melas2019plusglobal_vcf/melas_2019_plusglobal.2023_07_25.genotyped.vcf.gz

## make sure they are infexed
## copy in the tbi files too because they take forever to index.

## merge VCFs

bcftools merge gambiae_malariagen_GB_GM-ABC_Bijagos_merged.vcf.gz melas_2019_plusglobal.2023_07_25.genotyped.vcf.gz -Oz -o bijagos_gambiae_melas_malariagen_GB_GM-ABC_alignedtoAgamP4_merged.vcf.gz

## filter VCF using filter script
## create the samples_40_10.txt file to make sure you are happy with which samples are going in

bash /mnt/storage11/sophie/gitrepos/anopheles_melas/filter_combined_melas_gambiae_malariagen_vcf.sh

## create just mitochondria VCF
bcftools view -r Mt bijagos_gambiae_melas_malariagen_GB_GM-ABC_AgamP4_aligned_phased.vcf.gz -Oz -o mito_only_bijagos_gambiae_melas_malariagen_GB_GM-ABC_AgamP4_aligned_phased.vcf.gz

## use vcf2fasta to make just fasta sequences. must input multisample vcf
python /mnt/storage11/sophie/fastq2matrix/scripts/vcf2fasta_noiupac.py \
        --vcf "VGSC_only_2022gambiae_phased.vcf.gz" \
        --ref "/mnt/storage11/sophie/reference_genomes/A_gam_P4_ensembl/Mt_only_Anopheles_gambiae.AgamP4.dna.toplevel.fa" \
        --threads 10 \
        --whole-genome