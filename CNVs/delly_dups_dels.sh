# Use for loops to call duplications and deletions for every sample, and create bcf files (input files for delly)
# Starts with unfiltered bam files

for f in *.mkdup.bam; do delly call -t DUP -g Anopheles_gambiae.AgamP4.dna.toplevel.fa "$f" -o "${f%.*}.dup.bcf"; done
 
for f in *.mkdup.bam; do delly call -t DEL -g Anopheles_gambiae.AgamP4.dna.toplevel.fa "$f" -o "${f%.*}.del.bcf"; done

# and insertions?

# use Jody's script to merge and create a multi-sample vcf of deletions and duplications
# https://jodyphelan.gitbook.io/tutorials/ngs/fst-with-delly

# STEP 1
# create a multi-sample vcf of deletions and a multisample vcf of duplications

delly merge -o deletions_to_genotype.bcf *.del.bcf
delly merge -o duplications_to_genotype.bcf *.dup.bcf

# STEP 2
# re-genotype the bam files to look for specific deletions/duplications in the merged deletions_to_genotype.bcf

ls *.mkdup.bam | sed 's/.mkdup.bam//' > all_samples.txt

cat all_samples.txt | parallel -j 15 --bar "delly call -g Anopheles_gambiae.AgamP4.dna.toplevel.fa -v deletions_to_genotype.bcf -o {}.del.genotyped.bcf  {}.mkdup.bam"
cat all_samples.txt | parallel -j 15 --bar "delly call -g Anopheles_gambiae.AgamP4.dna.toplevel.fa -v duplications_to_genotype.bcf -o {}.dup.genotyped.bcf  {}.mkdup.bam"

# STEP 3 Next step is to merge the genotyped bcfs

bcftools merge -m id -O b -o dup_merged.bcf *.dup.genotyped.bcf
bcftools merge -m id -O b -o del_merged.bcf *.del.genotyped.bcf

# STEP 4 Filtering

bcftools view deletions_to_genotype.bcf -H | wc -l
# total number of deletion variants genotyped 136359

# remove samples with >20% missing data (136359*0.2, rounded up to the nearest 10)
bcftools +smpl-stats del_merged.bcf | grep ^FLT0 | awk '$12<27270'  | cut -f2 > del_good_samples.txt
# ADDITIONAL STEP UNIQUE FOR THIS DATASET: edit del_good_samples.txt to also remove any specific samples that were filtered out during my combined filtering process (to keep samples aligned)

# extract the good samples from the bcf
bcftools view -S del_good_samples.txt del_merged.bcf -Oz -o del_merged.sample_filt.vcf.gz
tabix -p vcf del_merged.sample_filt.vcf.gz

# additionally filter to mask calls that are heterozygous (i.e. have mixed ref/alt calls) and to remove deletions with >20% missing data
bcftools view del_merged.sample_filt.vcf.gz | bcftools filter -e 'GT="het"' -S .  | bcftools view -i 'F_PASS(GT!="mis")>0.8' -Oz -o del_merged.sample_filt.site_filt.vcf.gz
tabix -p vcf del_merged.sample_filt.site_filt.vcf.gz


# STEP 5 filter the VCF based on genes of interest
bcftools view -R genes_of_interest_+-1kb.bed del_merged.sample_filt.site_filt.vcf.gz -Oz -o del_genes_of_interest_merged.sample_filt.site_filt.vcf.gz
tabix -p vcf del_genes_of_interest_merged.sample_filt.site_filt.vcf.gz

# Query number of variants
bcftools query -f '%CHROM\n' del_genes_of_interest_merged.sample_filt.site_filt.vcf.gz | wc -l
# 117 deletions

# Some of these deletions are higher quality than others.

### Analysis Type 1: between Bijagos and Elsewhere Fst ###

# Will now look to see if there are any in the Bijagos melas (pop1) vs non-Bijagos melas (pop2)
# Made two populations.txt files pop1.txt (Bijagos samples) and pop2.txt (cameroon/gambia samples)
vcftools --gzvcf del_genes_of_interest_merged.sample_filt.site_filt.vcf.gz --weir-fst-pop pop1.txt --weir-fst-pop pop2.txt --out Bijagos_vs_Elsewhere

# There are four deletions with positive fst values between the two populations
# CHROM	POS	WEIR_AND_COCKERHAM_FST
# 2L	2367585	0.64433
# 2L	2394887	1
# 2L	25424721	0.2019
# 2L	25428858	0.928783


### Analysis Type 2: Deletions in the Bijagos population that look potentially important ###
# 117 in total in the candidate genes, which ones shall I report on?

### Updated delly version ###

## Trying with delly call ALL

# read the sample identifiers from all_samples.txt, process the .all.bcf files, in parallel using -P 10,
# then convert this to a vcf and then convert it to a text file without the header lines for easy viewing.
# Produces individual text files for each sample

ls *.mkdup.bam | sed 's/.bam//' > delly_all_samples.txt

for f in *.mkdup.bam; do delly call -t ALL -g Anopheles_gambiae.AgamP4.dna.toplevel.fa "$f" -o "${f%.*}.all.bcf"; done

delly merge -o structural_variants.bcf *.all.bcf

bcftools view structural_variants.bcf > structural_variants.vcf

bgzip structural_variants.vcf

## Have a look at the structural_variants.vcf.gz
zgrep -v ^"##" structural_variants.vcf.gz | less

## genotype the structural variants # Delly genotyping requires local SV assembly (INFO/CONSENSUS) and breakpoint (INFO/CONSBP) introduced in delly v1.1.7!

cat delly_all_samples.txt | parallel -j 15 --bar "delly call -g Anopheles_gambiae.AgamP4.dna.toplevel.fa -v structural_variants.bcf -o {}.all.genotyped.bcf  {}.bam"

# %% merge all of the genotyped bcf files together
bcftools merge -m id -O b -o merged.bcf *.all.genotyped.bcf

## Filter - keep samples that were used in downstream analysis - edit delly_good_samples.txt for this (removed samples with over 20% missing data)

bcftools view -S delly_good_samples.txt structural_variants.vcf.gz -Oz -o structural_variants.sample_filt.vcf.gz

## Re-annotate with snpeff



### Delly CNVs?

## creat mappability map

dicey chop Anopheles_gambiae.AgamP4.dna.toplevel.fa
bwa index Anopheles_gambiae.AgamP4.dna.toplevel.fa
bwa mem Anopheles_gambiae.AgamP4.dna.toplevel.fa car1002_Combined_1.fastq.gz car1002_Combined_2.fastq.gz | samtools sort -@ 8 -o srt.bam -
samtools index srt.bam 
dicey mappability2 srt.bam 
gunzip map.fa.gz && bgzip map.fa && samtools faidx map.fa.gz 