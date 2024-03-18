## download additional samples from ncbi using fasterq-dump
#conda activate fasterq-dump
#fasterq-dump SRR000001

### Put all of the samples in one directory for making database
### Put reference genome in the same directory

## to make a genomics database of sample VCFs, use the following

ls *.g.vcf.gz | sed 's/.g.vcf.gz//' > database_samples.txt #removed samples that do not pass the cut off

## make a list of the sample names that pass the cutoff to be included in the database
## Use merge_vcfs.py import

/mnt/storage11/sophie/fastq2matrix/scripts/merge_vcfs.py import --sample-file database_samples.txt --ref VectorBase-62_AmelasCM1001059_A_Genome.fasta --prefix melas_global_melas_bijagos_2019 --vcf-dir . --threads 25

## now merge VCF files
## use merge_vcfs.py genotype

/mnt/storage11/sophie/fastq2matrix/scripts/merge_vcfs.py genotype --ref VectorBase-62_AmelasCM1001059_A_Genome.fasta --prefix melas_global_melas_bijagos_2019 --threads 25 > melas_alignedmelas_mergevcf_log.txt 2>&1

# resulting vcf should be called melas_global_melas_bijagos_2019.genotyped.vcf.gz
# initial database and vcf made here: /mnt/storage11/sophie/bijagos_mosq_wgs/2019_melas_fq2vcf_melas_aligned
# so many files that cannot ls this folder
# when you ls for *_2019.genotyped.vcf.gz you see that multiple merged vcf files have been made


# concerned that this genotyping does not use a set of validated variants like I had for plasmodium
# this has a genotyping pipeline:
# https://github.com/malariagen/pipelines/blob/v0.0.4/docs/specs/snp-genotyping-vector.md
# pipeline being sent from MalariaGEN.

## Because the melas reference genome was so large and there were so many files, Jody helped me to run the merge_vcfs.py genotype separately:

bedtools makewindows -n 20 -g VectorBase-62_AmelasCM1001059_A_Genome.fasta.fai | awk '{print "melas_global_melas_bijagos_2019_"$1"_"$2+1"_"$3".genotyped.vcf.gz"}' > vcfs.txt
bcftools concat -f vcfs.txt -Oz -o merged.vcf.gz --threads 20

# cleanup step
cat vcfs.txt | xargs rm


### I used the anopheles gambiae database not the melas database for my paper ###