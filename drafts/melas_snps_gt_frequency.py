

cd /mnt/storage11/sophie/bijagos_mosq_wgs/2019_melas_fq2vcf_gambiae_aligned/genomics_database_melas2019plusglobal/genomics_database_melas2019plusglobal_vcf/melas_2019_plusglobal_filtering/ir_snps


# %%
import subprocess
import pandas as pd
import subprocess

# %% Define the input file with positions (CHR and POS)
vcf_filename = "FMISSING_MAF_AC0_DP5_GQ20_gatk_filtered_miss_40_mac_bi_snps_melas_2019_plusglobal.2023_07_25.genotyped.vcf.gz"  # Replace with your VCF file path

# %% Define the input and output files
positions_file = "genes_of_interest.txt"

# %% Define the genotype options
genotypes_of_interest = {"0/0", "0/1", "1/0", "1/1", "0|0", "0|1", "1|0", "1|1", "./."}

# %% Prepare a list to collect data
data = []

# %% Read positions from the input file
with open(positions_file, "r") as infile:
    next(infile) # Skip header line
    for line in infile:
        columns = line.strip().split("\t")
        if len(columns) >= 4:
            gene, chr, start, stop = columns

        # Define the bcftools command
        cmd = f'bcftools query -r {chr}:{start}-{stop} -f"%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t%INFO\t[%SAMPLE:%GT]\n" {vcf_filename}'

        # Execute the bcftools command
        try:
            result = subprocess.run(cmd, stdout=subprocess.PIPE, text=True, check=True, shell=True)
            bcftools_output = result.stdout.splitlines()

            for line in bcftools_output:
                fields = line.split('\t')
                position = fields[1]
                sample_genotypes = fields[9:]

                # Count the genotypes
                genotype_counts = {gt: sum(1 for sample_gt in sample_genotypes if sample_gt.split(':')[1] == gt) for gt in genotypes_of_interest}

                # Add the data to the list
                data.append([gene, position] + [genotype_counts[gt] for gt in sorted(genotypes_of_interest)])

        except subprocess.CalledProcessError as e:
            print(f"Error processing gene: {gene}, CHR: {chr}, range: {start}-{stop}. Error: {e}")

# %% Create DataFrame
df = pd.DataFrame(data, columns=["Gene", "Position"] + sorted(genotypes_of_interest))

# %% Optionally, save the DataFrame to a file
output_file = "filtered_vcf_snp_info.csv"
df.to_csv(output_file, index=False)

print("DataFrame created and saved.")

# %%
