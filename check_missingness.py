import subprocess
import pandas as pd
import argparse

# Set up argparse to handle command line arguments
parser = argparse.ArgumentParser(description='Calculate missingness from a VCF file.')
parser.add_argument('vcf_file', type=str, help='Path to the VCF file')

# Parse the command line arguments
args = parser.parse_args()

# Use the VCF file path from the command line argument
vcf_file_path = args.vcf_file

# Extract sample names from the VCF header
header_command = f"bcftools query -l {vcf_file_path}"
header_process = subprocess.run(header_command, shell=True, capture_output=True, text=True)
sample_names = header_process.stdout.strip().split('\n')

# Command to extract missing genotypes using bcftools
bcftools_command = f"bcftools query -f '%CHROM\t%POS\t[%GT\t]\n' {vcf_file_path}"

# Execute bcftools command and capture output
process = subprocess.run(bcftools_command, shell=True, capture_output=True, text=True)

# Split the output by lines and filter for missing genotypes (./.)
missing_entries = [line for line in process.stdout.split('\n') if './.' in line]

# Extract the sample index for missing genotypes
missing_samples = []
for entry in missing_entries:
    fields = entry.split('\t')
    for i, field in enumerate(fields[2:], start=1):
        if field == './.':
            missing_samples.append(sample_names[i-1])  # Append sample name

# Count missing entries per sample name
missing_counts = pd.Series(missing_samples).value_counts()

# Count total number of variants per sample
# Assuming one variant per line after header lines in the VCF
total_variants = sum(1 for line in open(vcf_file_path) if not line.startswith('#'))

# Calculate the proportion of missingness per sample
missing_proportion = missing_counts / total_variants

# Print and save the missing proportion results
with open('proportion_missingness.txt', 'w') as output_file:
    output_file.write("Sample Name\tMissing Proportion\n")
    print("Sample Name\tMissing Proportion")
    for sample_name, proportion in missing_proportion.items():
        output_line = f"{sample_name}\t{proportion}\n"
        print(output_line.strip())
        output_file.write(output_line)
