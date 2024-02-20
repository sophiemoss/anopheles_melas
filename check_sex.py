# Calvulate the coverage of Mt and 3R chromosomes for each sample
#ls *bam > bamlist.txt

import subprocess
import pandas as pd

# Function to read BAM file names from bamlist.txt
def read_bam_list(filename):
    with open(filename, 'r') as file:
        bamfiles = file.read().splitlines()
    return bamfiles

# Function to calculate average coverage using samtools depth
def calculate_average_coverage(bamfile, chromosome):
    command = f"samtools depth -r {chromosome} {bamfile} | awk '{{sum+=$3; count++}} END {{if (count > 0) print sum/count; else print 0}}'"
    process = subprocess.run(command, shell=True, capture_output=True, text=True)
    average_coverage = float(process.stdout.strip())
    return average_coverage

# Read BAM files from bamlist.txt
bamfiles = read_bam_list('bamlist.txt')

# Initialize an empty list to store results
results = []

# Process each BAM file
for bamfile in bamfiles:
    # Calculate average coverage for chromosomes X and 3R
    x_coverage = calculate_average_coverage(bamfile, 'X')
    r3_coverage = calculate_average_coverage(bamfile, '3R')

    # Calculate coverage ratio
    if r3_coverage > 0:  # Prevent division by zero
        coverage_ratio = x_coverage / r3_coverage
    else:
        coverage_ratio = 0

    # Determine sex based on coverage ratio
    if 0.4 <= coverage_ratio <= 0.6:
        sex = 'Male'
    elif 0.8 <= coverage_ratio <= 1.2:
        sex = 'Female'
    else:
        sex = 'Excluded'

    # Append results
    results.append([bamfile, sex])

# Convert results to a DataFrame and save to a file
df = pd.DataFrame(results, columns=['BamFile', 'Sex'])
df.to_csv('sex_from_bam.txt', index=False, sep='\t')

print('Finished processing. Results saved to sex_from_bam.txt.')
