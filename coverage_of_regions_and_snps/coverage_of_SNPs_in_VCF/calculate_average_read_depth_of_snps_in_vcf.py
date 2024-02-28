import gzip

def get_sample_names(vcf_file_path):
    with gzip.open(vcf_file_path, 'rt') as vcf:
        for line in vcf:
            if line.startswith('#CHROM'):
                parts = line.strip().split('\t')
                return parts[9:]  # Sample names are in the columns after the FORMAT field
    return []

def calculate_dp_per_sample(vcf_file_path, chrom, pos, sample_names):
    dp_values = {sample: 'NA' for sample in sample_names}  # Initialize DP values for each sample as 'NA'
    with gzip.open(vcf_file_path, 'rt') as vcf:
        for line in vcf:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            vcf_chrom, vcf_pos = parts[0], int(parts[1])
            if chrom == vcf_chrom and pos == vcf_pos:
                format_parts = parts[8].split(':')
                dp_index = format_parts.index('DP') if 'DP' in format_parts else None
                if dp_index is not None:
                    for sample_idx, sample_data in enumerate(parts[9:]):
                        sample_parts = sample_data.split(':')
                        dp_value = sample_parts[dp_index] if sample_parts[dp_index] != '.' else 'NA'
                        dp_values[sample_names[sample_idx]] = dp_value
                break
    return dp_values

def process_positions(input_file_path, vcf_file_path, output_file_path):
    sample_names = get_sample_names(vcf_file_path)
    with open(input_file_path, 'r') as positions_file, open(output_file_path, 'w') as output_file:
        header = "Chrom\tPos\tAverage_Read_Depth\t" + "\t".join(sample_names) + "\n"
        output_file.write(header)
        for line in positions_file:
            chrom, pos = line.strip().split()
            pos = int(pos)
            dp_values = calculate_dp_per_sample(vcf_file_path, chrom, pos, sample_names)
            average_dp = 'NA' if all(dp == 'NA' for dp in dp_values.values()) else sum(int(dp) for dp in dp_values.values() if dp != 'NA') / len([dp for dp in dp_values.values() if dp != 'NA'])
            output_file.write(f"{chrom}\t{pos}\t{average_dp}\t" + "\t".join(dp_values[sample] for sample in sample_names) + "\n")

# Specify your file paths
input_file_path = 'positions_calculate_average_read_depth.txt'
vcf_file_path = 'melas_2019_plusglobal.2023_07_25.genotyped.vcf.gz'
output_file_path = 'missense_snps_calculated_average_read_depth.txt'

# Run the process
process_positions(input_file_path, vcf_file_path, output_file_path)
