######################## SELECTION STATISTICS #########################

## SELECT SAMPLE POPULATION TO WORK WITH

# some selection tests only support biallelic variants, not multiallelic. 
# This filtered VCF should already be biallelic SNPs only.
# Note that the zarr file should have been made from a phased vcf
# You make the zarr file with allel.vcf_to_zarr('phased_vcf_file_name.vcf.gz', 'output_name.zarr', fields='*', overwrite=True)

# %%
import os
os.chdir('/mnt/storage11/sophie/bijagos_mosq_wgs/2019_melas_fq2vcf_gambiae_aligned/genomics_database_melas2019plusglobal/genomics_database_melas2019plusglobal_vcf/melas_2019_plusglobal_filtering')
os.getcwd()

# %%
## have now phased VCFs using beagle

import numpy as np
import allel
import zarr
import pandas as pd
import gffutils

## convert phased, filtered, VCF file to zarr file
# %% already converted to zarr
#allel.vcf_to_zarr('2019melasglobal_finalfiltered_gambiaealigned_phased.vcf.gz', '2019melasglobal_finalfiltered_gambiaealigned_phased.zarr', fields='*', overwrite=True)

# %%
callset = zarr.open('2019melasglobal_finalfiltered_gambiaealigned_phased.zarr', mode='r')
#callset.tree(expand=True)

# %%
## convert zarr file to genotype array
gt = allel.GenotypeDaskArray(callset['calldata/GT'])
print(gt.shape)

# %%
## import metadata
df_samples=pd.read_csv('metadata_melasplusglobal.csv',sep=',',usecols=['sample','country','year','species','island'])
df_samples.head()
df_samples.groupby(by=['country']).count

# %%
## working with Guinea-Bissau samples

sample_ids = callset['samples'][:]
# Get sample identifiers for Cameroon samples from df_samples
gb_sample_ids = df_samples[df_samples['country'] == 'Guinea-Bissau']['sample'].values
# Find indices of these samples in the genotype array
gb_indices = np.array([np.where(sample_ids == id)[0][0] for id in gb_sample_ids if id in sample_ids])
# Verify the indices are within the correct range
print("Max index:", gb_indices.max(), "Sample array size:", len(sample_ids))
# Select genotypes for Cameroon samples using the indices
gt_gb_samples = gt.take(gb_indices, axis=1)
gt_gb_samples

# %%
## select variants that are segregating within gb_samples as only these will be informative
## also some selection tests don't support multiallelic variants, so just keep biallelics
## for this pipeline the VCF is already filtered so should be no biallelic SNPs anyway

ac_gb = gt_gb_samples.count_alleles(max_allele=8).compute()
gb_seg_variants = ac_gb.is_segregating() & ac_gb.is_biallelic_01()
ac_gb_seg = ac_gb.compress(gb_seg_variants, axis=0)
gt_gb_seg = gt_gb_samples.compress(gb_seg_variants, axis = 0)
gt_gb_seg

# %%
## this is from a phased VCF so we can convert this genotype array to haplotype array

h_gb_seg = gt_gb_seg.to_haplotypes().compute()
h_gb_seg

# %%
# we need variant positions
pos = callset['variants/POS'][:]
pos_gb_seg = pos.compress(gb_seg_variants, axis=0)
pos_gb_seg

# %%
# some variants in 1000 genomes project have multiple variants at the same genomic position, 
# which causes problems for some selection tests in scikit-allel. 
# Let's check if there any of these.
count_multiple_variants = np.count_nonzero(np.diff(pos_gb_seg == 0))

if count_multiple_variants == 0:
    print("No cases where there are multiple variants at the same genomic position, script will continue")
else:
    print("There are multiple variants at the same genomic position. This causes problems with some selection tests using sci-kit allel.")
    sys.exit()  # This will stop the script. If you want the script to continue anyway, # out this line

# %%
# compute raw iHS

ihs_gb_raw = allel.ihs(h_gb_seg, pos_gb_seg, use_threads=True, include_edges=True)
ihs_gb_raw

# %%

%matplotlib inline
import matplotlib.pyplot as plt
from datetime import datetime

# %%
# ~np.isnan(ihs_gb_std[0]) is used to filter out NaN values
fig, ax = plt.subplots()
ax.hist(ihs_gb_raw[~np.isnan(ihs_gb_raw)], bins=20)
ax.set_xlabel('Raw IHS')
ax.set_ylabel('Frequency (no. variants)');

# %% Standardize iHS

ihs_gb_std = allel.standardize_by_allele_count(ihs_gb_raw, ac_gb_seg[:, 1])
ihs_gb_std

# %% Generate timestamp with current date and time for the figure you are about to make
timestamp = datetime.now().strftime('%Y-%m-%d_%H-%M-%S')

# Here we deviate from the Jupyter notebook and use ihs_res_std[0]
# ~np.isnan(ihs_gb_std[0]) is used to filter out NaN values
fig, ax = plt.subplots()
ax.hist(ihs_gb_std[0][~np.isnan(ihs_gb_std[0])], bins=20)
ax.set_xlabel('Standardised IHS')
ax.set_ylabel('Frequency (no. variants)');

# Save the figure as a file (e.g., PNG) in the current working directory
filename = f'standardised_ihs_histogram_{timestamp}.png'
plt.savefig(filename)

# show the plot (optional, could # this out)
plt.show()

# %% note that iHS has been calculated with unpolarized data, so only the magnitude of iHS
# is informative, not the sign.

# plot over the genome
# np.abs is converting all iHS vales to their absoltue values before plotting. This means that if ihs_gb_std[0]
# contains any negative valyes, those values will be made positive. It is plotting the magnitude of iHS without considering the
# direction of selection, which the sign of iHS could indicate

fig, ax = plt.subplots(figsize=(10, 3))
ax.plot(pos_gb_seg, np.abs(ihs_gb_std[0]), linestyle=' ', marker='o', mfc='none', mew=.5, mec='k')
ax.set_xlabel('Genomic position (bp) chromosome agnostic')
ax.set_ylabel('$|IHS|$')
ax.set_ylim(-2, 9);

# Save the figure as a file (e.g., PNG) in the current working directory
timestamp = datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
filename = f'ihs_manhattan_{timestamp}.png'
plt.savefig(filename)

# %% find the index of the variant with the highest iHS value
idx_hit_max = np.nanargmax(ihs_gb_std[0])

# %% genomic position of top hit
pos_gb_seg[idx_hit_max]
print(f'Genomic position with highest iHS value (chr agnostic):', pos_gb_seg[idx_hit_max])

# %% Visualise EHH decay around top hit with a Voight plot
# pull out haplotypes for region around top hit
flank_size = 2000
h_hit = h_gb_seg[idx_hit_max - flank_size:idx_hit_max + flank_size + 1]
h_hit

# %%
fig = allel.fig_voight_painting(h_hit[:, h_gb_seg[idx_hit_max] == 0], index=flank_size, height_factor=0.02)
fig.suptitle('Reference allele', y=1);

# %% 
fig = allel.fig_voight_painting(h_hit[:, h_gb_seg[idx_hit_max] == 1], index=flank_size, height_factor=0.02)
fig.suptitle('Alternate allele', y=1);

# %%
# Display which iHS values are significant

# Default red line for significance for iHS is 4. Default for XP-EHH is 5. Defauly for RSB is 5.
# I am using iHS significanc of 5 because otherwise it is too many SNPs at 4.

# Include red line in the plot showing significance

### PLOT 1 INCORRECT ###

# fig, ax = plt.subplots(figsize=(10, 3))
# ax.plot(pos_gb_seg, np.abs(ihs_gb_std[0]), linestyle=' ', marker='o', mfc='none', mew=.5, mec='k')
# ax.axhline(y=5, color='red', linestyle='--')
# ax.set_xlabel('Genomic position (bp)')
# ax.set_ylabel('$|IHS|$')
# ax.set_ylim(0, 9)
# ax.legend()
# 
# # Disable scientific notation for x-axis so that full numbers are printed
# ax.get_xaxis().get_major_formatter().set_scientific(False)
# ax.legend()
# 
# # save figure
# timestamp = datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
# filename = f'standardised_ihs_histogram_with_cutoff_{timestamp}.png'
# plt.savefig(filename)
# 
### Need to fix plot, this is not plotting correctly across the genome, it's going from 0 upwards for each chromosome so they're overlapping
# %% PLOT 2 ATTEMPT###

# %% define chromosome lengths and colours
# length of gb_seg_variants = 6785669, it is a numpy.ndarray
chromosomes = callset['variants/CHROM'][:]
chrom_gb_seg = chromosomes.compress(gb_seg_variants, axis = 0)
chrom_gb_seg
# length of chrom_gb_seg = 5673213, it is a numpy.ndarray
pos = callset['variants/POS'][:]
pos_gb_seg = pos.compress(gb_seg_variants, axis=0)
pos_gb_seg
# length of pos_gb_seg = 5673213, it is a numpy.ndarray

# %% 
chromosome_lengths = {
    '2L': 49364325,
    '2R': 61545105,
    '3L': 41963435,
    '3R': 53200684,
    'anop_mito': 15363,
    'anop_X': 24393108
}

# %% Calculate cumulative offsets for each chromosome
cumulative_lengths = {}
cumulative_length = 0
for chrom, length in chromosome_lengths.items():
    cumulative_lengths[chrom] = cumulative_length
    cumulative_length += length

# %% Define colors for each chromosome
chromosome_colors = {
    '2L': 'red',
    '2R': 'blue',
    '3L': 'green',
    '3R': 'orange',
    'anop_mito': 'purple',
    'anop_X': 'cyan'
}

# %% Assuming `pos`, `chromosomes`, and `ihs_gb_std[0]` arrays are already defined and aligned
# Define the threshold
threshold = 5

# Set up the plot
fig, ax = plt.subplots(figsize=(10, 6))

# Define colors for each chromosome (for illustration)
chromosome_colors = {
    '2L': 'red', '2R': 'blue', '3L': 'green', '3R': 'orange', 'anop_mito': 'purple', 'anop_X': 'cyan'
}

# Filter out 'Y_unplaced' or any other chromosomes not in your chromosome_lengths dictionary
filtered_chroms = [chrom for chrom in sorted(set(chrom_gb_seg)) if chrom in chromosome_lengths]

# Iterate through each chromosome to plot its variants
for chrom in filtered_chroms:
    mask = chrom_gb_seg == chrom
    chrom_positions = pos_gb_seg[mask]
    chrom_ihs_values = ihs_gb_std[0][mask]
    adjusted_positions = chrom_positions + cumulative_lengths[chrom]

    below_threshold_mask = chrom_ihs_values < threshold
    above_threshold_mask = chrom_ihs_values >= threshold

    ax.scatter(adjusted_positions[below_threshold_mask], chrom_ihs_values[below_threshold_mask], color=chromosome_colors[chrom], alpha=0.5, s=10)
    ax.scatter(adjusted_positions[above_threshold_mask], chrom_ihs_values[above_threshold_mask], color=chromosome_colors[chrom], alpha=1.0, s=10)

ax.set_xlabel('Adjusted Genomic Position')
ax.set_ylabel('$|IHS|$')
ax.axhline(y=threshold, color='black', linestyle='--', label='Significance Threshold')
ax.set_ylim(-2, max(ihs_gb_std[0]) + 1)
ax.legend(title='Chromosome', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()
plt.show()


############


















# %% list all positions with iHS value over certain threshold
threshold = 5
mask_above_threshold = ihs_gb_std[0] >= threshold
positions_above_threshold = pos_gb_seg[mask_above_threshold]
ihs_values_above_threshold = ihs_gb_std[0][mask_above_threshold]

chromosomes = callset['variants/CHROM'][:]
chromosomes_above_threshold = chromosomes[gb_seg_variants][mask_above_threshold]

# Save positions and corresponding iHS values above the threshold to a text file
# Adjusted file writing to include chromosome, position, and iHS value
with open(f"gb_ihs_positions_above_threshold_{threshold}.txt", "w") as file:
    for chrom, position, ihs_value in zip(chromosomes_above_threshold, positions_above_threshold, ihs_values_above_threshold):
        file.write(f"{chrom}\t{position}\t{ihs_value}\n")

# %% bring in the gff file to understand where each of these variants is

gff_file = '/mnt/storage11/sophie/reference_genomes/A_gam_P4_ensembl/Anopheles_gambiae.AgamP4.56.chr.gff3'
db_file = 'annotations.db'
# create database from gff file (already created annotations.db now)
# db = gffutils.create_db(gff_file, dbfn=db_file, force=True, keep_order=True, merge_strategy='merge', sort_attribute_values=True)

# Querythe database to identify genes associated with the significant positions
# Initialize an empty list to store the positions
positions_above_threshold = []

# Open and read the file
with open("gb_ihs_positions_above_threshold_5.txt", "r") as file:
    for line in file:
        parts = line.strip().split("\t")
        # Assuming the format is: chromosome, position, iHS value
        chromosome, position = parts[0], int(parts[1])
        positions_above_threshold.append((chromosome, position))

# Connect to the database
db = gffutils.FeatureDB(db_file)

# Open a new file to write the annotated positions
with open("gb_ihs_positions_above_threshold_5_gff_annotated.txt", "w") as outfile:
    # Write the header line
    outfile.write("Chromosome\tPosition\tiHS Value\tOverlapping Genes\n")
    
    # Read the original file again to keep the iHS values this time
    with open("gb_ihs_positions_above_threshold_5.txt", "r") as infile:
        for line in infile:
            parts = line.strip().split("\t")
            chromosome, position, ihs_value = parts[0], parts[1], parts[2]
            overlapping_genes = db.region(seqid=chromosome, start=int(position), end=int(position))
            
            # Collect all overlapping gene IDs
            gene_ids = [gene.id for gene in overlapping_genes]
            
            # Write to outfile
            outfile.write(f"{chromosome}\t{position}\t{ihs_value}\t{','.join(gene_ids)}\n")

# %% ########### Cross-population extended haplotype homozygosity (XPEHH) ###########

## Compute the unstandardized cross-population extended haplotype homozygosity score (XPEHH) for each variant.
## allel.xpehh(h1, h2, pos, map_pos=None, min_ehh=0.05, include_edges=False, gap_scale=20000, max_gap=200000, is_accessible=None, use_threads=True)
# create h1 and h2, selecting all variants instead of segregating variants only, which is what we did in iHS


# %%
## VCF is phased so we can convert genotype arrays made earlier to haplotype array
## Create arrays needed for Cameroon samples
sample_ids = callset['samples'][:]
# Get sample identifiers for Cameroon samples from df_samples
cam_sample_ids = df_samples[df_samples['country'] == 'Cameroon']['sample'].values
# Find indices of these samples in the genotype array
cam_indices = np.array([np.where(sample_ids == id)[0][0] for id in cam_sample_ids if id in sample_ids])
# Verify the indices are within the correct range
print("Max index:", cam_indices.max(), "Sample array size:", len(sample_ids))
# Select genotypes for Cameroon samples using the indices
gt_cam_samples = gt.take(cam_indices, axis=1)

# %% select variants that are segregating within cam_samples as only these will be informative
## also some selection tests don't support multiallelic variants, so just keep biallelics

ac_cam = gt_cam_samples.count_alleles(max_allele=8).compute()
cam_seg_variants = ac_cam.is_segregating() & ac_cam.is_biallelic_01()
ac_cam_seg = ac_cam.compress(cam_seg_variants, axis=0)
gt_cam_seg = gt_cam_samples.compress(cam_seg_variants, axis = 0)
gt_cam_seg

# %% this is from a phased VCF so we can convert this genotype array to haplotype array

h_cam_seg = gt_cam_seg.to_haplotypes().compute()
h_cam_seg

# %% we need variant positions
pos = callset['variants/POS'][:]
pos_cam_seg = pos.compress(cam_seg_variants, axis=0)
pos_cam_seg

# Let's check if there any of genomic positions with multiple variants.
count_multiple_variants = np.count_nonzero(np.diff(pos_cam_seg == 0))

if count_multiple_variants == 0:
    print("No cases where there are multiple variants at the same genomic position, script will continue")
else:
    print("There are multiple variants at the same genomic position. This causes problems with some selection tests using sci-kit allel.")
    sys.exit()  # This will stop the script. If you want the script to continue anyway, # out this line

# %% Continue with xp-ehh 




h_sus = gt_sus_samples.to_haplotypes().compute()
h_sus

h_res = gt_res_samples.to_haplotypes().compute()
h_res

ac_gt = gt.count_alleles(max_allele=8).compute()

# %%
# get variant positions

pos = callset['variants/POS'][:]
pos

# %% look at shapes of arrays
print("h_sus shape:", h_sus.shape)
print("h_res shape:", h_res.shape)
print("pos shape", pos.shape)

# %% compute xpehh
# xpehh_raw = allel.xpehh(h_sus, h_res, pos, map_pos=None, min_ehh=0.05, include_edges=False, gap_scale=20000, max_gap=20000, is_accessible=None, use_threads=True)

# xpehh_raw = allel.xpehh(h_sus, h_res, pos, map_pos=None, include_edges=True, use_threads=True)
xpehh_raw = allel.xpehh(h_sus, h_res, pos, use_threads=True)
xpehh_raw

# %% look for where the biggest signal is
xpehh_hit_max = np.nanargmax(xpehh_raw)
xpehh_hit_max

# %% genomic position of top hit
pos[xpehh_hit_max]

# %%
%matplotlib inline
import matplotlib.pyplot as plt

# %%

fig, ax = plt.subplots()
ax.hist(xpehh_raw[~np.isnan(xpehh_raw)], bins=20)
ax.set_xlabel('Raw XP-EHH')
ax.set_ylabel('Frequency (no. variants)');

# %% Standardize XP-EHH - do not think that we really need to do this

# xpehh_std = allel.standardize_by_allele_count(xpehh_raw, ac_gt[:, 1])

# %% 
#fig, ax = plt.subplots()
#ax.hist(xpehh_std[0][~np.isnan(xpehh_std[0])], bins=20)
#ax.set_xlabel('Standardized XP-EHH')
#ax.set_ylabel('Frequency (no. variants)');

# %% note that iHS has been calculated with unpolarized data, so only the magnitude of iHS
# is informative, not the sign.
# xpehh_std

# %% look at shapes

print("pos shape:", pos.shape)
print("xpehh_raw shape:", xpehh_raw.shape)

min_pos = pos.min()
max_pos = pos.max()

print("Minimum Genomic Position:", min_pos)
print("Maximum Genomic Position:", max_pos)


# %% plot on manhattan plot

fig, ax = plt.subplots(figsize=(10, 3))
ax.plot(pos, np.abs(xpehh_raw), linestyle=' ', marker='o', mfc='none', mew=3, mec='k', label='$|XP-EHH|$')
ax.axhline(y=4, color='red', linestyle='--')
ax.set_xlabel('Genomic position (bp)')
ax.set_ylabel('$|XP-EHH|$')
ax.set_ylim(0, 8)
ax.set_xlim(0,61542681)
ax.legend()

# %% list all positions with iHS value over certain threshold

xpehh_positions_above_threshold_4 = pos[xpehh_raw >= 4]

# Save positions_above_threshold to a text file
with open("xpehh_positions_above_threshold_4", "w") as file:
    for position in xpehh_positions_above_threshold_4:
        file.write(str(position) + "\n")



# %% ################################ TAJIMA'S D #####################################

# allel.moving_delta_tajima_d(ac1, ac2, size, start=0, stop=None, step=None
# compute the difference in Tajima's D between two populations in moving windows

# create allele counts arrays ac1 and ac2
# genotype arrays were made earlier in this script: gt_sus_samples and gt_res_samples

ac_sus_samples = gt_sus_samples.count_alleles()
ac_sus_samples

ac_res_samples = gt_res_samples.count_alleles()
ac_res_samples
