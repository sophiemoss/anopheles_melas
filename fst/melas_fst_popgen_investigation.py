## scikitallel_workflow
# /mnt/storage11/sophie/gitrepos/anopheles_melas/melas_fst_v1.py
# run this script using
# python /mnt/storage11/sophie/gitrepos/anopheles_melas/melas_fst_v1.py /pathway/to/dir /pathway/to/zarr chromosomename
# python /mnt/storage11/sophie/gitrepos/anopheles_melas/melas_fst_v2.py . /mnt/storage11/sophie/bijagos_mosq_wgs/2019_melas_fq2vcf_gambiae_aligned/genomics_database_melas2019plusglobal/genomics_database_melas2019plusglobal_vcf/melas_2019_plusglobal_filtering/2019melasglobal_finalfiltered_gambiaealigned_phased.zarr 2R
# You make the zarr file with allel.vcf_to_zarr('phased_vcf_file_name.vcf.gz', 'output_name.zarr', fields='*', overwrite=True)

# allel.vcf_to_zarr('2019_melas_phased.vcf.gz', '2019_melas_phased.zarr', fields='*', overwrite=True)
# python /mnt/storage11/sophie/gitrepos/anophelesmelas_popgen/scikit_allel_fst_calculate_and_plot_melasPCAinvestigation.py /mnt/storage11/sophie/bijagos_mosq_wgs/2019_melas_fq2vcf/genomics_database_melas2019plusglobal/genomics_database_melas2019plusglobal_vcf/melas2019plusglobal_vcf_filtering 2019_melas_phased.zarr 2L

#### This script is altered specifically for Mt and X becuase the names in the gff and the zarr file are different. Use the zarr file chr name for the command, then change this line [chromosome_for_gff = chromosome_arg.replace("anop_mito","Mt")] to determine what chromosome name should be found in the gff file for the annotation

######################## CALCULATING FST #########################
# %% Adapted jupyter notebook from http://alimanfoo.github.io/2015/09/21/estimating-fst.html
import os
import argparse
import zarr
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style('white')
sns.set_style('ticks')
import pandas as pd
import allel
import sys
from datetime import datetime

# %%
def main(args):
    print('scikit-allel', allel.__version__)
    
    # set working directory
    os.chdir(args.working_directory)
    print("Working directory:", os.getcwd())

    # Open the callset file
    callset = zarr.open(args.callset_file, mode='r')
    print("Callset file:", args.callset_file)

    # Filter by chromosome
    chromosome_arg = args.chromosome
    chromosome_filter = callset['variants/CHROM'][:] == chromosome_arg
    pos_all = callset['variants/POS'][np.where(chromosome_filter)[0]]
    print("Chromosome being analysed:", args.chromosome)

    # Path to your GFF3 file
    gff_file_path = '/mnt/storage11/sophie/reference_genomes/A_gam_P4_ensembl/Anopheles_gambiae.AgamP4.56.chr.gff3'
    chromosome_for_gff = chromosome_arg.replace("anop_X","X")
    # Read the GFF3 file
    gff_df = pd.read_csv(
        '/mnt/storage11/sophie/reference_genomes/A_gam_P4_ensembl/Anopheles_gambiae.AgamP4.56.chr.gff3',
        sep='\t',
        comment='#',
        names=['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes'],
        usecols=['seqid', 'start', 'end', 'attributes'],
        dtype={"start": int, "end": int}
    )

    # Define the extract_gene_ID function
    def extract_gene_id(attributes_str):
        """Extract and return the gene ID from the attributes string."""
        # Split the string into a list of "key=value" strings
        attributes = attributes_str.split(';')
        # Iterate over each "key=value" string
        for attribute in attributes:
            # Split each "key=value" string into a key and value
            key, _, value = attribute.partition('=')
            # Check if the key is the one we're interested in (e.g., 'ID')
            if key == 'ID':
                return value
        # Return a placeholder or None if no ID was found
        return None
    
    # Process the 'attributes' to split it into a dictionary
    def parse_attributes(attributes_str):
        # Split the string by semicolon and then by equal sign to create a dict
        return dict(item.split('=') for item in attributes_str.split(';') if '=' in item)

    def find_annotations_for_window(window, gff_df, chromosome):
        """Find gene IDs for annotations overlapping a given window on a specific chromosome."""
        # Filter GFF for overlapping annotations on the specified chromosome
        overlapping_annotations = gff_df[
            (gff_df['seqid'] == chromosome) &
            (gff_df['start'] <= window[1]) &
            (gff_df['end'] >= window[0])
        ]
        # Extract the gene ID from each annotation's attributes
        gene_ids = overlapping_annotations['attributes'].apply(lambda x: x.get('ID', 'NA') if isinstance(x, dict) else 'NA')
        # Join all non-None gene IDs into a single string; filter out any None values first
        annotations = '; '.join(filter(None, gene_ids))
        return annotations
    
    gff_df['attributes'] = gff_df['attributes'].apply(parse_attributes)

    gff_df.head()
    
    # %%  CONVERT ZARR FILE TO GENOTYPE ARRAY
    # whole genome
    # genotype_all = allel.GenotypeChunkedArray(callset['calldata/GT'])
    # genotype_all

    # create genotype array for just chrom
    genotype_all = allel.GenotypeChunkedArray(callset['calldata/GT'][np.where(chromosome_filter)[0]])

    # check length of pos_all and genotype_all are the same
    
    if len(pos_all)==len(genotype_all):
        print("Length of positions and genotypes in the genotype array are the same, script continuing")
    else:
        print("Something is wrong with the genotype_all array as the lenght of pos_all and genotype_all are different. Stopping script.")
        sys.exit()  # This will stop the script. If you want the script to continue anyway, # out this line

    # %%  IMPORT METADATA
    df_samples= pd.read_csv('pca_investigation_metadata_melasplusglobal.csv',sep=',',usecols=['sample','country','year','species','island','pcagroup'])
    df_samples.head()
    df_samples.groupby(by=['pcagroup']).count()
    print("Imported metadata")

    # %%  choose sample populations to work with
    pop1 = 'group1'
    pop2 = 'group2'
    n_samples_pop1 = np.count_nonzero(df_samples.pcagroup == pop1)
    n_samples_pop2 = np.count_nonzero(df_samples.pcagroup == pop2)
    print("Population 1:", pop1, "Number of samples in pop1:", n_samples_pop1, "Population 2:", pop2, "Number of samples in pop2:", n_samples_pop2)

    # %% dictonary mapping population names to sample indices
    subpops = {
        pop1: df_samples[df_samples.pcagroup == pop1].index,
        pop2: df_samples[df_samples.pcagroup == pop2].index,
    }

    # %% get allele counts
    acs = genotype_all.count_alleles_subpops(subpops)
    acs

    # Check for segregating variants in each population separately
    flt_pop1_segregating = acs[pop1].is_segregating()
    flt_pop2_segregating = acs[pop2].is_segregating()

    # Ensure variants are segregating in both populations
    flt_segregating_in_both = flt_pop1_segregating & flt_pop2_segregating

    # Get combined allele counts for further biallelic filtering
    acu = allel.AlleleCountsArray(acs[pop1][:] + acs[pop2][:])

    # Ensure variants are biallelic in the union of the two populations
    flt_biallelic = (acu.max_allele() == 1)

    # Combine filters: segregating in both populations and biallelic
    flt = flt_segregating_in_both & flt_biallelic

    print('Filtered out variants that are not segregating in both populations, or are not biallelic. Now retaining', np.count_nonzero(flt), 'SNPs')

    # %% create the new genotype array

    pos = pos_all.compress(flt)
    ac1 = allel.AlleleCountsArray(acs[pop1].compress(flt, axis=0)[:, :2])
    ac2 = allel.AlleleCountsArray(acs[pop2].compress(flt, axis=0)[:, :2])
    genotype = genotype_all.compress(flt, axis=0)
    genotype

    print("Created genotype array")

    # check that pos and genotype are the same size
    if len(pos)==len(genotype):
        print("Length of positions and genotypes in the genotype array are the same, script continuing")
    else:
        print("Something is wrong with the genotype_all array as the lenght of pos_all and genotype_all are different. Stopping script.")
        sys.exit()  # This will stop the script. If you want the script to continue anyway, # out this line

    # %% PLOT FST using windowed weird and cockerham fst
    print("Plotting Fst using allel.windowed_weir_cockerham_fst, window size 1000")

    subpoplist = [list(subpops['group1']),
                    list(subpops['group2'])]

    fst, windows, counts = allel.windowed_weir_cockerham_fst(pos, genotype, subpoplist, size=1000)

    # use the per-block average Fst as the Y coordinate
    y = fst

    x = [np.mean(w) for w in windows]

    # plot
    fig, ax = plt.subplots(figsize=(10, 4))
    sns.despine(ax=ax, offset=5)
    ax.plot(x, y, 'k-', lw=.5)
    ax.set_ylabel('$F_{ST}$')
    ax.set_xlabel(f'Chromosome {args.chromosome} position (bp)')
    ax.set_xlim(0, pos.max())

    print(len(x),len(y))

    # save fst figure
    timestamp = datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
    filename = f'fst_w&c_windowed_{args.chromosome}_{pop1}_{pop2}_{timestamp}.png'
    plt.savefig(filename)
    print("Saving windowed Fst plot")

    # INSPECT Fst values
    # %% Print the maximum FST value and its corresponding window
    
    print("Inspecting windowed Fst plot for maximum value")

    max_value = max(fst)
    max_index = np.argmax(fst)
    max_window = windows[max_index]

    print("Maximum FST Value:", max_value)
    print("Corresponding Window (Genomic bps):", max_window)

    # %% Print the mean and median Fst value for the chromosome

    # Filter the array to include only non-negative (positive and zero) and non-NaN values
    filtered_array_temp = fst[(~np.isnan(fst)) & (fst > 0)]
    # Calculate the mean of the filtered array
    mean_value = np.mean(filtered_array_temp)
    print("Mean fst value:", mean_value)
    # Calculate the median of the filtered array
    median_value = np.median(filtered_array_temp)
    print("Median fst value:", median_value)

    # %% save Fst values over a certain threshold
    fst_threshold = 0.8
    fst_over_threshold = [(window, value) for window, value in zip(windows, fst) if value >= fst_threshold]

    if fst_over_threshold:
        filename = f'fst_greater_than_{fst_threshold}_{pop1}_{pop2}_{args.chromosome}_window_size_1000_v2_ann.txt'
        with open(filename, 'w') as fst_file:
            fst_file.write("Window (Genomic bps)\tFST Value\tAnnotations\n")
            for window, value in fst_over_threshold:
                annotations = find_annotations_for_window(window, gff_df, chromosome_for_gff)
                fst_file.write(f"{window[0]}-{window[1]}\t{value}\t{annotations}\n")
        print(f"Saved FST values over {fst_threshold} to {filename} with annotations")
    else:
        print(f"No FST values over the threshold of {fst_threshold} were found")

    ## %% Calculte fst with blen = 1 so that fst is calculated for every variant. Note this takes a bit longer to run.
    ## PLOT FST using windowed weird and cockerham fst - try this, not sure if x = pos
    #print("Calculating fst with blen = 1, this takes a bit longer to run as fst is being calculating for each variant instead of in windows")
#
    #subpoplist = [list(subpops['group1']),
    #                list(subpops['group2'])]
#
#
    ## fst, se, vb, _ = allel.blockwise_hudson_fst(ac1, ac2, blen=blen)
    #a, b, c = allel.weir_cockerham_fst(genotype, subpoplist, blen=1)
#
    ## estimate Fst for each variant and each allele individually
    ## fst = a / (a + b + c)
#
    ## or estimate Fst for each variant individually, averaging over alleles
    #print("Estimating Fst for each variant individually, averaging over alleles")
#
    #fst_pervariant = (np.sum(a, axis=1) /
    #    (np.sum(a, axis=1) + np.sum(b, axis=1) + np.sum(c, axis=1)))
#
    ## use the per-block average Fst as the Y coordinate
    #print("Plotting fst with blen = 1")
#
    #y = fst_pervariant
#
    #x = pos
#
    ## plot
    #fig, ax = plt.subplots(figsize=(10, 4))
    #sns.despine(ax=ax, offset=5)
    #ax.plot(x, y, 'k-', lw=.5)
    #ax.set_ylabel('$F_{ST}$')
    #ax.set_xlabel(f'Chromosome {args.chromosome} position (bp)')
    #ax.set_xlim(0, pos.max())
#
    #print(len(x),len(y))
#
    ## save fst figure
    #timestamp = datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
    #filename = f'fst_w&c_blen=1_{pop1}_{pop2}_{timestamp}_{args.chromosome}_v2.png'
    #plt.savefig(filename)
    #print("Saving blen = 1 Fst plot - {args.chromosome}")
#
    ## Here you can see that fst values are much higher when calculated for each position individually
#
    #print("Inspecting Fst plot for maximum value")
    #max_value = max(fst_pervariant)
    #max_index = np.argmax(fst_pervariant)
    ##max_window = windows[max_index]
#
    #print("Maximum FST Value:", max_value)
    #print("Maximum index:", max_index)
#
    ## Filter FST values greater than or equal to 0.6
    #fst_threshold = 0.6
    #fst_over_threshold = [(pos[i], value) for i, value in enumerate(fst_pervariant) if value >= fst_threshold]
#
    #if fst_over_threshold:
    #    # Dynamically generate filename based on fst_threshold
    #    filename = f'fst_individual_variants_greater_than_{fst_threshold}_{args.chromosome}_blen1_v2_ann.txt'
    #    with open(filename, 'w') as fst_file:
    #        fst_file.write("Variant (Genomic bps)\tFST Value\n")
    #        for variant, value in fst_over_threshold:
    #            annotations = find_annotations_for_window(window, gff_df, chromosome_for_gff)
    #            fst_file.write(f"{window[0]}-{window[1]}\t{value}\t{annotations}\n")
    #    # Dynamically include fst_threshold in the print statement
    #    print(f"Saved individual variant FST values over {fst_threshold} to {filename} with annotations")
    #else:
    #    # Dynamically include fst_threshold in the print statement
    #    print(f"No individual variant FST values over the threshold (FST >= {fst_threshold}) were found.")
#
# arguments

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculate and plot FST using scikit-allel.")
    parser.add_argument("working_directory", type=str, help="Working directory for the script.")
    parser.add_argument("callset_file", type=str, help="Path to the callset file.")
    parser.add_argument("chromosome", type=str, help="Chromosome to analyze.")
    parser.set_defaults(func=main)
    args = parser.parse_args()
    args.func(args)
