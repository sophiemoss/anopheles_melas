## scikitallel_workflow
# /mnt/storage11/sophie/gitrepos/anopheles_melas/melas_fst_v1.py
# run this script using
# python /mnt/storage11/sophie/gitrepos/anopheles_melas/melas_fst_v1.py /pathway/to/dir /pathway/to/zarr chromosomename
# python /mnt/storage11/sophie/gitrepos/anopheles_melas/melas_fst_v2.py . 2019melasglobal_finalfiltered_gambiaealigned_phased.zarr 2L
# You make the zarr file with allel.vcf_to_zarr('phased_vcf_file_name.vcf.gz', 'output_name.zarr', fields='*', overwrite=True)

# allel.vcf_to_zarr('2019_melas_phased.vcf.gz', '2019_melas_phased.zarr', fields='*', overwrite=True)
# python /mnt/storage11/sophie/gitrepos/anophelesmelas_popgen/scikit_allel_fst_calculate_and_plot_melasPCAinvestigation.py /mnt/storage11/sophie/bijagos_mosq_wgs/2019_melas_fq2vcf/genomics_database_melas2019plusglobal/genomics_database_melas2019plusglobal_vcf/melas2019plusglobal_vcf_filtering 2019_melas_phased.zarr 2L

######################## CALCULATING FST #########################
# %% adapted jupyter notebook from http://alimanfoo.github.io/2015/09/21/estimating-fst.html

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
    chromosome_filter = callset['variants/CHROM'][:] == args.chromosome
    pos_all = callset['variants/POS'][np.where(chromosome_filter)[0]]
    print("Chromosome being analysed:", args.chromosome)

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
    df_samples= pd.read_csv('metadata_melasplusglobal.csv',sep=',',usecols=['sample','country','year','species','island','pcagroup'])
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
    print("Plotting Fst using allel.windowed_weir_cockerham_fst, window size 500")

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

    # %% save Fst values over a certain threshold
    fst_threshold = 0.9
    fst_over_threshold = [(window, value) for window, value in zip(windows, fst) if value >= fst_threshold]

    if fst_over_threshold:
    # Dynamically include fst_threshold in the filename
        filename = f'fst_greater_than_{fst_threshold}_{pop1}_{pop2}_{args.chromosome}_window_size_1000_v2.txt'
        with open(filename, 'w') as fst_file:
            fst_file.write("Window (Genomic bps)\tFST Value\n")
            for window, value in fst_over_threshold:
                fst_file.write(f"{window[0]}-{window[1]}\t{value}\n")
        print(f"Saved FST values over {fst_threshold} to {filename}")
    else:
        print(f"No FST values over the threshold of {fst_threshold} were found")

    # %% Calculte fst with blen = 1 so that fst is calculated for every variant. Note this takes a bit longer to run.
    # PLOT FST using windowed weird and cockerham fst - try this, not sure if x = pos
    print("Calculating fst with blen = 1, this takes a bit longer to run as fst is being calculating for each variant instead of in windows")

    subpoplist = [list(subpops['group1']),
                    list(subpops['group2'])]


    # fst, se, vb, _ = allel.blockwise_hudson_fst(ac1, ac2, blen=blen)
    a, b, c = allel.weir_cockerham_fst(genotype, subpoplist, blen=1)

    # estimate Fst for each variant and each allele individually
    # fst = a / (a + b + c)

    # or estimate Fst for each variant individually, averaging over alleles
    print("Estimating Fst for each variant individually, averaging over alleles")

    fst_pervariant = (np.sum(a, axis=1) /
        (np.sum(a, axis=1) + np.sum(b, axis=1) + np.sum(c, axis=1)))

    # use the per-block average Fst as the Y coordinate
    print("Plotting fst with blen = 1")

    y = fst_pervariant

    x = pos

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
    filename = f'fst_w&c_blen=1_{pop1}_{pop2}_{timestamp}_{args.chromosome}_v2.png'
    plt.savefig(filename)
    print("Saving blen = 1 Fst plot - {args.chromosome}")

    # Here you can see that fst values are much higher when calculated for each position individually

    print("Inspecting Fst plot for maximum value")
    max_value = max(fst_pervariant)
    max_index = np.argmax(fst_pervariant)
    #max_window = windows[max_index]

    print("Maximum FST Value:", max_value)
    print("Maximum index:", max_index)

    # Filter FST values greater than or equal to 0.9
    fst_threshold = 0.9
    fst_over_threshold = [(pos[i], value) for i, value in enumerate(fst_pervariant) if value >= fst_threshold]

    if fst_over_threshold:
        # Dynamically generate filename based on fst_threshold
        filename = f'fst_individual_variants_greater_than_{fst_threshold}_blen1_v2.txt'
        with open(filename, 'w') as fst_file:
            fst_file.write("Variant (Genomic bps)\tFST Value\n")
            for variant, value in fst_over_threshold:
                fst_file.write(f"{variant}\t{value}\n")
        # Dynamically include fst_threshold in the print statement
        print(f"Saved individual variant FST values over {fst_threshold} to {filename}")
    else:
        # Dynamically include fst_threshold in the print statement
        print(f"No individual variant FST values over the threshold (FST >= {fst_threshold}) were found.")

# arguments

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculate and plot FST using scikit-allel.")
    parser.add_argument("working_directory", type=str, help="Working directory for the script.")
    parser.add_argument("callset_file", type=str, help="Path to the callset file.")
    parser.add_argument("chromosome", type=str, help="Chromosome to analyze.")
    parser.set_defaults(func=main)
    args = parser.parse_args()
    args.func(args)

    ## A. Miles used Hudson Fst estimator