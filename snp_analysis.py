'''
This file contains the functions to run the SNP analysis, above the phylogenetic analysis.
'''
from commons import *

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


'''
@pairs_in: The initial pairs file
@pairs_out: Pairs file with SNP distances
@distances: SNP distance matrix, output from snp-dists software.
This function adds the SNP distances in a new column 'd' for each pair in @pairs_in, based on @distances file.
It outputs the result into @pairs_out.
'''
def add_snp_to_pairs(pairs_in, pairs_out, distances):
    pairs = pd.read_csv(pairs_in)

    distmat = pd.read_csv(distances, delimiter='\t', index_col=0)

    pairs.insert(loc=len(pairs.columns), column='d', value=None)

    for i in pairs.index:
        seq1 = pairs.loc[i, 'strain.x']
        seq2 = pairs.loc[i, 'strain.y']

        d = distmat.loc[distmat.index == seq1, seq2]

        pairs.loc[i, 'd'] = d[0]

    pairs.to_csv(pairs_out, index=False)


'''
@ba1_pairs: case-contact pairs CSV file where both patients are infected with Omicron BA.1.
@ba2_pairs: case-contact pairs CSV file where both patients are infected with Omicron BA.2.
@ba1ba2_pairs: case-contact pairs CSV file where patients are infected with different variants.
@all_pairs_phylo: all case-contact pairs merged together in a CSV file, with phylogenetic information.
This function merges all BA.1, BA.2, and distinct variants data sets together and saves the result into a new data frame.
'''
def merge_all_pairs(ba1_pairs, ba2_pairs, ba1ba2_pairs, all_pairs_phylo):
    ba1_pairs = pd.read_csv(ba1_pairs)
    ba2_pairs = pd.read_csv(ba2_pairs)
    ba1ba2_pairs = pd.read_csv(ba1ba2_pairs)

    ba1ba2_pairs.insert(loc=len(ba1ba2_pairs.columns), column="phylogeny_validation", value=False)

    df = pd.merge(ba1_pairs, ba2_pairs, how='outer')
    df = pd.merge(df, ba1ba2_pairs, how='outer')

    assert not any(df.duplicated(subset=['CV.x', 'CV.y'])), "Some duplicate pairs were found within the whole data set."
    assert not any((df['variant.x'] != df['variant.y']) & (df['phylogeny_validation'] == True)), "Some pairs infected with different variants were marked as phylogenetically validated."
    assert len(df.index) == len(ba1_pairs.index) + len(ba2_pairs.index) + len(ba1ba2_pairs.index), "Some pairs were duplicated."

    df.to_csv(all_pairs_phylo, index=False)

    print('')


'''
@ct_pairs: case-contact pairs file with 'phylogeny_validation' (based on clustering) and 'd' (SNP distances) columns.
@out_graph: coloured distribution of SNP distances.
This function plots the distribution of SNP distances among all case-contact pairs.
The bars are coloured according to phylogenetic clustering.
'''
def clustering_on_snp_graph(ct_pairs, out_graph, colors=[DARK_GREEN, LIGHT_GREEN], show_snp_correction=True, plot_legend=True):
    pairs = pd.read_csv(ct_pairs)

    pairs = pairs[pairs['variant.x'] == pairs['variant.y']]

    v = pairs['d'].values

    phylo_ct_dist = pairs.loc[pairs['phylogeny_validation'] == True, 'd'].values
    non_phylo_ct_zero_dist = pairs.loc[(pairs['phylogeny_validation'] == False) & (pairs['d'] == 0), 'd'].values
    non_phylo_ct_non_zero_dist = pairs.loc[(pairs['phylogeny_validation'] == False) & (pairs['d'] != 0), 'd'].values

    bin_edges = np.arange(v.min(), v.max() + 2)

    n_phylo, _ = np.histogram(phylo_ct_dist, bins=bin_edges)
    n_non_phylo_non_zero, _ = np.histogram(non_phylo_ct_non_zero_dist, bins=bin_edges)
    n_non_phylo_zero, _ = np.histogram(non_phylo_ct_zero_dist, bins=bin_edges)

    fig, ax = plt.subplots(1)
    fig.set_figheight(6)
    fig.set_figwidth(8.5)
    plt.rcParams['legend.title_fontsize'] = 18

    plt.title("Distribution of [same lineage] case-contact pairs distances", fontsize=18)
    ax.set_xlabel("SNP distances", fontsize=18)
    ax.set_ylabel("Number of pairs", fontsize=18)

    bin_width = np.diff(bin_edges)[0]

    bar_width = bin_width * 0.8

    ax.bar(bin_edges[:-1], n_phylo, width=bar_width, align='edge', color=colors[0],
           label='Same phylogenetic cluster')

    cumulative_heights = n_phylo

    ax.bar(bin_edges[:-1], n_non_phylo_non_zero, width=bar_width, align='edge',
           bottom=cumulative_heights, color=colors[1], label='Different phylogenetic clusters')

    cumulative_heights += n_non_phylo_non_zero

    if show_snp_correction:
        bars_non_phylo_zero = ax.bar(bin_edges[:-1], n_non_phylo_zero, width=bar_width, align='edge',
                                     bottom=cumulative_heights, color='none', label='Non-clustered identical sequences')
    else:
        bars_non_phylo_zero = ax.bar(bin_edges[:-1], n_non_phylo_zero, width=bar_width, align='edge',
                                     bottom=cumulative_heights, color='none')

    bars_non_phylo_zero[0].set_facecolor(colors[1])
    if show_snp_correction:
        bars_non_phylo_zero[0].set_edgecolor(colors[0])
        bars_non_phylo_zero[0].set_hatch('///')

    if plot_legend:
        ax.legend(title='Phylogenetic clustering classification', fontsize=18)

    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    ax.set_xticks(bin_centers)
    ax.set_xticklabels([int(x) for x in bin_centers], rotation=0, fontsize=14)
    ax.tick_params(axis='y', labelsize=14)  # Change y-tick font size

    plt.tight_layout()

    plt.savefig(out_graph, dpi=300)


'''
This function computes precision based on the pairs_phylo file.
This file must have a 'phylogeny_validation' column.
This function is for instance useful to compute precision from the output of the ace analysis, which source code
is stored within ml_dta.R.
'''
def ct_precision(pairs_clusters_snp):
    pairs = pd.read_csv(pairs_clusters_snp)

    n_correct_pairs = len(pairs[(pairs['phylogeny_validation'] == True) | (pairs['d'] == 0)].index)
    n_tot = len(pairs.index)

    # Agresti-Coull Method (according to ChatGPT)
    # Description: This method adjusts the Wald interval to improve coverage accuracy, especially for smaller sample sizes.
    # Use Case: Suitable for small to moderate sample sizes and when a more accurate interval is desired compared to the Wald method.
    # Pros: Provides better coverage properties than the Wald interval.
    # Cons: Slightly more complex calculation than the Wald method.
    prec, ci_low, ci_upp = ci_ratio_num_denom(n_correct_pairs, n_tot, meth='agresti_coull')

    print(f'CT precision: {n_correct_pairs}/{n_tot} = {prec:.2f}% [{ci_low:.2f}%, {ci_upp:.2f}%]')

    return prec, ci_low, ci_upp
