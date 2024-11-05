'''
This file contains general functions, in use for this project.
'''
import os

import pandas as pd
from Bio import SeqIO
from Bio import Phylo
import dendropy
import sys
from warnings import warn
from statsmodels.stats.proportion import proportion_confint
from datetime import datetime, timedelta
from dateutil import parser
from ete3 import Tree
from matplotlib import pyplot as plt


# Colors list

DEFAULT_COLORS = plt.rcParams['axes.prop_cycle'].by_key()['color']
DARK_GREEN = '#009E73'
LIGHT_GREEN = '#F0E442'
GRAY = '#BFBFBF'
WHITE = '#FFFFFF'

# Folders

ANALYSES = 'analyses'
DATA = 'data'
FIGURES = 'figures'

GENOMICS = 'genomics'
NETWORK = 'network'

ARTICLE_FIGURES = 'article_figures_intermediate_files'
GENOMIC_VALIDATION = 'genomic_validation'
PASTML = 'pastml'
SNP = 'snp'
TEMPEST = 'tempest'
THORNEY_BEAST = 'thorney_beast'
TOPOLOGY_INVESTIGATION = 'topology_investigation'
TREES = 'trees'
TREETIME = 'treetime'

BA1 = 'ba1'
BA2 = 'ba2'
BA1_POLYTOMIES = 'ba1_polytomies'
BA2_POLYTOMIES = 'ba2_polytomies'
BA1_THORNEYBEAST = 'ba1_tb'
BA2_THORNEYBEAST = 'ba2_tb'

ZERO_SNP_NONCLUSTERED = '0snp_non_clustered'
ONE_SNP_NONCLUSTERED = '1snp_non_clustered'
TWO_SNP_NONCLUSTERED = '2snp_non_clustered'
THREE_FIVE_SNP_CLUSTERED_TB = '3-5snp_clustered_tb'
PAIR15_INVESTIGATION = 'pair15_investigation'
PAIR33_INVESTIGATION = 'pair33_investigation'


'''
@df: pandas DataFrame.
@columns: columns based on which duplicated values are removed.
@keep: Which duplicated occurrence do we keep?
@skipna [Boolean]: ignores NA values in @na_check_col when removing duplicates
@na_check_col [String]: if @skipna is True, and the value in @na_check_col is None, do not remove that line.
This function removes duplicates in a data frame.
'''
def remove_duplicated(df, columns, keep='first', skipna=False, na_check_col=None):
    assert not (skipna and na_check_col is None), 'Yous must specify on which column to check na values'

    if skipna:
        is_duplicated = df.duplicated(subset=columns, keep=keep) & ~df[na_check_col].isna()
    else:
        is_duplicated = df.duplicated(subset=columns, keep=keep)

    df = df[~is_duplicated]

    return df


'''
@metadata_path: Reference TSV file.
@sequences_path_in: FASTA file that includes all sequences from metadata, and more.
@sequences_path_out: FASTA file in which all sequences from metadata will be writen. 
This function outputs within sequences_path_out, within a FASTA format, each sequence included within metadata.
Those sequences must be included within sequences_path_in (FASTA).
'''
def keep_sequences_from_metadata(metadata_path, sequences_path_in, sequences_path_out):
    df = pd.read_csv(metadata_path, delimiter='\t')

    set_to_keep = set(df['name'].values)

    ffile = SeqIO.parse(sequences_path_in, "fasta")

    with open(sequences_path_out, 'w') as f:
        for seq_record in ffile:
            try:
                set_to_keep.remove(seq_record.name)
            except KeyError:
                # The sequence is not included within set_to_keep. Just skip it.
                continue
            else:
                f.write(seq_record.format("fasta"))

        if len(set_to_keep) != 0:
            print(f'{len(set_to_keep)} of the headers from list were not identified in the input fasta file.')

    cmd = 'grep -c ">" {0:s}'.format(sequences_path_out)
    nseq = int(os.popen(cmd).read())
    nrow = len(df.index)

    assert (nseq == nrow), f'The number of rows {(nrow)} in the METADATA file does not correspond to the number of sequences {(nseq)} in the FASTA file'

    print(f'We kept {nseq} sequences')


'''
@sequences_to_exclude: CSV files with all tips to exclude.
@tree_in: Newick tree with outliers to remove.
@tree_out: Newick tree after removing outliers.
@column_name: Name of the column in CSV file that contains sequences to remove.
This function prunes the tree by excluding all sequences listed in @sequences_to_exclude.
'''
def drop_tree_sequences(sequences_to_exclude, tree_in, tree_out, column_name='outliers'):
    df = pd.read_csv(sequences_to_exclude)

    tree0 = dendropy.Tree.get(
        path=tree_in,
        schema="newick",
        preserve_underscores=True
    )

    num_tips0 = len(tree0.leaf_nodes())

    labels = set(df[column_name].values)
    tree1 = tree0.extract_tree_without_taxa_labels(labels=labels)

    num_tips1 = len(tree1.leaf_nodes())

    assert num_tips1 == num_tips0 - len(labels), "The number of excluded sequences does not match the number of sequences specified in @sequences_to_exclude file."

    tree1.write(path=tree_out, schema='newick', unquoted_underscores=True)


'''
@sequences_to_exclude: CSV files with all tips to exclude.
@metadata_in: Metadata file with outliers to remove.
@metadata_out: Metadata file after removing outliers.
@column_name: Name of the column in CSV file that contains sequences to remove.
This function excludes all sequences listed in @sequences_to_exclude from @metadata_in.
'''
def drop_metadata_sequences(sequences_to_exclude, metadata_in, metadata_out, column_name='outliers'):
    sequences = pd.read_csv(sequences_to_exclude).rename(columns={column_name: 'name'})
    sequences.insert(loc=1, column='to_exclude', value=True)

    metadata = pd.read_csv(metadata_in, sep='\t')
    metadata = pd.merge(metadata, sequences, how='left')

    assert len(metadata[metadata['to_exclude'] == True].index) == len(sequences.index), "The number of excluded sequences does not match the number of sequences specified in @sequences_to_exclude file."

    metadata = metadata[metadata['to_exclude'].isna()].drop('to_exclude', axis=1)

    metadata.to_csv(metadata_out, sep='\t', index=False)


'''
@sequences_to_exclude: CSV files with all tips to exclude.
@fasta_in: FASTA file with outliers to remove.
@fasta_out: FASTA file after removing outliers.
@column_name: Name of the column in CSV file that contains sequences to remove.
This function excludes all sequences listed in @sequences_to_exclude from @fasta_in.
'''
def drop_fasta_sequences(sequences_to_exclude, fasta_in, fasta_out, column_name='outliers'):
    exclusion_sequences = pd.read_csv(sequences_to_exclude)[column_name].values
    seqiter = SeqIO.parse(fasta_in, 'fasta')

    sys.stdout = open(fasta_out, 'w')
    seq_to_write = list()
    n_written_seq = 0
    n_seq = 0
    for seq in seqiter:
        n_seq += 1
        if seq.id not in exclusion_sequences:
            n_written_seq += 1
            seq_to_write.append(seq)

    assert n_written_seq == n_seq - len(exclusion_sequences), "The number of excluded sequences does not match the number of sequences specified in @sequences_to_exclude file."

    SeqIO.write(seq_to_write, sys.stdout, "fasta")

    sys.stdout.close()


'''
@sequences_to_keep: TSV file with all sequences to keep, writen within the @column_name column.
@fasta_in: FASTA file.
@fasta_out: FASTA file with only sequences writen in @sequences_to_keep file.
@column_name: Name of the column in TSV file that contains sequences to keep.
This function extracts the sequences included within @sequences_to_keep from the FASTA @fasta_in and outputs them into a new FASTA file @fasta_out.
'''
def fasta_subset(sequences_to_keep, fasta_in, fasta_out, column_name='node'):
    wanted_sequences = list(pd.read_csv(sequences_to_keep, delimiter='\t')[column_name].values)

    seqiter = SeqIO.parse(fasta_in, 'fasta')
    sys.stdout = open(fasta_out, 'w')
    seq_to_write = [None] * len(wanted_sequences)
    n_wanted_seq = len(wanted_sequences)

    n_written_seq = 0
    for seq in seqiter:
        if seq.id in wanted_sequences:
            idx = wanted_sequences.index(seq.id)
            n_written_seq += 1
            seq_to_write[idx] = seq

    not_found_indexes = [index for index, value in enumerate(seq_to_write) if value is None]
    not_found_sequences = [wanted_sequences[i] for i in not_found_indexes]

    assert n_written_seq == n_wanted_seq, f"We did not find the {not_found_sequences} wanted sequences within the provided FASTA file."

    SeqIO.write(seq_to_write, sys.stdout, "fasta")

    sys.stdout.close()


'''
@sequences_to_keep: TSV file with all sequences to keep, writen within the @column_name column.
@metadata_in: metadata TSV file.
@metadata_out: metadata TSV file only containing sequences to keep.
@column_name: Name of the column in @sequences_to_keep file that contains sequences to keep.
This function extracts the sequences included within @sequences_to_keep from the FASTA @fasta_in and outputs them into a new FASTA file @fasta_out.
'''
def metadata_subset(sequences_to_keep, metadata_in, metadata_out, column_name='node'):
    wanted_sequences = pd.read_csv(sequences_to_keep, delimiter='\t').rename(columns={column_name: 'name'})
    metadata = pd.read_csv(metadata_in, delimiter='\t')

    metadata = pd.merge(metadata, wanted_sequences, how='inner')

    metadata.to_csv(metadata_out, index=False, sep='\t')


'''
@list_sequences_to_keep: list object including all sequences to keep.
@fasta_in: FASTA file.
@fasta_out: FASTA file with only sequences writen in @sequences_to_keep file.
This function extracts the sequences included within @list_sequences_to_keep from the FASTA @fasta_in and outputs them into a new FASTA file @fasta_out.
'''
def fasta_subset_from_list(list_sequences_to_keep, fasta_in, fasta_out):
    seqiter = SeqIO.parse(fasta_in, 'fasta')
    sys.stdout = open(fasta_out, 'w')
    seq_to_write = [None] * len(list_sequences_to_keep)
    n_wanted_seq = len(list_sequences_to_keep)

    n_written_seq = 0
    for seq in seqiter:
        if seq.id in list_sequences_to_keep:
            idx = list_sequences_to_keep.index(seq.id)
            n_written_seq += 1
            seq_to_write[idx] = seq

    not_found_indexes = [index for index, value in enumerate(seq_to_write) if value is None]
    not_found_sequences = [list_sequences_to_keep[i] for i in not_found_indexes]

    assert n_written_seq == n_wanted_seq, f"We did not find the {not_found_sequences} wanted sequences within the provided FASTA file."

    SeqIO.write(seq_to_write, sys.stdout, "fasta")

    sys.stdout.close()


'''
@treetime_std_out: Treetime's stdout file.
@excluded_sequences: CSV file in which that contains all Treetime outliers.
@white_list_pattern: All sequences which name contains @white_list_pattern won't be included within @excluded_sequences file
@column_name: Name of the column that contains all excluded sequences name, within @excluded_sequences file.
This function reads Treetime output file and creates the @excluded_sequences file that contains sequences considered as outliers by Treetime.
If treetime defines as an outlier a sequence which name contains @white_list_pattern, do not write it in the @excluded_sequences file, and gives a warning.
'''
def output_treetime_outliers(treetime_std_out, excluded_sequences, white_list_pattern='rega', column_name='outliers'):
    sequences_list = []

    with open(treetime_std_out) as fp:
        for l in fp.readlines():
            if "input date" in l:
                seq = l.split(' ')[0].split('\t')[1][:-1]
                if white_list_pattern in seq:
                    warn(f'WARNING: TreeTime set {seq} sequence as an outlier.')
                else:
                    sequences_list.append(seq)

    df = pd.DataFrame(data={column_name: sequences_list})
    df.to_csv(excluded_sequences, index=False)


'''
@nexus_tree: Input tree in NEXUS format.
@newick_tree: Output tree in NEWICK format.
This function converts a NEXUS tree into a NEWICK tree.
'''
def convert_tree(nexus_tree, newick_tree):
    tree = dendropy.Tree.get(
        path=nexus_tree,
        schema="nexus",
        preserve_underscores=True
    )
    tree.write(path=newick_tree, schema='newick', suppress_rooting=True)


'''
@tree: Tree in NEWICK format.
This function returns the number of internal nodes within a tree.
'''
def count_internal_nodes(tree):
    tree = Phylo.read(tree, "newick")
    internal_node_count = 0
    for clade in tree.find_clades():
        if not clade.is_terminal():
            internal_node_count += 1
    return internal_node_count


'''
@tree_in: Input tree in newick format.
@tree_out: Tree with renamed internal nodes, in newick format.
This function renames internal nodes within @tree_in and outputs the result into @tree_out.
'''
def rename_internal_node(tree_in, tree_out):
    tree = Phylo.read(tree_in, "newick")
    n_internal_nodes = count_internal_nodes(tree_in)
    n_digits = len(str(n_internal_nodes))
    node_count = 0
    for clade in tree.find_clades(order="postorder"):
        if not clade.is_terminal():
            clade.name = f"NODE_{node_count:0{n_digits}d}"
            node_count += 1

    Phylo.write(tree, tree_out, "newick")


'''
@tree: Tree in newick format.
@tip1: Tip label 1 within the @tree.
@tip2: Tip label 2 within the @tree.
This function outputs the name of the most recent common ancestor between @tip1 and @tip2.
'''
def get_MRCA(tree, tip1, tip2):
    tree = Phylo.read(tree, "newick")
    mrca = tree.common_ancestor([tip1, tip2])

    return mrca.name


'''
@tree_in: [PATH] NEWICK input tree, with 0-length branches.
@tree_out: [PATH] NEWICK output tree, without 0-length branches.
This function replaces 0-length branches with 0.00001 branches and writes the result into a new file.
'''
def remove_zero_branches(tree_in, tree_out):
    with open(tree_in, 'r') as file:
        tree_string = file.read()

    corrected_tree_string = tree_string.replace('-0.0', '0.00001')

    with open(tree_out, 'w') as file:
        file.write(corrected_tree_string)


'''
@n_num: [DOUBLE] Fraction numerator.
@n_denom: [DOUBLE] Fraction denominator.
@meth: [STRING] Method used to compute the confidence interval.
This function computes a 95% confidence interval for a proportion, expressed in [%].
'''
def ci_ratio_num_denom(n_num, n_denom, meth='wilson'):
    if n_denom == 0:
        ratio = 0,
        ci_low = 0
        ci_upp = 0
    else:
        ratio = 100 * n_num / n_denom
        ci_low, ci_upp = proportion_confint(count=n_num, nobs=n_denom, alpha=0.05, method=meth)
        ci_low *= 100
        ci_upp *= 100

    return ratio, ci_low, ci_upp


'''
@date_string: [STRING] Date represented as a string in format "YYY-MM-DD".
This function returns the decimal date, computed from a date written as a string.
'''
def decimal_date(date_string):
    date = parser.parse(date_string)

    year_start = datetime(date.year, 1, 1)
    next_year_start = datetime(date.year + 1, 1, 1)

    days_since_year_start = (date - year_start).days
    total_days_in_year = (next_year_start - year_start).days
    decimal_year = date.year + days_since_year_start / total_days_in_year

    return decimal_year


'''
@tree: [PATH] Timetree in newick format
@metadatafile: [PATH] The corresponding metadata file with all leaves dates information.
This function returns the inferred time of the root of @tree in decimal date format.
'''
def root_time(tree, metadatafile):
    metadata = pd.read_csv(metadatafile, delimiter='\t')
    tree = Phylo.read(tree, 'newick')

    root = tree.root
    leaf = tree.get_terminals()[0]
    leaf_sampling_date = metadata.loc[metadata['name'] == leaf.name, 'date'].values[0]
    leaf_decimal_sampling_date = decimal_date(leaf_sampling_date)
    root_age =  leaf_decimal_sampling_date - tree.distance(leaf, root)

    return root_age


'''
@file_path: [PATH] The file to check.
This function checks whether a the file @file_path exists, and raises an error if it does.
'''
def check_file(file_path):
    if os.path.exists(file_path):
        raise FileExistsError(f"Error: The file '{file_path}' already exists.")


'''
@tree_in: [PATH] Input NEWICK tree file, containing @node.
@tree_out: [PATH] Outputs Newick tree file only containing subtree descending from @node.
@node: [String] Root node of the tree to be extracted.
@format: [Integer] Tree format according to ete3.Tree() object definition.
This function extracts the clade descending from @node within @tree_in file and outputs the result within @tree_out file.
'''
def extract_subtree(tree_in, tree_out, node, format=1):
    t = Tree(tree_in, format=format, quoted_node_names=True)

    for n in t.traverse():
        if n.name == node:
            n.write(format=format, outfile=tree_out)


def get_mrca_time_old(newick_file, metadata_file, tip1, tip2):
    """
    Returns the time of a node in a time-scaled phylogenetic tree parsed from a Newick file.

    :param newick_file: Path to the Newick file containing the phylogenetic tree.
    :param target_node_name: The name of the node to compute the time for (tip or internal node).
    :param metadata: A dictionary mapping tip names to their sampling dates (in "YYYY-MM-DD" format).

    :return: The computed time (in years) of the node relative to the present.
    """
    target_node_name = get_MRCA(newick_file, tip1, tip2)
    metadata = pd.read_csv(metadata_file, delimiter='\t')

    # Load the tree from the Newick file
    tree = Phylo.read(newick_file, "newick")

    # Convert metadata sampling dates to datetime objects
    date_format = "%Y-%m-%d"
    sampling_dates = {row['name']: datetime.strptime(row['date'], date_format) for _, row in metadata.iterrows()}
    # Find the most recent sampling date to consider as the "present" time reference
    most_recent_date = max(sampling_dates.values())

    # Recursive function to calculate the time of a node relative to the most recent sampling date
    def compute_time_from_root(node):
        if node.is_terminal():  # leaf node
            # Get the sampling date of the leaf
            leaf_name = node.name
            if leaf_name in sampling_dates:
                sampling_date = sampling_dates[leaf_name]
                # Calculate time difference (in years) from the most recent date
                time_diff_years = (most_recent_date - sampling_date).days / 365.25
                return time_diff_years
            else:
                raise ValueError(f"Sampling date not found for tip: {leaf_name}")
        else:
            # Internal node: calculate time recursively based on child node distances
            child_times = [compute_time_from_root(child) for child in node.clades]
            # Add the branch length to the root (which is in years)
            return min(child_times) + node.branch_length


    # Find the target node by name
    for clade in tree.find_clades():
        if clade.name == target_node_name:
            # Calculate and return the time for the requested node
            time_in_years = compute_time_from_root(clade)
            print(time_in_years)
            # Convert time (years ago) into a date by subtracting from the most recent date
            divergence_date = most_recent_date - timedelta(days=time_in_years * 365.25)
            return divergence_date

    raise ValueError(f"Node with name {target_node_name} not found in the tree.")


def get_mrca_time(newick_file, metadata_file, tip1, tip2):
    """
    Returns the time of a node in a time-scaled phylogenetic tree parsed from a Newick file.

    :param newick_file: Path to the Newick file containing the phylogenetic tree.
    :param target_node_name: The name of the node to compute the time for (tip or internal node).
    :param metadata: A dictionary mapping tip names to their sampling dates (in "YYYY-MM-DD" format).

    :return: The computed time (in years) of the node relative to the present.
    """
    target_node_name = get_MRCA(newick_file, tip1, tip2)
    metadata = pd.read_csv(metadata_file, delimiter='\t')

    # Load the tree from the Newick file
    tree = Phylo.read(newick_file, "newick")

    # Convert metadata sampling dates to datetime objects
    date_format = "%Y-%m-%d"
    sampling_dates = {row['name']: datetime.strptime(row['date'], date_format) for _, row in metadata.iterrows()}
    # Find the most recent sampling date to consider as the "present" time reference

    # Recursive function to calculate the time of a node relative to the most recent sampling date
    def compute_time_to_leaf(node):
        if node.is_terminal():  # leaf node
            return node.branch_length, node.name
        else:
            child_branch_length, found_leaf = compute_time_to_leaf(node.clades[0])
            return  child_branch_length + node.branch_length, found_leaf


    # Find the target node by name
    for clade in tree.find_clades():
        if clade.name == target_node_name:
            # Calculate and return the time for the requested node
            time_in_years, reference_leaf = compute_time_to_leaf(clade[0])
            reference_date = sampling_dates[reference_leaf]
            # Convert time (years ago) into a date by subtracting from the most recent date
            divergence_date = reference_date - timedelta(days=time_in_years * 365.25)
            return divergence_date

    raise ValueError(f"Node with name {target_node_name} not found in the tree.")
