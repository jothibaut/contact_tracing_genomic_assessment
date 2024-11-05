'''
This file contains the functions used to run the DTA, with PastML Python package.
'''
from pastml.acr import pastml_pipeline
from commons import *
import time
import pandas as pd


'''
@metadata: Must only contains tips in @tree
@traits: TSV file containing all traits.
'''
def run_pastml(metadata, traits, all_tips_traits, tree, html_compressed, html, prediction_method='MPPA'):
    assert prediction_method == 'MPPA' or prediction_method == 'MAP', f"The {prediction_method} prediction method is not supported."

    # We start from tempest metadata file since some tips were excluded. We will merge traits to this file
    metadata = pd.read_csv(metadata, delimiter='\t')[['name']].rename(columns={'name': 'node'})
    traits = pd.read_csv(traits, delimiter='\t')[['node', 'component']]

    metadata = pd.merge(metadata, traits, how='left')

    assert not any(metadata.duplicated(['node'])), "Some metadata lines got duplicated during the merging process."
    assert len(metadata.loc[~metadata['component'].isna()].index) == len(traits.index), "Not all remaining sequences with an associated trait were found in the metadata."

    metadata.loc[metadata['component'].isna(), 'component'] = 'bg'
    metadata['component'] = metadata['component'].astype(str)
    metadata['component'] = metadata['component'].apply(lambda x: x.split('.')[0] if '.' in x else x)
    metadata.to_csv(all_tips_traits, sep='\t', index=False)

    start_time = time.time()

    # By default, this function uses the MPPA prediction method.
    # It sometimes keeps multiple state predictions per node but only when they have similar and high probabilities.
    # See https://pastml.pasteur.fr/help (MPPA section).
    # If we decided to only keep the most likely state for each internal node, we would specify prediction_method='MAP'
    pastml_pipeline(data=all_tips_traits,
                    data_sep='\t',
                    columns=['component'],
                    name_column='component',
                    tree=tree,
                    html_compressed=html_compressed,
                    html=html,
                    verbose=True,
                    prediction_method=prediction_method)

    end_time = time.time()
    elapsed_time = end_time - start_time
    elapsed_time_hour = elapsed_time / 3600
    print(f"Elapsed time: {elapsed_time_hour:.5f} hours")


def get_pastml_clusters_mppa(traits, ancestral_states, pairs, tree, pairs_phylo, strain_col='strain'):
    traits = pd.read_csv(traits, delimiter='\t')
    traits['component'] = traits['component'].astype(str)
    ancestral_states = pd.read_csv(ancestral_states, delimiter='\t')
    pairs = pd.read_csv(pairs)

    n_pairs_init = len(pairs.index)

    pairs['mrca'] = pairs.apply(lambda row: get_MRCA(tree, row[f'{strain_col}.x'], row[f'{strain_col}.y']), axis=1)

    assert not any(
        pairs.duplicated(subset=['strain.x', 'strain.y'])), "Duplicates initially found within the pairs df."

    pairs = pd.merge(pairs, traits.rename(columns={'node': 'strain.x'}),
                     how='left')  # By definition, strain.x and strain.y belong to the same component --> We only need to merge it to strain.x

    assert not any(pairs.duplicated(subset=['strain.x', 'strain.y'])), "Duplicated found after merging tips"

    pairs = pd.merge(pairs, ancestral_states.rename(columns={'node': 'mrca', 'component': 'intro_component'}),
                     how='left')

    # Some pairs here are duplicated: that's normal. Since we use the MPPA prediction method, the algo sometimes decides
    # to keep multiple states by ancestral node. It does so only when node have similar and high probabilities.
    # Here we assume that if the MRCA of a contact pair contains the state of components of those sequences among
    # its states, those sequences are clustered together.

    print(
        f'Did PastML attribute several states for some internal nodes ? {any(pairs.duplicated(subset=["source_id", "target_id"], keep="first"))}')

    pairs['phylogeny_validation'] = pairs['component'] == pairs['intro_component']

    pairs = pairs.sort_values(by=['component', 'source_id', 'target_id', 'phylogeny_validation'])
    pairs = remove_duplicated(pairs, ['source_id', 'target_id'], keep='last')

    assert n_pairs_init == len(
        pairs.index), "The number of pairs after processing does not match the initial number of pairs."

    pairs.to_csv(pairs_phylo, index=False)

