'''
This file is the main analysis script.
In this script, we call all the functions necessary to run the genomic assessment analysis of contact tracing networks.
'''
from pastml_analysis import *
from snp_analysis import *
from thorney_beast_analysis import *

import subprocess


# TODO: Update those 'executor' variables according to your installation.
# Link to your executable files.
FASTTREE_EXE = '/YOUR_PATH/FastTreeMP'
IQTREE_EXE = '/YOUR_PATH/iqtree'
SNP_DISTS_EXE = '/YOUR_PATH/snp-dists'
TEMPEST_EXE = '/YOUR_PATH/TempEst_v1.5.3/bin/tempest'
THORNEY_BEAST_EXE = '/YOUR_PATH/BEASTv1.10.5pre_thorney_0.1.2/bin/beast'
TRACER_EXE = '/YOUR_PATH/Tracer_v1.7.2/bin/tracer'
TREE_ANNOTATOR_EXE = '/YOUR_PATH/BEASTv1.10.5pre_thorney_0.1.2/bin/treeannotator'

if __name__ == '__main__':
    print('')

    # -----------------------------------
    # (D) Selecting informative sequences
    # -----------------------------------

    # Infer Fasttree
    os.system(f"{FASTTREE_EXE} -gtr -nt -gamma -log {os.path.join(ANALYSES, TREES, 'ba1_fasttree.log')} {os.path.join(DATA, GENOMICS, 'ba1_all.fasta')} > {os.path.join(ANALYSES, TREES, 'ba1_fasttree.nwk')}")
    os.system(f"{FASTTREE_EXE} -gtr -nt -gamma -log {os.path.join(ANALYSES, TREES, 'ba2_fasttree.log')} {os.path.join(DATA, GENOMICS, 'ba2_all.fasta')} > {os.path.join(ANALYSES, TREES, 'ba2_fasttree.nwk')}")


    # Down-sample Fasttree
    subprocess.run(["Rscript", "downsample_monophyletic_cluster.R"], check=True)
    keep_sequences_from_metadata(os.path.join(DATA, GENOMICS, "ba1_subsampled.tsv"),
                                 os.path.join(DATA, GENOMICS, "ba1_all.fasta"),
                                 os.path.join(DATA, GENOMICS, "ba1_subsampled.fasta"))
    keep_sequences_from_metadata(os.path.join(DATA, GENOMICS, "ba2_subsampled.tsv"),
                                 os.path.join(DATA, GENOMICS, "ba2_all.fasta"),
                                 os.path.join(DATA, GENOMICS, "ba2_subsampled.fasta"))

    # ----------------------------------------------------------------
    # E. Assessing case-contact pairs based on phylogenetic clustering
    # ----------------------------------------------------------------

    # Run IQ-tree on down-sampled genomic data set.
    subprocess.run(["iqtree", "-s", os.path.join(DATA, GENOMICS, "ba1_subsampled.fasta"), "-m", "GTR", "-alrt", "1000", "-B", "1000", "-T", "AUTO"])
    subprocess.run(["iqtree", "-s", os.path.join(DATA, GENOMICS, "ba2_subsampled.fasta"), "-m", "GTR", "-alrt", "1000", "-B", "1000", "-T", "AUTO"])

    # Exclude outlier tips, using tempEST
    # See tutorial here: https://beast.community/tempest_tutorial
    os.system(f"{TEMPEST_EXE}")

    drop_tree_sequences(os.path.join(ANALYSES, TEMPEST, 'ba1_tempest_outliers.csv'),
                        os.path.join(ANALYSES, TREES, 'ba1_iqtree.nwk'),
                        os.path.join(ANALYSES, TREES, 'ba1_iqtree_no_outliers.nwk'))
    drop_metadata_sequences(os.path.join(ANALYSES, TEMPEST, 'ba1_tempest_outliers.csv'),
                            os.path.join(DATA, GENOMICS, 'ba1_subsampled.tsv'),
                            os.path.join(DATA, GENOMICS, 'ba1_subsampled_no_tempest_outliers.tsv'))
    drop_fasta_sequences(os.path.join(ANALYSES, TEMPEST, 'ba1_tempest_outliers.csv'),
                         os.path.join(DATA, GENOMICS, 'ba1_subsampled.fasta'),
                         os.path.join(DATA, GENOMICS, 'ba1_subsampled_no_tempest_outliers.fasta'))

    drop_tree_sequences(os.path.join(ANALYSES, TEMPEST, 'ba2_tempest_outliers.csv'),
                        os.path.join(ANALYSES, TREES, 'ba2_iqtree.nwk'),
                        os.path.join(ANALYSES, TREES, 'ba2_iqtree_no_outliers.nwk'))
    drop_metadata_sequences(os.path.join(ANALYSES, TEMPEST, 'ba2_tempest_outliers.csv'),
                            os.path.join(DATA, GENOMICS, 'ba2_subsampled.tsv'),
                            os.path.join(DATA, GENOMICS, 'ba2_subsampled_no_tempest_outliers.tsv'))
    drop_fasta_sequences(os.path.join(ANALYSES, TEMPEST, 'ba2_tempest_outliers.csv'),
                         os.path.join(DATA, GENOMICS, 'ba2_subsampled.fasta'),
                         os.path.join(DATA, GENOMICS, 'ba2_subsampled_no_tempest_outliers.fasta'))

    # Fit tree in time. See Alternative section below for results with Thorney BEAST inference.
    os.system(f"treetime --aln {DATA}/{GENOMICS}/ba1_subsampled_no_tempest_outliers.fasta --dates {DATA}/{GENOMICS}/ba1_subsampled_no_tempest_outliers.tsv --tree {ANALYSES}/{TREES}/ba1_iqtree_no_outliers.nwk --clock-rate 0.0008 --clock-std-dev 0.0004 --clock-filter 8 --coalescent skyline --outdir {ANALYSES}/{TREETIME}/{BA1} > {ANALYSES}/{TREETIME}/{BA1}/treetime_stdout.txt 2>&1")
    os.system(f"treetime --aln {DATA}/{GENOMICS}/ba2_subsampled_no_tempest_outliers.fasta --dates {DATA}/{GENOMICS}/ba2_subsampled_no_tempest_outliers.tsv --tree {ANALYSES}/{TREES}/ba2_iqtree_no_outliers.nwk --clock-rate 0.0008 --clock-std-dev 0.0004 --clock-filter 8 --coalescent skyline --outdir {ANALYSES}/{TREETIME}/{BA2} > {ANALYSES}/{TREETIME}/{BA2}/treetime_stdout.txt 2>&1")

    # Remove Treetime outliers.
    for variant in [BA1, BA2]:
        output_treetime_outliers(os.path.join(ANALYSES, TREETIME, variant, 'treetime_stdout.txt'),
                                 os.path.join(ANALYSES, TREETIME, variant, 'treetime_outliers.csv'))

        convert_tree(os.path.join(ANALYSES, TREETIME, variant, 'timetree.nexus'),
                     os.path.join(ANALYSES, TREETIME, variant, 'timetree.nwk'))
        convert_tree(os.path.join(ANALYSES, TREETIME, variant, 'divergence_tree.nexus'),
                     os.path.join(ANALYSES, TREETIME, variant, 'divergence_tree.nwk'))

        drop_tree_sequences(os.path.join(ANALYSES, TREETIME, variant, 'treetime_outliers.csv'),
                            os.path.join(ANALYSES, TREETIME, variant, 'timetree.nwk'),
                            os.path.join(ANALYSES, TREETIME, variant, 'timetree_no_outliers.nwk'))
        drop_tree_sequences(os.path.join(ANALYSES, TREETIME, variant, 'treetime_outliers.csv'),
                            os.path.join(ANALYSES, TREETIME, variant, 'divergence_tree.nwk'),
                            os.path.join(ANALYSES, TREETIME, variant, 'divergence_tree_no_outliers.nwk'))
        drop_metadata_sequences(os.path.join(ANALYSES, TREETIME, variant, 'treetime_outliers.csv'),
                                os.path.join(DATA, GENOMICS, f'{variant}_subsampled_no_tempest_outliers.tsv'),
                                os.path.join(DATA, GENOMICS, f'{variant}_subsampled_no_treetime_outliers.tsv'))
        drop_fasta_sequences(os.path.join(ANALYSES, TREETIME, variant, 'treetime_outliers.csv'),
                             os.path.join(DATA, GENOMICS, f'{variant}_subsampled_no_tempest_outliers.fasta'),
                             os.path.join(DATA, GENOMICS, f'{variant}_subsampled_no_treetime_outliers.fasta'))

    # Clean Timetree: rename internal nodes
    for variant in [BA1, BA2]:
        rename_internal_node(os.path.join(ANALYSES, TREETIME, variant, 'timetree_no_outliers.nwk'),
                             os.path.join(ANALYSES, TREETIME, variant, 'timetree_renamed_nodes.nwk'))

    # Run PastML on clean Treetime's output.
    # BA.1
    run_pastml(metadata=os.path.join(DATA, GENOMICS, 'ba1_subsampled_no_treetime_outliers.tsv'),
               traits=os.path.join(DATA, NETWORK, 'ba1_network_components.tsv'),
               tree=os.path.join(ANALYSES, TREETIME, BA1, 'timetree_renamed_nodes.nwk'),
               html_compressed=os.path.join(ANALYSES, PASTML, BA1, "traits.html"),
               html=os.path.join(ANALYSES, PASTML, BA1, "traits_tree.html"),
               all_tips_traits=os.path.join(ANALYSES, PASTML, BA1, 'all_tips_traits.tsv'),
               prediction_method='MPPA')
    os.system(f"mv {os.path.join(ANALYSES, TREETIME, BA1, 'timetree_renamed_nodes_pastml', '*')} {os.path.join(ANALYSES, PASTML, BA1)}")
    os.system(f"rm -r {os.path.join(ANALYSES, TREETIME, BA1, 'timetree_renamed_nodes_pastml')}")

    get_pastml_clusters_mppa(traits=os.path.join(DATA, NETWORK, 'ba1_network_components.tsv'),
                             ancestral_states=os.path.join(ANALYSES, PASTML, BA1, 'combined_ancestral_states.tab'),
                             pairs=os.path.join(DATA, NETWORK, 'ba1_network_pairs.csv'),
                             tree=os.path.join(ANALYSES, TREETIME, BA1, 'timetree_renamed_nodes.nwk'),
                             pairs_phylo=os.path.join(ANALYSES, GENOMIC_VALIDATION, 'ba1_network_pairs_pastml.csv'))

    # BA.2
    run_pastml(metadata=os.path.join(DATA, GENOMICS, 'ba2_subsampled_no_treetime_outliers.tsv'),
               traits=os.path.join(DATA, NETWORK, 'ba2_network_components.tsv'),
               tree=os.path.join(ANALYSES, TREETIME, BA2, 'timetree_renamed_nodes.nwk'),
               html_compressed=os.path.join(ANALYSES, PASTML, BA2, "traits.html"),
               html=os.path.join(ANALYSES, PASTML, BA2, "traits_tree.html"),
               all_tips_traits=os.path.join(ANALYSES, PASTML, BA2, 'all_tips_traits.tsv'),
               prediction_method='MPPA')
    os.system(f"mv {os.path.join(ANALYSES, TREETIME, BA2, 'timetree_renamed_nodes_pastml', '*')} {os.path.join(ANALYSES, PASTML, BA2)}")
    os.system(f"rm -r {os.path.join(ANALYSES, TREETIME, BA2, 'timetree_renamed_nodes_pastml')}")

    get_pastml_clusters_mppa(traits=os.path.join(DATA, NETWORK, 'ba2_network_components.tsv'),
                             ancestral_states=os.path.join(ANALYSES, PASTML, BA2, 'combined_ancestral_states.tab'),
                             pairs=os.path.join(DATA, NETWORK, 'ba2_network_pairs.csv'),
                             tree=os.path.join(ANALYSES, TREETIME, BA2, 'timetree_renamed_nodes.nwk'),
                             pairs_phylo=os.path.join(ANALYSES, GENOMIC_VALIDATION, 'ba2_network_pairs_pastml.csv'))

    # --------------------------------------------------------
    # F. Correcting phylogenetic clustering with SNPs analysis
    # --------------------------------------------------------

    fasta_subset(os.path.join(DATA, NETWORK, 'ba1_network_components.tsv'),
                 os.path.join(DATA, GENOMICS, 'ba1_all.fasta'),
                 os.path.join(DATA, GENOMICS, 'ba1_ct.fasta'))
    fasta_subset(os.path.join(DATA, NETWORK, 'ba2_network_components.tsv'),
                 os.path.join(DATA, GENOMICS, 'ba2_all.fasta'),
                 os.path.join(DATA, GENOMICS, 'ba2_ct.fasta'))

    os.system(f"{SNP_DISTS_EXE} {os.path.join(DATA, GENOMICS, 'ba1_ct.fasta')} > {os.path.join(ANALYSES, SNP, 'ba1_ct.dist')}")
    os.system(f"{SNP_DISTS_EXE} {os.path.join(DATA, GENOMICS, 'ba2_ct.fasta')} > {os.path.join(ANALYSES, SNP, 'ba2_ct.dist')}")

    # Creating assessment files with Treetime runs.

    add_snp_to_pairs(pairs_in=os.path.join(ANALYSES, GENOMIC_VALIDATION, 'ba1_network_pairs_pastml.csv'),
                     pairs_out=os.path.join(ANALYSES, GENOMIC_VALIDATION, 'ba1_network_pairs_pastml_snp.csv'),
                     distances=os.path.join(ANALYSES, SNP, 'ba1_ct.dist'))
    add_snp_to_pairs(pairs_in=os.path.join(ANALYSES, GENOMIC_VALIDATION, 'ba2_network_pairs_pastml.csv'),
                     pairs_out=os.path.join(ANALYSES, GENOMIC_VALIDATION, 'ba2_network_pairs_pastml_snp.csv'),
                     distances=os.path.join(ANALYSES, SNP, 'ba2_ct.dist'))

    merge_all_pairs(ba1_pairs=os.path.join(ANALYSES, GENOMIC_VALIDATION, 'ba1_network_pairs_pastml_snp.csv'),
                    ba2_pairs=os.path.join(ANALYSES, GENOMIC_VALIDATION, 'ba2_network_pairs_pastml_snp.csv'),
                    ba1ba2_pairs=os.path.join(DATA, NETWORK, 'diff_network_pairs.csv'),
                    all_pairs_phylo=os.path.join(ANALYSES, GENOMIC_VALIDATION, 'all_network_pairs_pastml_snp.csv'))

    clustering_on_snp_graph(ct_pairs=os.path.join(ANALYSES, GENOMIC_VALIDATION, 'all_network_pairs_pastml_snp.csv'),
                            out_graph=os.path.join(FIGURES, "distances_distribution.svg"))





    #===================================================================================================================



    #######################################
    # ALTERNATIVE PIPELINE: THORNEY BEAST #
    #######################################

    # ----------------------------------------------------------------
    # E. Assessing case-contact pairs based on phylogenetic clustering
    # ----------------------------------------------------------------

    # First, we need to build the divergence tree with polytomies
    os.system(f"treetime --aln {DATA}/{GENOMICS}/ba1_subsampled_no_tempest_outliers.fasta --dates {DATA}/{GENOMICS}/ba1_subsampled_no_tempest_outliers.tsv --tree {ANALYSES}/{TREES}/ba1_iqtree_no_outliers.nwk --clock-rate 0.0008 --clock-std-dev 0.0004 --clock-filter 8 --coalescent skyline --keep-polytomies --outdir {ANALYSES}/{TREETIME}/{BA1_POLYTOMIES} > {ANALYSES}/{TREETIME}/{BA1_POLYTOMIES}/treetime_stdout.txt 2>&1")
    os.system(f"treetime --aln {DATA}/{GENOMICS}/ba2_subsampled_no_tempest_outliers.fasta --dates {DATA}/{GENOMICS}/ba2_subsampled_no_tempest_outliers.tsv --tree {ANALYSES}/{TREES}/ba2_iqtree_no_outliers.nwk --clock-rate 0.0008 --clock-std-dev 0.0004 --clock-filter 8 --coalescent skyline --keep-polytomies --outdir {ANALYSES}/{TREETIME}/{BA2_POLYTOMIES} > {ANALYSES}/{TREETIME}/{BA2_POLYTOMIES}/treetime_stdout.txt 2>&1")

    # Remove outliers based on Treetime inference, with polytomies resolution (See code above).
    for variant in [BA1, BA2]:
        variant_polytomies = f'{variant}_polytomies'
        convert_tree(os.path.join(ANALYSES, TREETIME, variant_polytomies, 'timetree.nexus'),
                     os.path.join(ANALYSES, TREETIME, variant_polytomies, 'timetree.nwk'))
        convert_tree(os.path.join(ANALYSES, TREETIME, variant_polytomies, 'divergence_tree.nexus'),
                     os.path.join(ANALYSES, TREETIME, variant_polytomies, 'divergence_tree.nwk'))

        drop_tree_sequences(os.path.join(ANALYSES, TREETIME, variant, 'treetime_outliers.csv'),
                            os.path.join(ANALYSES, TREETIME, variant_polytomies, 'timetree.nwk'),
                            os.path.join(ANALYSES, TREETIME, variant_polytomies, 'timetree_no_outliers.nwk'))
        drop_tree_sequences(os.path.join(ANALYSES, TREETIME, variant, 'treetime_outliers.csv'),
                            os.path.join(ANALYSES, TREETIME, variant_polytomies, 'divergence_tree.nwk'),
                            os.path.join(ANALYSES, TREETIME, variant_polytomies, 'divergence_tree_no_outliers.nwk'))
        remove_zero_branches(os.path.join(ANALYSES, TREETIME, variant, 'timetree_no_outliers.nwk'),
                             os.path.join(ANALYSES, TREETIME, variant, 'timetree_no_outliers_no_zero_branches.nwk'))

    # Run Thorney BEAST
    # BA.1
    configure_xml(template=os.path.join(ANALYSES, THORNEY_BEAST, 'covid_skygrid_template.xml'),
                  resolved_timetree=os.path.join(ANALYSES, TREETIME, BA1, 'timetree_renamed_nodes.nwk'),
                  unresolved_divergence_tree=os.path.join(ANALYSES, TREETIME, BA1_POLYTOMIES, 'divergence_tree_no_outliers.nwk'),
                  metadatafile=os.path.join(DATA, GENOMICS, 'ba1_subsampled_no_treetime_outliers.tsv'),
                  run_name='ba1_long',
                  xml_out=os.path.join(ANALYSES, THORNEY_BEAST, BA1, 'ba1_long.xml'),
                  root_time_guess_file=os.path.join(ANALYSES, THORNEY_BEAST, BA1, 'root_time_guess.log'))

    os.chdir(os.path.join(ANALYSES, THORNEY_BEAST, BA1))
    os.system(f"{THORNEY_BEAST_EXE} -beagle_gpu -beagle_double -beagle_order 2 ba1_long.xml")
    os.chdir("../../..")

    os.system(f"{TRACER_EXE} {os.path.join(ANALYSES, THORNEY_BEAST, BA1, 'ba1_long.log')}")

    os.system(f"{TREE_ANNOTATOR_EXE} -burnin 100000000 {os.path.join(ANALYSES, THORNEY_BEAST, BA1, 'ba1_long.trees')} {os.path.join(ANALYSES, THORNEY_BEAST, BA1, 'ba1_long_mcc.nexus')}")


    # BA.2
    configure_xml(template=os.path.join(ANALYSES, THORNEY_BEAST, 'covid_skygrid_template.xml'),
                  resolved_timetree=os.path.join(ANALYSES, TREETIME, BA2, 'timetree_renamed_nodes.nwk'),
                  unresolved_divergence_tree=os.path.join(ANALYSES, TREETIME, BA2_POLYTOMIES, 'divergence_tree_no_outliers.nwk'),
                  metadatafile=os.path.join(DATA, GENOMICS, 'ba2_subsampled_no_treetime_outliers.tsv'),
                  run_name='ba2_long',
                  xml_out=os.path.join(ANALYSES, THORNEY_BEAST, BA2, 'ba2_long.xml'),
                  root_time_guess_file=os.path.join(ANALYSES, THORNEY_BEAST, BA2, 'root_time_guess.log'))

    os.chdir(os.path.join(ANALYSES, THORNEY_BEAST, BA2))
    os.system(f"{THORNEY_BEAST_EXE} -beagle_gpu -beagle_double -beagle_order 1 ba2_long.xml")
    os.chdir("../../..")

    os.system(f"{TRACER_EXE} {os.path.join(ANALYSES, THORNEY_BEAST, BA2, 'ba2_long.log')}")

    os.system(f"{TREE_ANNOTATOR_EXE} -burnin 100000000 {os.path.join(ANALYSES, THORNEY_BEAST, BA2, 'ba2_long.trees')} {os.path.join(ANALYSES, THORNEY_BEAST, BA2, 'ba2_long_mcc.nexus')}")

    for variant in [BA1, BA2]:
        convert_tree(os.path.join(ANALYSES, THORNEY_BEAST, variant, f'{variant}_long_mcc.nexus'),
                     os.path.join(ANALYSES, THORNEY_BEAST, variant, f'{variant}_long_mcc.nwk'))
        rename_internal_node(os.path.join(ANALYSES, THORNEY_BEAST, variant, f'{variant}_long_mcc.nwk'),
                             os.path.join(ANALYSES, THORNEY_BEAST, variant, f'{variant}_long_renamed_nodes_mcc.nwk'))
        # If the actual root time is > the guessed root time, it means that we did not infer the SkyGrid model over the
        # whole tree period. We must therefore re-infer Thorney BEAST with a better root time guess, < actual root time.
        verify_root_time(os.path.join(ANALYSES, THORNEY_BEAST, variant, 'root_time_guess.log'),
                         os.path.join(ANALYSES, THORNEY_BEAST, variant, f'{variant}_long_renamed_nodes_mcc.nwk'),
                         os.path.join(DATA, GENOMICS, f'{variant}_subsampled_no_treetime_outliers.tsv'))

    # Run PastML on clean Thorney BEAST's output.
    # BA.1
    run_pastml(metadata=os.path.join(DATA, GENOMICS, 'ba1_subsampled_no_treetime_outliers.tsv'),
               traits=os.path.join(DATA, NETWORK, 'ba1_network_components.tsv'),
               tree=os.path.join(ANALYSES, THORNEY_BEAST, BA1, 'ba1_long_renamed_nodes_mcc.nwk'),
               html_compressed=os.path.join(ANALYSES, PASTML, BA1_THORNEYBEAST, "traits.html"),
               html=os.path.join(ANALYSES, PASTML, BA1_THORNEYBEAST, "traits_tree.html"),
               all_tips_traits=os.path.join(ANALYSES, PASTML, BA1_THORNEYBEAST, 'all_tips_traits.tsv'),
               prediction_method='MPPA')
    os.system(f"mv {os.path.join(ANALYSES, THORNEY_BEAST, BA1, 'ba1_long_renamed_nodes_mcc_pastml', '*')} {os.path.join(ANALYSES, PASTML, BA1_THORNEYBEAST)}")
    os.system(f"rm -r {os.path.join(ANALYSES, THORNEY_BEAST, BA1, 'ba1_long_renamed_nodes_mcc_pastml')}")

    get_pastml_clusters_mppa(traits=os.path.join(DATA, NETWORK, 'ba1_network_components.tsv'),
                             ancestral_states=os.path.join(ANALYSES, PASTML, BA1_THORNEYBEAST, 'combined_ancestral_states.tab'),
                             pairs=os.path.join(DATA, NETWORK, 'ba1_network_pairs.csv'),
                             tree=os.path.join(ANALYSES, THORNEY_BEAST, BA1, 'ba1_long_renamed_nodes_mcc.nwk'),
                             pairs_phylo=os.path.join(ANALYSES, GENOMIC_VALIDATION, 'ba1_network_pairs_pastml_tb.csv'))

    # BA.2
    run_pastml(metadata=os.path.join(DATA, GENOMICS, 'ba2_subsampled_no_treetime_outliers.tsv'),
               traits=os.path.join(DATA, NETWORK, 'ba2_network_components.tsv'),
               tree=os.path.join(ANALYSES, THORNEY_BEAST, BA2, 'ba2_long_renamed_nodes_mcc.nwk'),
               html_compressed=os.path.join(ANALYSES, PASTML, BA2_THORNEYBEAST, "traits.html"),
               html=os.path.join(ANALYSES, PASTML, BA2_THORNEYBEAST, "traits_tree.html"),
               all_tips_traits=os.path.join(ANALYSES, PASTML, BA2_THORNEYBEAST, 'all_tips_traits.tsv'),
               prediction_method='MPPA')
    os.system(f"mv {os.path.join(ANALYSES, THORNEY_BEAST, BA2, 'ba2_long_renamed_nodes_mcc_pastml', '*')} {os.path.join(ANALYSES, PASTML, BA2_THORNEYBEAST)}")
    os.system(f"rm -r {os.path.join(ANALYSES, THORNEY_BEAST, BA2, 'ba2_long_renamed_nodes_mcc_pastml')}")

    get_pastml_clusters_mppa(traits=os.path.join(DATA, NETWORK, 'ba2_network_components.tsv'),
                             ancestral_states=os.path.join(ANALYSES, PASTML, BA2_THORNEYBEAST, 'combined_ancestral_states.tab'),
                             pairs=os.path.join(DATA, NETWORK, 'ba2_network_pairs.csv'),
                             tree=os.path.join(ANALYSES, THORNEY_BEAST, BA2, 'ba2_long_renamed_nodes_mcc.nwk'),
                             pairs_phylo=os.path.join(ANALYSES, GENOMIC_VALIDATION, 'ba2_network_pairs_pastml_tb.csv'))

    # --------------------------------------------------------
    # F. Correcting phylogenetic clustering with SNPs analysis
    # --------------------------------------------------------

    # Creating assessment files with Thorney BEAST runs.

    add_snp_to_pairs(pairs_in=os.path.join(ANALYSES, GENOMIC_VALIDATION, 'ba1_network_pairs_pastml_tb.csv'),
                     pairs_out=os.path.join(ANALYSES, GENOMIC_VALIDATION, 'ba1_network_pairs_pastml_snp_tb.csv'),
                     distances=os.path.join(ANALYSES, SNP, 'ba1_ct.dist'))
    add_snp_to_pairs(pairs_in=os.path.join(ANALYSES, GENOMIC_VALIDATION, 'ba2_network_pairs_pastml_tb.csv'),
                     pairs_out=os.path.join(ANALYSES, GENOMIC_VALIDATION, 'ba2_network_pairs_pastml_snp_tb.csv'),
                     distances=os.path.join(ANALYSES, SNP, 'ba2_ct.dist'))

    merge_all_pairs(ba1_pairs=os.path.join(ANALYSES, GENOMIC_VALIDATION, 'ba1_network_pairs_pastml_snp_tb.csv'),
                    ba2_pairs=os.path.join(ANALYSES, GENOMIC_VALIDATION, 'ba2_network_pairs_pastml_snp_tb.csv'),
                    ba1ba2_pairs=os.path.join(DATA, NETWORK, 'diff_network_pairs.csv'),
                    all_pairs_phylo=os.path.join(ANALYSES, GENOMIC_VALIDATION, 'all_network_pairs_pastml_snp_tb.csv'))

    clustering_on_snp_graph(ct_pairs=os.path.join(ANALYSES, GENOMIC_VALIDATION, 'all_network_pairs_pastml_snp_tb.csv'),
                            out_graph=os.path.join(FIGURES, "distances_distribution_tb.png"))

