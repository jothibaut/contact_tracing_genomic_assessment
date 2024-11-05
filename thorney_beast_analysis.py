'''
This file contains the function to parse the ThorneyBEAST template, which aims to fit the tree in time.
It also contains a function to check whether we correctly guessed the tree root age, within the XML file.
'''
from commons import *

import pandas as pd
import math


'''
@template: XML template file for Thorney BEAST analysis.
@resolved_timetree: Path to the resolved timetree NEWICK file.
@unresolved_divergence_tree: Path to the unresolved divergence tree NEWICK file.
@metadata: Metadata TSV file only containing sequences within the trees.
@run_name: [String] The name you give to the run.
@xml_out: The configured XML file, ready to run with Thorney BEAST, with the specified data.
@root_time_guess_file: LOG file in which the guessed root time will be printed.
@chain_length: [Integer] Length of Markov Chain.
@stats_logEvery: [Integer]: Number of iterations within Markov Chain before displaying priors and parameters on the screen and in log files.
@trees_logEvery: [Integer]: Number of iterations within Markov Chain before printing the tree in trees file.
This function creates a Thorney BEAST configuration file, starting from the @template XML file.
'''
def configure_xml(template, resolved_timetree, unresolved_divergence_tree, metadatafile, run_name, xml_out, root_time_guess_file, chain_length=2000000000, stats_logEvery=20000, trees_logEvery=1000000):
    metadata = pd.read_csv(metadatafile, delimiter='\t')

    # Evaluate the timespan of the tree that will be inferred.
    # Here, we have decided to discretize the skygrid model by intervals of one week.
    r_time_resolved_timetree = root_time(resolved_timetree, metadatafile)
    r_time_unresolved_divergence_tree = root_time(unresolved_divergence_tree, metadatafile)
    r_time = min(r_time_resolved_timetree, r_time_unresolved_divergence_tree) # We take the time of the oldest tree root.
    r_time = float(f"{r_time:.2f}") - 0.1 # We lower the root age by 0.1, as a security factor, and only keep 2 decimals.
    os.system(f"echo {r_time} > {root_time_guess_file}")
    most_recent_date = decimal_date(metadata['date'].max()) # Age of the youngest leaf.
    time_window = most_recent_date - r_time # Time interval covering the whole tree, including security factor.
    n_weeks = int(math.ceil(time_window * 53)) # There are 52.142857143 weeks in a non-leap year, we took 53 to get a time interval of a bit less than a week (6.89 days).

    # Get the NEWICK trees as strings.
    with open(resolved_timetree, 'r') as file:
        resolved_timetree = file.read()
    with open(unresolved_divergence_tree, 'r') as file:
        unresolved_divergence_tree = file.read()

    # Configure the XML file.
    with open(template, 'r') as infile, open(xml_out, 'w') as outfile:
        for line in infile:
            if 'TEMPLATE' in line:
                line = line.replace('TEMPLATE', run_name)
            if 'CHAIN_LENGTH' in line:
                line = line.replace('CHAIN_LENGTH', str(chain_length))
            if 'STATS_LOGEVERY' in line:
                line = line.replace('STATS_LOGEVERY', str(stats_logEvery))
            if 'TREES_LOGEVERY' in line:
                line = line.replace('TREES_LOGEVERY', str(trees_logEvery))
            if 'SKYGRID_LOGPOPSIZE' in line:
                line = line.replace('SKYGRID_LOGPOPSIZE', str(n_weeks)) # Number of intervals for SkyGrid = number of weeks in the study period.
            if 'SKYGRID_NUMGRIDPOINTS' in line:
                line = line.replace('SKYGRID_NUMGRIDPOINTS', str(n_weeks-1)) # Number of intervals to compute the population size == skygrid.logPopSize - 1, because we will use the most recent date for the last one
            if 'SKYGRID_CUTOFF' in line:
                line = line.replace('SKYGRID_CUTOFF', str(time_window)) # Time window from approximated root to most recent date

            # Those condition gates work since, if one of the (7) conditions above is met, it will always enter the else gate below.
            # (Based on the XML template file structure.)
            if "<!--insert taxon block here-->" in line:
                for row in metadata.itertuples():
                    outfile.write(f'\t\t<taxon id="{row.name}">\n')
                    outfile.write(f'\t\t\t<date value="{decimal_date(row.date)}" direction="forwards" units="years"/>\n')
                    outfile.write(f'\t\t</taxon>\n')
            elif "<!--insert resolved timetree in newick here-->" in line:
                outfile.write(f'\t\t{resolved_timetree}')
            elif "<!--insert unresolved divergence tree in newick here--> " in line:
                outfile.write(f'\t\t{unresolved_divergence_tree}')
            else:
                outfile.write(line)

    print(f"Processing complete. Modified file saved as '{xml_out}'.")


'''
@root_time_guess_file: log file only containing the guessed time of Thorney BEAST's tree root, during the configuration of the corresponding XML file.
@tb_tree: Thorney BEAST's time-scaled tree, in NEWICK format.
@metadatafile: Metadata which sequences match both trees.
This function checks whether the root age computed in the configure_xml() function falls before the actual root age of the tree inferred with Thorney BEAST.
'''
def verify_root_time(root_time_guess_file, tb_tree, metadatafile):
    with open(root_time_guess_file, 'r') as file:
        root_time_guess = float(file.read())

    tb_root_time = root_time(tb_tree, metadatafile)

    is_tb_analysis_valid = root_time_guess < tb_root_time

    print(f"Is the guessed Thorney BEAST root time before its actual root time ({root_time_guess:.3f} < {tb_root_time:.3f})? {is_tb_analysis_valid}")
    if not is_tb_analysis_valid:
        print('Thorney BEAST inference must be rerun with a better guess for root time.')
    else:
        print('The root time was correctly guessed when configuring the Thorney BEAST XML file.')

