import time
import sys
import numpy as np
import pandas as pd
import gfapy
from gfa_utils import generate_matrix, node_position_handler, compute_coreness, compute_coreness_stats, create_heatmap


def main(gfa_path, output_path):
    """
    This script reads a GFA file and generates a heatmap of the nodes and genomes.
    It also outputs csv files containing the sequences of the nodes and the coreness statistics.

    Args:
        gfa_path (str): Path to the GFA file

    Returns:
        None
    """

    st = time.time()

    # Read GFA file
    print("Reading GFA file...")
    gfa = gfapy.Gfa.from_file(gfa_path)

    # Generate node presence matrix
    print("Generating data...")
    node_presence_matrix = generate_matrix(gfa)
    node_presence_matrix.to_csv(output_path + "matrix.csv", sep=";")

    # Get the total sum of columns
    total_node_occurrence = np.sum(np.abs(node_presence_matrix), axis=0)

    # Create a dataframe of all nodes and their sequence
    segment_list = pd.DataFrame({'name': [int(s.name) for s in gfa.segments],
                                 'sequence': [s.sequence for s in gfa.segments]})
    segment_list.set_index('name', inplace=True)

    # Output a csv of the nodes and their sequences
    segment_list.to_csv(output_path + "nodes.csv", sep=';')

    # Create a list of all the lengths of all the nodes in the gfa file
    sequence_lengths = [len(sequence) for sequence in segment_list.sequence]

    # Get the names of the genomes and the names of the nodes
    genomes = gfa.path_names
    nodes = gfa.segment_names

    # Calculate and save coreness
    coreness = compute_coreness(total_node_occurrence, sequence_lengths, len(genomes))
    coreness_stats = compute_coreness_stats(coreness, gfa, segment_list)
    coreness_stats.to_csv(output_path + "coreness_stats.csv")

    # Only execute this part if the dataset is small; File size of heatmap is pretty much directly correlated to this
    if len(nodes)*len(genomes) < 100000:

        # Create heatmap and save as HTML and PNG
        print("Creating heatmap...")

        # Get node positions
        start_pos_matrix, end_pos_matrix = node_position_handler(segment_list, gfa)

        # Plot everything
        fig = create_heatmap(genomes, nodes, sequence_lengths, node_presence_matrix, total_node_occurrence, coreness,
                             start_pos_matrix, end_pos_matrix)
        fig.write_html(output_path + "heatmap.html")
        fig.write_image(output_path + "heatmap.png", width=1800)

    # Print execution time and number of nodes plotted
    et = time.time()
    elapsed_time = et - st
    print(f"Done! \nExecution time: {str(elapsed_time)}\nNodes plotted: {str(len(nodes) * len(genomes))}")


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
