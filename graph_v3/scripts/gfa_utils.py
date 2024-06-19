import gfapy
import pandas as pd
import numpy as np
import concurrent.futures
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from joblib import Parallel, delayed


def create_heatmap(genomes, nodes, sequence_lengths, data, col_totals, coreness, start_pos_matrix, end_pos_matrix):
    """
    Creates a Plotly heatmap of node presence in genomes.

    Args:
        genomes (list): List of genome names.
        nodes (list): List of node names.
        sequence_lengths (list): List of integers indicating the length of each node.
        data (pandas.DataFrame): A binary matrix where each row represents a genome and each column represents a node.
            A 1 indicates that the node is present in the genome, while a 0 indicates that it is not present.
        col_totals (list): A list of integers representing the total number of genomes in which each node is present.
        coreness (list): A list of integers representing the coreness of each node.
        start_pos_matrix (pandas.DataFrame): A matrix of start positions for each node in each genome.
        end_pos_matrix (pandas.DataFrame): A matrix of end positions for each node in each genome.

    Returns:
        A Plotly Figure object representing the heatmap of node presence in genomes.
    """
    # Define the tile positions and width
    tile_position = 0

    tile_positions = []
    tile_widths = []

    for length in sequence_lengths:
        tile_positions.append(tile_position)
        tile_widths.append(length)
        tile_position += length

    tile_positions.append(tile_position)  # Add the last position

    # Define a discrete color scale
    colors = ["#32CD32", "#D3D3D3", "#1f77b4"]

    if -1 in data.values:
        colorscale = [[0, colors[0]], [0.5, colors[1]], [1, colors[2]]]  # -1 is green, 0 is white, 1 is blue
    else:
        colorscale = [[0, colors[1]], [1, colors[2]]]  # -1 is green, 0 is white, 1 is blue

    # Create hovertext for every point on the heatmap
    hovertext = np.empty((len(genomes), len(nodes)), dtype=object)  # hovertext should be same size as data
    for i, row in enumerate(data.to_numpy()):
        for j, val in enumerate(row):
            if abs(val) == 1:  # Absolute value because reverse nodes are -1
                hovertext[i][j] = \
                    f'Genome: {genomes[i]}<br>Node: {nodes[j]}<br>Length: {sequence_lengths[j]}<br>' + \
                    f'Coreness: {coreness[j]}<br>Start pos: {start_pos_matrix[j+1][i]} bp<br>' + \
                    f'End pos: {end_pos_matrix[j+1][i]} bp'  # +1 because this dataframe has the nodes and therefore
                # does not start with 0.
            else:
                hovertext[i][j] = \
                    f'Genome: {genomes[i]}<br>Node: {nodes[j]}<br>Length: {sequence_lengths[j]} (not present)'

    # Create the main heatmap of 1s and 0s representing node presence in the genomes
    fig1 = go.Figure(data=go.Heatmap(
        z=data,
        x=tile_positions,
        y=genomes,
        colorscale=colorscale,
        hovertext=hovertext,  # Need to reference this again or it doesn't work
        hovertemplate='%{hovertext}',
        name='',  # To get rid of 'trace0' next to hovertext
        showscale=False  # A legend for the main heatmap is useless as it is binary
    ))
    # Create the second heatmap of total presence of nodes in genomes
    fig2 = go.Figure(data=go.Heatmap(
        z=[col_totals],
        x=tile_positions,
        y=[''],
        colorscale='Reds',
        colorbar=dict(title='Total'),
        hovertext=[
            [f'Total: {t}<br>Node: {nodes[i]}<br>Length: {sequence_lengths[i]}<br>Coreness: {coreness[i]}' for i, t in
             enumerate(col_totals)]],
        hovertemplate='%{hovertext}',
        name=''  # To get rid of 'trace0' next to hovertext
    ))

    # Add both heatmaps to a subplot
    height2 = 1 / (len(genomes) + 1)  # The height of the total row should be proportional to the number of genomes
    height1 = 1 - height2  # The main heatmap takes up the rest of the space

    fig = make_subplots(rows=2, cols=1, shared_xaxes=True, vertical_spacing=0, row_heights=(height1, height2),
                        subplot_titles=['Genomes and Sequences', 'Total node presence'])

    fig.append_trace(fig1.data[0], 1, 1)
    fig.append_trace(fig2.data[0], 2, 1)

    # Update the layout
    fig.update_layout(
        title='Node presence in genomes',
        xaxis_title='Total node presence',  # This is actually what you see in the total row
        yaxis_title='Genomes',
        font=dict(family='Open Sans, sans-serif', size=16),
        margin=dict(l=100, r=50, t=100, b=50),
        plot_bgcolor='#f7f7f7',
        xaxis=dict(showgrid=False, ticks='', showticklabels=True),
        yaxis=dict(showgrid=False, ticks='', showticklabels=True),
        showlegend=True,
        legend=dict(
            title='',
            yanchor="top",
            y=0.95,
            xanchor="right",
            x=1.05
        ),
        annotations=[
            dict(
                x=1.15,
                y=0.5,
                xref="paper",
                yref="paper",
                text="Total node presence",
                showarrow=False
            )
        ]
    )
    return fig


def get_node_positions(segment_list, path, start_pos_matrix, end_pos_matrix):
    """
    Calculates the start and end positions of nodes in a given path.

    Args:
    - segment_list: A Pandas DataFrame containing information about DNA segments.
    - path: An object representing a path of segments, with captured_segments and segment_names attributes.
    - start_pos_matrix: A Pandas DataFrame to store the start positions of each node.
    - end_pos_matrix: A Pandas DataFrame to store the end positions of each node.

    Returns:
    - None. The start_pos_matrix and end_pos_matrix are updated in-place.
    """

    node_list = []
    for segment in path.segment_names:
        node_list.append(int(segment.name))

    if not all(node_id in segment_list.index for node_id in node_list):
        raise ValueError("Some nodes in path not present in segment_list")

    prev_end = 0
    for segment in path.captured_segments:
        node_id = int(segment.name)
        start = prev_end + 1
        start_pos_matrix.loc[path.name, node_id] = start
        length = len(segment_list.loc[node_id].sequence)
        end = start + length - 1
        end_pos_matrix.loc[path.name, node_id] = end
        prev_end = end


def node_position_handler(segment_list, gfa):
    """
    Given a segment list and a GFA object, this function generates two matrices with the starting and ending positions
    of the nodes for each genome.

    Args:
        segment_list (pd.DataFrame): A dataframe representing the segments in the GFA file.
        gfa (GFA): A GFA object representing the graph.

    Returns:
        Tuple[pd.DataFrame, pd.DataFrame]: A tuple containing two dataframes: the first one has the starting position of
        each node for each genome, while the second one has the ending position of each node for each genome.
    """

    # Create empty matrices to hold the node positions for each genome
    start_pos_matrix = pd.DataFrame(columns=[int(x) for x in gfa.segment_names], index=gfa.path_names, dtype='int32')
    start_pos_matrix.columns.name = 'nodes'
    start_pos_matrix.index.name = 'genomes'
    end_pos_matrix = start_pos_matrix.copy()

    # Iterate over the paths in the GFA object
    for path in gfa.paths:
        # Generate a dataframe with the starting and ending positions of each node in the path
        get_node_positions(segment_list, path, start_pos_matrix, end_pos_matrix)

    return start_pos_matrix, end_pos_matrix



def convert_path_to_binary_list(node_count, path):
    """Converts a path to a binary list indicating which nodes are present in the path.

    Args:
        node_count (int): Total number of nodes.
        path (gfapy.Gfa.Path): Path to convert.

    Returns:
        tuple: A tuple containing the path name and binary list indicating which nodes are present in the path.
    """

    # Extract the node IDs and orientations from the path
    nodes_with_orientation = path.to_list()[2].split(",")
    node_ids = {int(node.replace("+", "").replace("-", "")) for node in nodes_with_orientation}


    # Create a binary list where 1 indicates the presence of a node in the path
    binary_list = np.array([int(i + 1 in node_ids) for i in range(node_count)], dtype=np.int8)

    # Invert the values of nodes that have a "-" orientation in the path
    inverted_nodes = np.array([int(node.replace("-", "")) for node in nodes_with_orientation if "-" in node],
                              dtype=np.int32)
    binary_list[inverted_nodes - 1] *= -1

    return path.name, binary_list


def generate_matrix(gfa: gfapy.Gfa) -> pd.DataFrame:
    """Generates a heatmap from GFA object.

    Args:
        gfa (gfapy.Gfa): GFA object.

    Returns:
        pandas.DataFrame: Binary heatmap with the present nodes for each path.
    """
    node_presence_matrix = pd.DataFrame(columns=gfa.segment_names, index=gfa.path_names, dtype='int8')
    node_presence_matrix.columns.name = 'nodes'
    node_presence_matrix.index.name = 'genomes'
    node_count = len(gfa.segment_names)

    # Multithread the binary list creation
    with concurrent.futures.ThreadPoolExecutor() as executor:
        futures = [executor.submit(convert_path_to_binary_list, node_count, path) for path in gfa.paths]
        for future in concurrent.futures.as_completed(futures):
            if future.exception() is not None:  # if broke
                print(f'An exception occurred: {future.exception()}')
            else:
                path_name, binary_list = future.result()
                node_presence_matrix.loc[path_name] = binary_list

    return node_presence_matrix


def compute_coreness(total_occurrence: pd.Series, sequence_lengths: pd.Series, genome_count: int) -> pd.DataFrame:
    """
    Computes the coreness of each node in a graph based on its total occurrence in all genomes.

    Parameters:
    - total_occurrence (pandas.Series): A Series containing the total occurrence of each node in all genomes.
    - sequence_lengths (pandas.Series): A Series containing the length of each genome sequence.
    - genome_count (int): The total number of genomes.

    Returns:
    - pandas.DataFrame: A DataFrame containing the total sequence length and percentage of nodes for each coreness state.
                        Coreness states are 'core', 'soft_core', 'accessory', and 'unique'.
    """

    # Define the thresholds for core, soft core, and unique nodes
    unique_threshold = 1
    soft_core_threshold = genome_count - 1
    core_threshold = genome_count

    # Compute the coreness state for each node based on its total occurrence
    coreness = total_occurrence.apply(
        lambda count: "unique" if count == unique_threshold else
        "soft_core" if count == soft_core_threshold else
        "core" if count == core_threshold else
        "accessory"
    )

    return coreness


def coreness_per_genomes(coreness, path, segment_list):
    """
    Computes the percentage of nodes with each coreness state in a given path.

    Parameters:
    - coreness (list): A list of integers indicating the coreness state of each node in a graph.
                      State can be one of 'core', 'soft_core', 'accessory', 'unique'.
    - path (gfapy.Gfa.path): A path object representing a connected sequence of nodes in a graph.

    Returns:
    - dict: A dictionary with the percentage of nodes in each coreness state.
           Keys are 'core', 'soft_core', 'accessory', and 'unique'.
    """

    # Extract node IDs from the path object
    nodes_with_orientation = path.to_list()[2].split(",")
    node_ids = [int(node.replace("+", "").replace("-", "")) for node in nodes_with_orientation]

    # Count the number of nodes in each coreness state
    counts = {
        'core': 0,
        'soft_core': 0,
        'accessory': 0,
        'unique': 0
    }
    for node in node_ids:
        state = coreness[node - 1]  # State is either core, accessory, etc.
        node_length = len(segment_list.loc[node].sequence)
        counts[state] += node_length  # Count up corresponding state in dict

    # Compute the percentage of nodes in each coreness state
    total = sum(counts.values())
    for key in counts:
        counts[key] = counts[key] / total * 100

    df = pd.DataFrame(counts.values(), index=counts.keys(), columns=['percentage']).T

    return df


def compute_coreness_stats(coreness, gfa, segment_list):
    """
    Returns a dictionary in the format expected by a specific tool called MultiQC.
    The MultiQC tool generates quality control reports based on the input data.
    The dictionary contains the coreness counts for each input path.

    Args:
        coreness (pandas.Series): A series object containing the coreness state for each node in the graph.
        paths (gfapy.Gfa.paths): A list of Path objects representing the input files.

    Returns:
        dict: A dictionary in the format expected by the MultiQC tool.
    """

    coreness_types = ['core', 'soft_core', 'accessory', 'unique']
    data = pd.DataFrame(index=gfa.path_names, columns=coreness_types)
    data.index.name = "genome"

    for path in gfa.paths:
        df = coreness_per_genomes(coreness, path, segment_list)
        data.loc[path.name] = df.iloc[0]

    return data
