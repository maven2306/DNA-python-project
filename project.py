
"""
author: Matteo Venturini
student number: 2469579
"""

import pandas as pd
import networkx as nx
import graphviz
import sys
import os
import json

# ------------------- PATH SETUP (BEST PRACTICE) -------------------

# Get the absolute path of the directory where this script is located # <-- ADDED
script_dir = os.path.dirname(os.path.abspath(__file__))

# TASK 1 Read the csv file given its name

# Default behaviour is to allow the user to enter the DNA_[x]_[x].csv they want
try:
    filename_only = sys.argv[1] # <-- MODIFIED (renamed for clarity)
# Handle pytest
except IndexError:
    filename_only = 'DNA_1_5.csv' # <-- MODIFIED (renamed for clarity)

# Create the full, absolute path to the input file # <-- ADDED
full_input_path = os.path.join(script_dir, filename_only)

# Obtain the DNA sequence number and the lenght of the k-mers
base = filename_only.split('_') # <-- MODIFIED (use the filename, not the full path)

K_length = base[2][0]

new_filename = base[0] + "_" + base[1]

# Create full paths for the output files # <-- ADDED
full_png_path = os.path.join(script_dir, f'{new_filename}.png')
full_text_path = os.path.join(script_dir, f'{new_filename}.text')


# Read the data using the full path
try: # <-- ADDED (Good practice to handle file not found)
    data = pd.read_csv(full_input_path)
except FileNotFoundError:
    print(f"Error: Could not find the file at the specified path: {full_input_path}")
    sys.exit(1) # Exit the script if the file doesn't exist


# Rename the columns with the correct names
data = data.rename(columns={'1': 'SegmentNr', '1.1': 'Position', '0': 'A',
                            '0.1': 'C', '0.2': 'G', '1.2': 'T'})


# TASK 2  Clean the data of the dataframe

def clean_data(data: pd.DataFrame) -> pd.DataFrame:
    """ This function cleans the dataset.
    - If there are duplicated segments, one is removed
    - If there are missing positions within a segment, the segment is removed
    - If there are misreadings in the nucleotides (e.g more no 0s or
      more than one 1), the segment is removed
       """

    # If one position is duplicated (and share the same output), drop
    # all but one
    data = data.drop_duplicates()

    segs_to_remove = []

    unique_sequences = []

    # Find which segments to remove
    for seg_num in sorted(data['SegmentNr'].unique()):

        # Filter rows by segment number
        segment_rows = data[data['SegmentNr'] == seg_num]

        # Order the column Position
        segment_rows = segment_rows.sort_values('Position')

        # Find correct position sequence
        correct_position = ''
        position = ''

        # check whether there are missing positions
        # check whether each valid positon has a valid output
        segment_nucleotide_sequence = []

        for i in range(len(segment_rows)):
            correct_position += str(i + 1)

            pos_value = segment_rows['Position'].iloc[i]

            position += str(pos_value)

            # Extract the output
            nucleotide_values = segment_rows[segment_rows['Position'] == pos_value][['A', 'C', 'G', 'T']].values[0].tolist()

            # If the position is wrong, remove it
            if correct_position != position:
                segs_to_remove.append(seg_num)
                break
            # If the output is wrong, remove it
            elif nucleotide_values.count(1) != 1 or nucleotide_values.count(0) > 3:
                segs_to_remove.append(seg_num)
                break
            else:
                # if the positon and the output are valid, add the nucleotide to the sequence
                segment_nucleotide_sequence.append(str(nucleotide_values))

        # Check whether the segment's nucleotide sequence already exist
        if segment_nucleotide_sequence and segment_nucleotide_sequence not in unique_sequences:
            unique_sequences.append(segment_nucleotide_sequence)
        else:
            # Remove it if it exists
            segs_to_remove.append(seg_num)

    # Filter out the segments that need to be removed
    data_filtered = data.query("SegmentNr not in @segs_to_remove")

    # Reset the index, otherwise pytest will throw an error
    data_filtered = data_filtered.reset_index(drop=True)

    return data_filtered


data = clean_data(data)


# Task 3 Generate JSON sequences from the dataframe

def generate_sequences(data: pd.DataFrame) -> str:
    """
    This function generates a JSON string from the readings
    """

    segment_sequences = {}

    for seg_num in sorted(data['SegmentNr'].unique()):

        # Filter rows by segment number
        segment_rows = data[data['SegmentNr'] == seg_num]

        # Order the column Position
        segment_rows = segment_rows.sort_values('Position')

        segment_sequence = ''

        for i in range(len(segment_rows)):

            pos_value = segment_rows['Position'].iloc[i]

            # Extract the output
            nucleotide_values = segment_rows[segment_rows['Position'] == pos_value][['A', 'C', 'G', 'T']].values[0].tolist()

            # Add the nucleotide that corresponds to the value '1'
            if nucleotide_values[0] == 1:
                segment_sequence += 'A'
            elif nucleotide_values[1] == 1:
                segment_sequence += 'C'
            elif nucleotide_values[2] == 1:
                segment_sequence += 'G'
            else:
                segment_sequence += 'T'

        # Add sequences to the dictionary
        segment_sequences[str(seg_num)] = segment_sequence

    # Generate the json string
    json_str = json.dumps(segment_sequences, indent=2)

    return json_str


json_sequences = generate_sequences(data)
print(json_sequences)

# Task 4 Construct de Bruijn graph

# Helper
def _kmers_rec(s, k, index=0, k_mers=None):
    """ This function finds the k-mers of a given sequence.
        It takes the argument k which specifies the length of the k-mers.
    """

    if k_mers is None:
        k_mers = []

        if k > len(s):
            raise ValueError('Chosen k is longer than the shortest segment')

    if k-1 == len(s):
        return k_mers

    k_mers.append(s[index:k])

    return _kmers_rec(s, k+1, index+1, k_mers)


# Helper
def _LRmers(k_mers, k, index=0, i=0, LR_mers=None):
    """ This function finds the left and right -mers of each k-mer """

    if LR_mers is None:
        LR_mers = []

    for i in k_mers:
        left = i[index:k-1]
        right = i[index+1:k]
        LR_mers.append((left, right))

    return LR_mers


# Make graph
def construct_graph(s: str, k: int) -> nx.MultiDiGraph:
    """ This function creates a multidirected graph object from the
    json string, with k being the lenght of k-mers """

    D = nx.MultiDiGraph()

    # Read the json file as a dict
    sequences = json.loads(s).values()
    k_mers = []
    for seq in sequences:
        seq_kmers = _kmers_rec(seq, k)
        k_mers.extend(seq_kmers)

    LR_mers = _LRmers(k_mers, k)

    for i in LR_mers:
        l, r = i
        D.add_edge(l, r)

    return D


G = construct_graph(json_sequences, int(K_length))


# Task 5 Plot the de Bruijn graph

def plot_graph(graph: nx.MultiDiGraph, filename: str) -> None:
    """ 
    This function creates a plot and saves it as a png file
    """
    # Create a new directed graph
    dot = graphviz.Digraph(comment='De Bruijn Graph')
    dot.attr('node', shape='circle') # Optional: make nodes look nice

    # Add all nodes from the networkx graph
    for n in graph.nodes():
        dot.node(str(n))

    # Add all edges from the networkx graph
    for u, v in graph.edges():
        dot.edge(str(u), str(v))

    # Render the graph to a file. 
    # .render() will create a file named 'filename.png'
    # and a source file 'filename'. We clean up the source file.
    try:
        dot.render(filename, format='png', view=False, cleanup=True)
        print(f"Successfully created graph image at: {filename}.png")
    except graphviz.backend.execute.ExecutableNotFound:
        print("ERROR: Graphviz executable not found.")
        print("Please ensure you have installed it in your Conda environment with:")
        print("conda install -c conda-forge python-graphviz")


# The call to the function remains the same, but we pass the base filename
# without the extension.
plot_graph(G, new_filename) 


# Task 6 Check whether the de Bruijn graph can be sequenced

def is_valid_graph(graph: nx.MultiDiGraph) -> bool:
    """ This function checks whether a graph is Eulerian or has an
    Eulerian's walk. If all the nodes have degree_out == degree_in, then it's
    Eulerian. If one node has degree_out > degree_in and one node with
    degree_in > degree_out, then it has an Eulerian's walk. The obvious
    requirement is that the graph is at least weakly connected.
    """

    start = []
    end = []
    null = []
    for n in graph:
        d_in = graph.in_degree(n)
        d_out = graph.out_degree(n)
        if d_out - d_in == 1:
            start.append(n)
        if d_in - d_out == 1:
            end.append(n)
        if d_in == d_out:
            null.append(n)

    # Check if one and only one node is a starting node and a ending node, and whether it's weakly connected
    if len(start) == 1 and len(end) == 1 and len(null) == len(list(graph.nodes))-2 and nx.is_weakly_connected(G):
        return True
    # Check if there are 0 starting nodes and 0 ending nodes (full circle)
    elif len(start) == 0 and len(end) == 0 and nx.is_weakly_connected(G):
        return True
    else:
        return False


is_valid_graph(G)


# Task 7 Construct DNA sequence

def _helper_construct_dna_sequence(graph, v=None, sequence=None, eulerian=None):

    # Find the starting node
    if v is None:
        v = ''
        for n in graph:
            d_in = graph.in_degree(n)
            d_out = graph.out_degree(n)
            # If there is a starting node, then set eulerian to false
            if d_in < d_out:
                v += n
                eulerian = False
        # If there is no starting node, then use the first node created,
        # and set eulerian to true
        if v == '':
            eulerian = True
            v += list(nx.nodes(graph))[0]

    # Add the starting node to the sequence
    if sequence is None:
        sequence = ''
        sequence += v

    # End case: return the sequence when the eulerian walk is terminated
    if graph.out_degree(v) == 0:
        if eulerian is True:
            return sequence[:-1]
        else:
            return sequence

    neibs = []
    # In case there are multiple nbr, the first nbr is the one with the
    # OLDEST edge created.
    for nbr in graph[v]:
        neibs.append((nbr, graph.out_degree(nbr)))

    # Make sure that the nbr with the highest out degrees is chosen first
    sorted_neibs = sorted(neibs, key=lambda x: x[1], reverse=True)

    # From ordered tuple list, choose the first tuple, first value (nbr)
    right_neighbour = sorted_neibs[0][0]

    sequence += right_neighbour[-1]

    # Print eulerian walk
    print(f'Visiting directed edge: {v} -> {right_neighbour}')

    # With key = None, it removes the OLDEST edge created
    graph.remove_edge(v, right_neighbour, key=None)

    return _helper_construct_dna_sequence(graph, right_neighbour, sequence, eulerian)


def construct_dna_sequence(real_graph: nx.MultiDiGraph) -> str:
    """ This function takes a multidirected graph as input,
    and checks whether it is Eulerian or has a Euler's walk.
    If the graph doesn't have either, it returns a string that says
    the DNA sequence cannot be constructed. Otherwise, it returns the DNA
    sequence. """

    graph = nx.MultiDiGraph(real_graph)

    if is_valid_graph(graph):
        return _helper_construct_dna_sequence(graph)
    else:
        return 'DNA sequence can not be constructed.'


sequence = construct_dna_sequence(G)


# Task 8 Save DNA sequence (or save the error message if the sequence
# was not valid)

def save_output(s: str, full_output_path: str) -> None: # <-- MODIFIED
    with open(full_output_path, 'w') as file: # <-- MODIFIED
        file.write(s)


save_output(sequence, full_text_path) # <-- MODIFIED
