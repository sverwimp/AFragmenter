import json
from pathlib import Path
from typing import Union, Tuple
import os
import sys
from io import StringIO
import csv

from rich.console import Console
from rich.table import Table
import numpy as np
import igraph
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker import FormatStrFormatter
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.patches as patches
from matplotlib import image, axes

from .sequence_extractor import read_sequence


FilePath = Union[str, Path]


def validate_pae(pae: np.ndarray) -> None:
    """
    Validate some properties of the Predicted Aligned Error (PAE) matrix.

    Parameters:
    - pae (np.ndarray): The PAE matrix.

    Returns:
    - None

    Raises:
    - TypeError: If the PAE matrix is not a numpy array.
    - ValueError: If the PAE matrix is not 2D, not square, or if it contains negative values.
    """
    if not isinstance(pae, np.ndarray):
        raise TypeError("pae must be a numpy array")
    if pae.ndim != 2:
        raise ValueError("PAE matrix must be 2D")
    if pae.shape[0] != pae.shape[1]:
        raise ValueError("PAE matrix must be square")
    if np.min(pae) < 0:
        raise ValueError("PAE values must be non-negative")


def load_pae(json_file: FilePath) -> np.ndarray:
    """
    Load the Predicted Aligned Error (PAE) data from a JSON file.

    Parameters:
    - json_file (FilePath): The path to the JSON file. (str or Path)

    Returns:
    - np.ndarray: The PAE matrix.

    Raises:
    - FileNotFoundError: If the JSON file does not exist.
    - ValueError: If the PAE data is not found in the JSON file.
    """
    if not os.path.exists(json_file):
        raise FileNotFoundError(f"{json_file} not found")
    
    with open(json_file, "r") as f:
        pae_data = json.load(f)

    # AF2 format loads as a list containing a dictionary, AF3 and colabfold directly load the dictionary
    if isinstance(pae_data, list):
        pae_data = pae_data[0]

    # AFDB v1 and v2 have different keys for the PAE data
    if "distance" in pae_data:
        nrows = max(pae_data.get('residue1'))
        pae_matrix = np.zeros((nrows + 1, nrows + 1))
        for r, c, v in zip (pae_data.get('residue1'), pae_data.get('residue2'), pae_data.get('distance')):
            pae_matrix[r, c] = v
    else:
        pae = pae_data.get("predicted_aligned_error") or pae_data.get("pae")
        if pae is None:
            raise ValueError("PAE data not found in JSON file")
        pae_matrix = np.stack(pae, axis=0)

    validate_pae(pae_matrix)
    return pae_matrix
    

def _create_graph(weights_matrix: np.ndarray) -> igraph.Graph:
    """
    Create an igraph Graph object from a given weights matrix.

    Parameters:
    - weights_matrix (np.ndarray): A square matrix containing the edge weights.

    Returns:
    - igraph.Graph: The graph object.
    """
    num_vertices = weights_matrix.shape[0]
    indices = np.triu_indices(num_vertices, k=1)
    edges = list(zip(*indices))
    weights = [weights_matrix[i, j] for i, j in edges]
    
    graph = igraph.Graph()
    graph.add_vertices(num_vertices)
    graph.add_edges(edges)
    graph.es["weight"] = weights

    return graph


def cluster(graph: igraph.Graph, 
            resolution: Union[float, None] = None, 
            n_iterations: int = -1, 
            objective_function: str = "modularity") -> igraph.VertexClustering:
    #TODO: add some info about objective_function
    """
    Cluster the graph using the Leiden algorithm.
    
    The resolution parameter determines the scale of the partition cluster. It can be thought of as the coarseness of the clustering.
    A higher resolution will lead to more and smaller clusters, while a lower resolution will lead to fewer and larger clusters.

    The number of iterations determines the number of times the Leiden algorithm will be run. A negative value will run the algorithm until the partition stabilizes.

    Parameters:
    - graph (igraph.Graph): The graph object.
    - resolution (float, optional): The resolution parameter for the Leiden algorithm.
    - n_iterations (int, optional): The number of iterations for the Leiden algorithm.
    - objective_function (str, optional): The objective function for the Leiden algorithm.

    Returns:
    - igraph.VertexClustering: The resulting vertex clustering.
    """
    default_resolutions = {"modularity": 0.8, "cpm": 0.3}
    if objective_function.lower() not in default_resolutions:
        raise ValueError("Objective function must be either 'modularity' or 'CPM'")
    resolution = resolution or default_resolutions[objective_function.lower()]

    partition = graph.community_leiden(
        weights="weight",
        resolution=resolution, 
        n_iterations=n_iterations, 
        objective_function=objective_function
    )
    return partition


def plot_matrix(matrix: np.ndarray, 
                ax = None, 
                tick_top: bool = True, 
                colorbar: bool = True, 
                colorbar_label: str = "", 
                colorbar_decimals: int = None, 
                **kwargs) -> Tuple[image.AxesImage, axes.Axes]:
    """
    Plot a matrix as a heatmap.

    Parameters:
    - matrix (np.ndarray): The matrix to be plotted.
    - ax (axes.Axes, optional): The matplotlib axes object to plot the matrix on.
    - tick_top (bool, optional): Whether to place the ticks on the top of the plot.
    - colorbar (bool, optional): Whether to display a colorbar.
    - colorbar_label (str, optional): The label for the colorbar.
    - colorbar_decimals (int, optional): The number of decimal places to display in the colorbar.
    - **kwargs: Additional keyword arguments to be passed to the matplotlib.pyplot.imshow function.

    Returns:
    - Tuple[image.AxesImage, axes.Axes]: The image and axes objects.
    """
    paegreen = ["#1E461E", "#249224", "#37B137", "#56C956", "#82DD82", "#BAEFBA", "#FFFFFF"]
    cm = LinearSegmentedColormap.from_list("custom", paegreen, N=256)

    if np.all((matrix >= 0) & (matrix <= 1)):
        kwargs.setdefault("vmax", 1)
        if colorbar_decimals is None:
            colorbar_decimals = 1
    else:
        if colorbar_decimals is None:
            colorbar_decimals = 0       
    
    kwargs.setdefault("aspect", "equal")
    kwargs.setdefault("cmap", cm)
    kwargs.setdefault("vmin", 0)
    
    if ax is None:
        _, ax = plt.subplots()
    image = ax.imshow(matrix, **kwargs)

    if tick_top:
        ax.xaxis.tick_top()
    if colorbar:
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        cbar = plt.colorbar(image, cax=cax)
        cbar.set_label(colorbar_label, rotation=270, labelpad=15)
        cbar.ax.yaxis.set_major_formatter(FormatStrFormatter(f"%.{colorbar_decimals}f"))
    return image, ax


# TODO: Change the name of this function to something more descriptive / intuitive
def find_cluster_intervals(clusters: igraph.VertexClustering) -> dict:
    """
    Create a dictionary of cluster intervals from the vertex clustering object.

    Parameters:
    - clusters (igraph.VertexClustering): The vertex clustering object.

    Returns:
    - dict: A dictionary where the keys are the cluster indices and the values are lists of tuples representing the cluster intervals.

    Example output:
    {
        0: [(0, 20), (35, 40)],
        1: [(21, 34)]
    }
    """
    results = {}
    for i, cluster in enumerate(clusters):
        region = []
        if cluster:
            start = cluster[0]
            for j in range(1, len(cluster)):
                if cluster[j] != cluster[j-1] + 1:
                    region.append((start, cluster[j-1]))
                    start = cluster[j]
            region.append((start, cluster[-1]))
        results[i] = region
    return results


def filter_size_intervals(intervals, min_size):
    """
    Filter out intervals that are smaller than a given size.
    
    Parameters:
    - intervals (dict): A dictionary where the keys are the cluster indices and the values are lists of 
                        tuples representing the cluster intervals.
    - min_size (int): The minimum total size of the intervals to keep.

    Returns:
    - dict: A dictionary with the same structure as the input, but with intervals smaller than min_size
            removed and the keys renumbered starting from zero.
    """
    filtered_intervals = {}
    for key, value in intervals.items():
        total_length = sum(interval[1] - interval[0] + 1 for interval in value)
        if total_length >= min_size:
            filtered_intervals[key] = value
    
    # Remove empty entries
    filtered_intervals = {k: v for k, v in filtered_intervals.items() if v}
    
    # Renumber the keys starting from zero
    renumbered_intervals = {i: v for i, (k, v) in enumerate(filtered_intervals.items())}
    return renumbered_intervals


#def plot_chopping(background_matrix: np.ndarray, cluster_intervals: dict, ax = None, linewidth: int = 2, edgecolor = 'r', **kwargs):
def plot_result(background_matrix: np.ndarray, cluster_intervals: dict, ax = None, linewidth: int = 2, edgecolor = 'r', **kwargs):
    """
    Plot the cluster intervals as rectangles on top of a background matrix.

    Parameters:
    - background_matrix (np.ndarray): The background matrix to plot.
    - cluster_intervals (dict): A dictionary where the keys are the cluster indices and the values are lists of
                                tuples representing the cluster intervals.
    - ax (axes.Axes, optional): The matplotlib axes object to plot the matrix on.
    - linewidth (int, optional): The width of the rectangle edges.
    - edgecolor (str, optional): The color of the rectangle edges.
    - **kwargs: Additional keyword arguments to be passed to the matplotlib.pyplot.imshow function.

    Returns:
    - Tuple[image.AxesImage, axes.Axes]: The image and axes objects.
    """

    image, ax = plot_matrix(background_matrix, ax=ax, **kwargs)
    
    for intervals in cluster_intervals.values():
        # Plot diagonal square for each interval
        for start, stop in intervals:
            rect = patches.Rectangle((start, start), stop - start, stop - start, 
                                        linewidth=linewidth, edgecolor=edgecolor, facecolor='none')
            ax.add_patch(rect)
        # Plot off-diagonal rectangles for discontinuous intervals
        for i, (start1, stop1) in enumerate(intervals[:-1]):
            for start2, stop2 in intervals[i + 1:]:
                # Off-diagonal intersection rectangles
                width1, height1 = stop1 - start1, stop2 - start2
                width2, height2 = stop2 - start2, stop1 - start1
                
                # Add the first off-diagonal rectangle
                ax.add_patch(
                    patches.Rectangle((start1, start2), width1, height1,
                                    linewidth=linewidth, edgecolor=edgecolor, facecolor='none')
                )
                # Add the second off-diagonal rectangle (symmetric to the first)
                ax.add_patch(
                    patches.Rectangle((start2, start1), width2, height2,
                                    linewidth=linewidth, edgecolor=edgecolor, facecolor='none')
                )
    return image, ax


def _table_format(intervals: dict, **kwargs) -> str:
    """
    Print the cluster intervals in a table format.
    If the output is a terminal or a jupyter notebook, the table will be printed using rich.
    If the output is a file, the table will be printed using csv.

    Parameters:
    - intervals (dict): A dictionary where the keys are the cluster indices and the values are lists of tuples representing the cluster intervals.
    - **kwargs: Additional keyword arguments to be passed to the rich.table.Table constructor.

    Returns:
    - str: The formatted table as a string.
    """

    def is_notebook():
        try:
            # Check if the 'ipykernel' module is loaded
            if 'ipykernel' in sys.modules:
                return True
            # Check if the 'JPY_PARENT_PID' environment variable is set
            if 'JPY_PARENT_PID' in os.environ:
                return True
            return False
        except NameError:
            return False
    
    # If the output is a terminal or a jupyter notebook, use rich to print the table
    if sys.stdout.isatty() or is_notebook():
        table = Table(**kwargs)
        
        table.add_column("Domain", justify="right") # , style="magenta"
        table.add_column("Number of Residues", justify="right") # , style="green"
        table.add_column("Chopping", justify="right") # , style="yellow"

        for i, cluster in enumerate(intervals.values()):
            chopping = '_'.join([f'{start + 1}-{end + 1}' for start, end in cluster])
            nres = sum(end - start + 1 for start, end in cluster)
            table.add_row(str(i+1), str(nres), chopping)
        
        console = Console()
        with console.capture() as capture:
            console.print(table)
            console.quiet = True # This is a weird workaround to avoid rendering an empty output cell in jupyter notebooks
        console.quiet = False
        return capture.get().rstrip('\n')

    # If the output is a file, use csv to print the table (e.g. when redirecting the output to a file, for easier parsing)
    else:
        output = StringIO()
        writer = csv.writer(output, delimiter=',')
        writer.writerow(["Domain", "Number of Residues", "Chopping"])
        for i, cluster in enumerate(intervals.values()):
            chopping = '_'.join([f'{start + 1}-{end + 1}' for start, end in cluster])
            nres = sum(end - start + 1 for start, end in cluster)
            writer.writerow([str(i+1), str(nres), chopping])
        return output.getvalue().rstrip('\n')

    

    
"""
TODO:
- Implement visualization method using py3Dmol
- Implement method to save the results to a file
- Save parameters used to create the AFragmenter object in a class attribute (e.g. self.parameters) as a dictionary
- Implement a method to filter out fasta sequences based on the cluster intervals
- Implement a method to save the filtered fasta sequences to a file
"""


class AFragmenter:

    def __init__(self, pae_matrix: Union[np.ndarray, FilePath], threshold: float = 5.0):

        if isinstance(pae_matrix, (str, Path)):
            pae_matrix = load_pae(pae_matrix)
        self.pae_matrix = pae_matrix
        self.edge_weights_matrix = self._logistic_transform(pae_matrix, threshold)
        self.graph = _create_graph(self.edge_weights_matrix)
        

    def _logistic_transform(self, pae_matrix: np.ndarray, threshold: float) -> np.ndarray:
        if threshold < 0 or threshold > 31.75:
            raise ValueError("Threshold must be between 0 and 31.75")
        return 1 / (1 + np.exp(1.0 * (pae_matrix - threshold)))
    

    def cluster(self, resolution: Union[float, None] = None, n_iterations: int = -1, objective_function: str = "modularity", min_size: int = 0):
        if min_size < 0 or min_size > self.pae_matrix.shape[0]:
            raise ValueError("Minimum cluster size must be between 0 and the number of residues")

        clusters = cluster(self.graph, resolution=resolution, n_iterations=n_iterations, objective_function=objective_function)
        cluster_intervals = find_cluster_intervals(clusters)
        self.cluster_intervals = filter_size_intervals(cluster_intervals, min_size)
        return self
    

    # TODO: remove this function
    def print_results_nr(self, chain_id: str, sep="\t"):
        """
        Print the results of the FragmentFold analysis.

        Parameters:
        - chain_id (str): The ID of the chain.
        - sep (str, optional): The separator used between columns in the output. Defaults to "\t".

        Raises:
        - ValueError: If the cluster intervals are not defined.
        """
        if not hasattr(self, 'cluster_intervals'):
            raise ValueError("Intervals not defined. Run the cluster method first.")

        print("chain_id", "ndom", "nres", "chopping", sep=sep)
        for i, cluster in enumerate(self.cluster_intervals.values()):
            chopping = '_'.join([f'{start + 1}-{end + 1}' for start, end in cluster])
            nres = sum(end - start + 1 for start, end in cluster)
            print(chain_id, i+1, nres, chopping, sep=sep) # TODO: Change chain_id


    def plot_pae(self, **kwargs) -> Tuple[image.AxesImage, axes.Axes]:
        return plot_matrix(self.pae_matrix, **kwargs)


    def plot_results(self, **kwargs) -> Tuple[image.AxesImage, axes.Axes]:
        if not hasattr(self, "cluster_intervals"):
            raise ValueError("No clustering results found, please run the cluster method first")
        return plot_result(self.pae_matrix, self.cluster_intervals, **kwargs)
        

    def print_results(self, **kwargs) -> str:
        table_string = _table_format(self.cluster_intervals, **kwargs)
        print(table_string)


    def __str__(self) -> str:
        if not hasattr(self, "cluster_intervals"):
            raise ValueError("No clustering results found, please run the cluster method first")
        return _table_format(self.cluster_intervals)


    def visualize_py3Dmol(self, structure_file: str, color_range: list = None, add_surface: bool = False, surface_opacity: float = 0.7):
        """
        Visualize the 3D strucutre of the protein using py3Dmol. Color the residues based on the clusters.

        Parameters:
        - pdb_file (str): The path to the PDB file.
        - color_range (list, optional): A list of colors to use for the clusters. Defaults to None.
        - add_surface (bool, optional): Whether to add a surface to the structure. Defaults to False.
        - surface_opacity (float, optional): The opacity of the surface. Defaults to 0.7.

        Raises:
        - ImportError: If the py3Dmol library is not installed.
        - ValueError: If the cluster intervals are not defined.
        """
        try:
            import py3Dmol
        except ImportError:
            raise ImportError(
                "The py3Dmol library is required for this function. "
                "Please install it using 'pip install py3Dmol'."
            )
        
        if not hasattr(self, "cluster_intervals"):
            raise ValueError("No clustering results found, please run the cluster method first")
        
        if color_range is None:
            color_range = ['red', 'blue', 'green', 'yellow', 'orange', 'purple', 'cyan', 'magenta']

        # TODO: Do something about this
        if len(self.cluster_intervals) > len(color_range):
            print("Warning: More clusters than available colors. Some clusters will have the same color.")
        
        view = py3Dmol.view(js="https://3dmol.org/build/3Dmol.js")
        
        if structure_file.lower().endswith('pdb'):
            view.addModel(open(structure_file, 'r').read(), 'pdb')
        else:
            view.addModel(open(structure_file, 'r').read(), 'cif')

        view.setStyle({'cartoon': {'color': 'grey'}}) # spectrum
        for cluster, ranges in self.cluster_intervals.items():
            col = color_range[cluster % len(color_range)]
            for start, end in ranges:
                view.setStyle({'resi': f'{start}-{end}'}, {'cartoon': {'color': col}})
        if add_surface:
            view.addSurface(py3Dmol.VDW,{'opacity':surface_opacity,'color':'white'})

        view.zoomTo()
        view.show()


    def _clusters_to_fasta(self, sequence_file: FilePath, query_chain: str = 'A') -> dict:
        """
        Parse the sequence file and return the sequences corresponding to each cluster.
        """

        if not hasattr(self, 'cluster_intervals'):
            raise ValueError("Intervals not defined. Run the cluster method first.")
        
        sequence = read_sequence(sequence_file, seq_length=self.pae_matrix.shape[0] , query_chain=query_chain)
        parsed_sequences = {}
        for cluster_id, interval in self.cluster_intervals.items():
            cluster_sequence = []
            for i, (start, end) in enumerate(interval):
                if i > 0: # If it isn't the first interval, check for gaps
                    previous_end = interval[i-1][1]
                    if start > previous_end:
                        gap_size = start - previous_end
                        cluster_sequence.extend(['-'] * gap_size)
                
                cluster_sequence.extend(sequence[start:end+1])
            parsed_sequences[cluster_id] = ''.join(cluster_sequence)
        return parsed_sequences


    def print_fasta(self, sequence_file: FilePath, prefix: str = "", width: int = 60):
        parsed_sequences = self._clusters_to_fasta(sequence_file)
        
        if not prefix:
            prefix = Path(sequence_file).stem

        for i, seq in parsed_sequences.items():
            interval = self.cluster_intervals[i]
            interval_str = "_".join([f"{start}-{end}" for start, end in interval])
            print(f">{prefix}_{i+1} {interval_str}")
            for i in range(0, len(seq), width):
                print(seq[i:i+width])


    def save_fasta(self, sequence_file: FilePath, output_file: FilePath, prefix: str = "", width: int = 60):
        parsed_sequences = self._clusters_to_fasta(sequence_file)
        
        if not prefix:
            prefix = Path(sequence_file).stem

        output_file = Path(output_file)
        
        with open(output_file, 'w') as f:
            for i, seq in parsed_sequences.items():
                interval = self.cluster_intervals[i]
                interval_str = "_".join([f"{start}-{end}" for start, end in interval])
                f.write(f">{prefix}_{i+1} {interval_str}\n")
                for i in range(0, len(seq), width):
                    f.write(seq[i:i+width] + "\n")
