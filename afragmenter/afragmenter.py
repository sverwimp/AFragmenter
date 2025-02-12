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
from matplotlib import image, axes

from .data_reader import PAEHandler, SequenceReader
from afragmenter import plotting


FilePath = Union[str, Path]
default_resolutions = {"modularity": 0.8, "cpm": 0.3}


def create_graph(weights_matrix: np.ndarray) -> igraph.Graph:
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
    
    graph = igraph.Graph(n=num_vertices, edges=edges, edge_attrs={'weight': weights}, directed=False)
    #graph.add_vertices(num_vertices)
    #graph.add_edges(edges)
    #graph.es["weight"] = weights

    return graph


def cluster_graph(graph: igraph.Graph,
            resolution: Union[float, None] = None, 
            n_iterations: int = -1, 
            objective_function: str = "modularity",
            **kwargs) -> igraph.VertexClustering:
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

    objective_funtions = ["modularity", "cpm"]
    if objective_function.lower() not in objective_funtions:
        raise ValueError(f"Objective function must be one of {objective_funtions}, got {objective_function}")

    if resolution is None:
        resolution = default_resolutions[objective_function.lower]

    partition = graph.community_leiden(
        weights="weight",
        resolution=resolution, 
        n_iterations=n_iterations, 
        objective_function=objective_function,
        **kwargs
    )
    return partition


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


def remove_min_size_intervals(intervals: dict, min_size: int) -> dict:
    """
    Remove clusters that are smaller than a given size based on the *total* size of the intervals (values) per group (keys).

    In example 1, clusters 0, 2 and 4 are removed because the intervals are too small.
    In example 2, clusters 0 and 3 are removed because the intervals are too small. Cluser 1 has intervals that
    are smaller than the min_size threshold, but the total size of the intervals is larger than min_size, so it is kept.

    Example 1:
        intervals = {0: [(0, 2)],                      -> total size interval = 3
                    1: [(3, 33)],                      -> total size interval = 31
                    2: [(34, 35)], 
                    3: [(36, 133), (143, 269)],        -> total size interval = 225
                    4: [(134, 138)]}
        min_size = 10
        merge_intervals(intervals, min_size) -> {1: [(3, 33)], 
                                                 3: [(36, 133), (143, 269)]}
    
    Example 2:
        intervals = {0: [(0, 5)],                                               -> total size interval = 6
                    1: [(6, 11), (807, 807), (810, 811), (813, 933)],           -> total size interval = 130
                    2: [(12, 179), (419, 518), (519, 603), (604, 806), 
                        (808, 809), (812, 812)],
                    3: [(180, 182)], 
                    4: [(183, 229), (230, 418)]}
        min_size = 10
        merge_intervals(intervals, min_size) -> {1: [(6, 11), (807, 807), (810, 811), (813, 933)],
                                                 2: [(12, 179), (419, 518), (519, 603), (604, 806), 
                                                     (808, 809), (812, 812)],
                                                 4: [(183, 229), (230, 418)]}

    Parameters:
    - intervals (dict): A dictionary where the keys are the cluster indices and the values are lists of 
                        tuples representing the cluster intervals.
    - min_size (int): The minimum total size of the intervals to keep.

    Returns:
    - dict: A dictionary with the same structure as the input, but with intervals smaller than min_size removed.
    """
    filtered_intervals = {}
    for key, value in intervals.items():
        total_size = sum(interval[1] - interval[0] + 1 for interval in value)
        if total_size >= min_size:
            filtered_intervals[key] = value

    # Remove any empty entries
    filtered_intervals = {k: v for k, v in filtered_intervals.items() if v}
    
    return filtered_intervals


def merge_intervals(intervals: dict, min_size: int) -> dict:
    """
    Iterates over all intervals and removes those that are smaller than the min_size threshold regardless of 
    the total size of the intervals in a cluster/group (only the size of the individual intervals is considered).
    Merges adjacent intervals of the same cluster while skipping intervals below the min_size threshold.

    Example:
        intervals = {0: [(0, 2)], 1: [(3, 33)], 2: [(34, 35)], 3: [(36, 133), (143, 269)], 4: [(134, 138)], 
                    5: [(139, 142)], 6: [(270, 272)]}
        min_size = 10
        merge_intervals(intervals, min_size) -> {1: [(3, 33)], 3: [(36, 269)]}

    In the example, cluster 0, 2 and 4 are removed because the intervals are too small. The intervals of cluster 3 
    are merged because there is no other cluster with a size >= min_size between its intervals.

    Parameters:
    - intervals (dict): A dictionary where the keys are the cluster indices and the values are lists of 
                        tuples representing the cluster intervals.
    - min_size (int): The minimum total size of the intervals to keep.

    Returns:
    - dict: A dictionary with the same structure as the input, but with intervals smaller than min_size 
            removed and the keys renumbered starting
    """
    class Interval:
        def __init__(self, group: int, start: int, end: int):
            self.group = group
            self.start = start
            self.end = end
            self.size = end - start + 1
        
        def __repr__(self):
            return f"Interval(group={self.group}, start={self.start}, end={self.end}, size={self.size})"

    interval_objects = []
    for key, value in intervals.items():
        for results in value:
            start, end = results
            interval_objects.append(Interval(key, start, end))
    interval_objects.sort(key=lambda x: x.start)

    previous_valid = None
    results = []
    for current in interval_objects:
        # Skip intervals that are too small
        if current.size < min_size:
            continue

        # First interval above the min_size threshold
        if previous_valid is None:
            previous_valid = current
            continue

        # If the current interval is adjacent to the previous one and of the same group, merge them
        if previous_valid.group == current.group:
            previous_valid = Interval(previous_valid.group, previous_valid.start, current.end)
        else:
            results.append(previous_valid)
            previous_valid = current
        
    # Add the last interval
    if previous_valid is not None:
        results.append(previous_valid)

    # Group the merged intervals by cluster
    merged_intervals = {}
    for interval in results:
        if interval.group not in merged_intervals:
            merged_intervals[interval.group] = []
        merged_intervals[interval.group].append((interval.start, interval.end))
    
    return merged_intervals


def filter_cluster_intervals(intervals: dict, min_size: int, attempt_merge: bool = True) -> dict:
    """
    Filter out intervals that are smaller than a given size. 
    Optionally attempt to merge smaller clusters with adjacent larger ones.

    Parameters:
    - intervals (dict): A dictionary where the keys are the cluster indices and the values are lists of 
                        tuples representing the cluster intervals.
    
    - min_size (int): The minimum total size of the intervals to keep.
    - attempt_merge (bool, optional): Whether to try to merge smaller clusters with adjacent larger ones. 
                                      Defaults to True.

    Returns:
    - dict: A dictionary with the same structure as the input, but with intervals smaller than min_size
            removed and the keys renumbered starting from zero.
    """

    if min_size < 0:
        raise ValueError("Minimum cluster size must be greater than or equal to 0")
    if min_size == 0:
        return intervals

    if attempt_merge:
        intervals = merge_intervals(intervals, min_size)
    else:
        intervals = remove_min_size_intervals(intervals, min_size)

    # Renumber the keys starting from zero
    renumbered_intervals = {i: v for i, (_, v) in enumerate(intervals.items())}    
    return renumbered_intervals


def format_results_table(intervals: dict, **kwargs) -> str:
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
    def is_notebook() -> bool:
        try:
            from IPython import get_ipython
            if get_ipython() is not None:
                return True
        except ImportError:
            return False
        return False
    
    def format_as_rich_table(intervals: dict, **kwargs) -> str:
        table = Table(**kwargs)
        table.add_column("Domain", justify="right")
        table.add_column("Number of Residues", justify="right")
        table.add_column("Chopping", justify="right")

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
    
    def format_as_csv(intervals: dict) -> str:
        output = StringIO()
        writer = csv.writer(output, delimiter=',')
        writer.writerow(["domain", "nres", "chopping"])
        for i, cluster in enumerate(intervals.values()):
            chopping = '_'.join([f'{start + 1}-{end + 1}' for start, end in cluster])
            nres = sum(end - start + 1 for start, end in cluster)
            writer.writerow([str(i+1), str(nres), chopping])
        return output.getvalue().rstrip('\n')

    if sys.stdout.isatty() or is_notebook():
        return format_as_rich_table(intervals, **kwargs)
    else:
        return format_as_csv(intervals)


class AFragmenter:
    """
    AFragmenter class for clustering protein domains based on Predicted Aligned Error (PAE) values.

    Parameters:
    - pae_matrix (Union[np.ndarray, FilePath, list, dict, StringIO]): The Predicted Aligned Error matrix.
    - threshold (float, optional): The threshold for the sigmoid function used to transform the PAE values into graph edge weights.

    Attributes:
    - pae_matrix (np.ndarray): The Predicted Aligned Error matrix.
    - edge_weights_matrix (np.ndarray): The matrix of edge weights. This is created by transforming the PAE matrix using a 
                                        logistic function to increase the contrast between high and low PAE values.
    - graph (igraph.Graph): The graph object created from the edge_weights_matrix.

    Methods:
    - cluster: Cluster the graph using the Leiden algorithm.
    - plot_pae: Plot the Predicted Aligned Error matrix as a heatmap.
    - plot_result: Plot the clustering results on top of the Predicted Aligned Error matrix.
    - print_result: Print the clustering results in a table format.
    - visualize_py3Dmol: Visualize the 3D structure of the protein using py3Dmol. (Requires the py3Dmol library to be installed)
    - print_fasta: Print the sequences corresponding to each cluster in FASTA format.
    - save_fasta: Save the sequences corresponding to each cluster in FASTA format to a file.
    """

    def __init__(self, pae_matrix: Union[np.ndarray, FilePath, list, dict, StringIO], threshold: float = 5.0):
        if isinstance(pae_matrix, (list, dict, StringIO)):
            pae_matrix = PAEHandler.process_pae_data(pae_matrix)
        elif isinstance(pae_matrix, (str, Path)):
            pae_matrix = PAEHandler.load_pae(pae_matrix)
        elif not isinstance(pae_matrix, np.ndarray):
            raise TypeError("pae_matrix must be a numpy array, a file path, or a dictionary containing PAE data")
        
        self.pae_matrix = pae_matrix
        self.edge_weights_matrix = AFragmenter._pae_transform(pae_matrix, threshold)
        self.graph = create_graph(self.edge_weights_matrix)
        self.params = {"threshold": threshold}
        self.sequence_reader = None
        

    @staticmethod
    def _pae_transform(pae_matrix: np.ndarray, threshold: float) -> np.ndarray:
        """
        Transform the PAE matrix into a matrix of edge weights using a sigmoid function.
        This creates a larger contrast between the high and low PAE values.

        Parameters:
        - pae_matrix (np.ndarray): The Predicted Aligned Error matrix.
        - threshold (float): The threshold for the sigmoid function.

        Returns:
        - np.ndarray: The edge weights matrix.

        Raises:
        - ValueError: If the threshold is less than 0 or greater than 31.75.
        """

        if threshold < 0 or threshold > 31.75:
            raise ValueError("Threshold must be between 0 and 31.75")
        return 1 / (1 + np.exp(1.0 * (pae_matrix - threshold)))
    

    def cluster(self, resolution: Union[float, None] = None, 
                objective_function: str = "modularity", 
                n_iterations: int = -1, 
                min_size: int = 10,
                attempt_merge: bool = True) -> 'AFragmenter':
        """
        Create a graph from the edge_weights_matrix and cluster it using the Leiden algorithm.

        Parameters:
        - resolution (float, optional): The resolution parameter for the Leiden algorithm. Recommended to be between 0.0 and 1.0
                                        A higher resolution will lead to more and smaller clusters, while a lower resolution will lead to fewer and larger clusters.
        - objective_function (str, optional): The objective function for the Leiden algorithm. [modularity, CPM]
        - n_iterations (int, optional): The number of iterations for the Leiden algorithm. 
                                        If a negative value is given, the algorithm will run until a stable iteration is reached.
        - min_size (int, optional): The minimum size of the clusters to keep. Must be between 0 and the number of residues.
        - attempt_merge (bool, optional): Whether to attempt to merge smaller clusters with adjacent larger ones. Defaults to True.

        Returns:
        - AFragmenter: The AFragmenter object with the cluster intervals stored in the cluster_intervals attribute.

        Raises:
        - ValueError: If the minimum cluster size is less than 0 or greater than the number of residues.
        """

        if min_size < 0 or min_size > self.pae_matrix.shape[0]:
            raise ValueError(f"Minimum cluster size must be between 0 and the number of residues in the protein ({self.pae_matrix.shape[0]})")
        
        objective_functions = ["modularity", "cpm"]
        if objective_function.lower() not in objective_functions:
            raise ValueError(f"Objective function must be one of {objective_functions}, got {objective_function}")
        
        if resolution is None:
            resolution = default_resolutions[objective_function.lower()]

        self.params.update({
            "resolution": resolution,
            "objective_function": objective_function,
            "n_iterations": n_iterations,
            "min_size": min_size,
            "attempt_merge": attempt_merge
        })
        
        clusters = cluster_graph(self.graph, resolution=resolution, n_iterations=n_iterations, objective_function=objective_function)
        cluster_intervals = find_cluster_intervals(clusters)
        self.cluster_intervals = filter_cluster_intervals(cluster_intervals, min_size, attempt_merge=attempt_merge)
        return self


    def plot_pae(self, **kwargs) -> Tuple[image.AxesImage, axes.Axes]:
        """
        Plot the Predicted Aligned Error matrix as a heatmap.

        Parameters:
        - **kwargs: Additional keyword arguments to be passed to the matplotlib.pyplot.imshow function.

        Returns:
        - Tuple[image.AxesImage, axes.Axes]: The image and axes objects.
        """
        kwargs.setdefault("colorbar_label", "Predicted Aligned Error (Å)")
        image, ax = plotting.plot_matrix(self.pae_matrix, **kwargs)
        ax.set_xlabel("Scored residue")
        ax.set_ylabel("Aligned residue")

        return image, ax


    def plot_result(self, **kwargs) -> Tuple[image.AxesImage, axes.Axes]:
        """
        Plot the clustering results on top of the Predicted Aligned Error matrix.

        Parameters:
        - **kwargs: Additional keyword arguments to be passed to the matplotlib.pyplot.imshow function.

        Returns:
        - Tuple[image.AxesImage, axes.Axes]: The image and axes objects.
        """
        if not hasattr(self, "cluster_intervals"):
            raise ValueError("No clustering results found, please run the cluster method first")
        
        kwargs.setdefault("colorbar_label", "Predicted Aligned Error (Å)")
        image, ax = plotting.plot_cluster_intervals(self.pae_matrix, self.cluster_intervals, **kwargs)
        ax.set_xlabel("Scored residue")
        ax.set_ylabel("Aligned residue")

        return image, ax
        

    def print_result(self, **kwargs) -> None:
        """
        Print the clustering results in a table format. Will use rich if the output is a terminal or a jupyter notebook, else will use csv.

        Parameters:
        - **kwargs: Additional keyword arguments to be passed to the rich.table.Table constructor.

        Returns:
        - None

        Raises:
        - ValueError: If the cluster intervals are not defined.
        """
        if not hasattr(self, "cluster_intervals"):
            raise ValueError("No clustering results found, please run the cluster method first")
        table_string = format_results_table(self.cluster_intervals, **kwargs)
        print(table_string)


    def __repr__(self) -> str:
        return (f"AFragmenter(pae_matrix=np.ndarray(shape={self.pae_matrix.shape}), "
                f"params={self.params} "
                f"cluster_intervals={self.cluster_intervals if hasattr(self, 'cluster_intervals') else None})")


    def __str__(self) -> str:
        return (f"AFragmenter(pae_matrix=np.ndarray(shape={self.pae_matrix.shape}), "
                f"params={self.params} "
                f"cluster_intervals={self.cluster_intervals if hasattr(self, 'cluster_intervals') else None})")


    def visualize_py3Dmol(self, structure_file: str, 
                          color_range: list = ['red', 'blue', 'green', 'yellow', 'purple', 'orange', 'cyan', 'magenta', 
                                               'lime', 'pink', 'teal', 'lavender', 'brown', 'apricot', 'maroon', 'mint', 'olive', 
                                               'beige', 'navy', 'grey', 'white', 'black'],
                          size: Tuple[int, int] = (800, 600),
                          style: str = 'cartoon',
                          add_surface: bool = False, 
                          surface_opacity: float = 0.7) -> None:
        """
        Visualize the 3D strucutre of the protein using py3Dmol. Color the residues based on the clusters.

        Parameters:
        - pdb_file (str): The path to the PDB or mmcif file.
        - color_range (list, optional): A list of colors to use for the clusters, expects color names or hex-codes.
        - size (Tuple[int, int], optional): The width and height of the viewer. Defaults to (800, 600).
        - style (str, optional): The style to use. Defaults to 'cartoon'.
        - add_surface (bool, optional): Whether to add a surface to the structure. Defaults to False.
        - surface_opacity (float, optional): The opacity of the surface. Defaults to 0.7.

        Raises:
        - ImportError: If the py3Dmol library is not installed.
        - ValueError: If the cluster intervals are not defined.
        - ValueError: If the structure file is not a PDB or mmCIF file.
        """
        try:
            import py3Dmol
        except ImportError:
            raise ImportError(
                "The py3Dmol library is required for the visualize_py3Dmol function. "
                "Please install it using 'pip install py3Dmol'."
            )
        
        if not hasattr(self, "cluster_intervals"):
            raise ValueError("No clustering results found, please run the cluster method first")

        # TODO: Do something about this
        if len(self.cluster_intervals) > len(color_range):
            print("Warning: More clusters than available colors. Some clusters will have the same color.")
    
        view = py3Dmol.view(width=size[0], height=size[1])

        if os.path.isfile(structure_file):
            content = open(structure_file, 'r').read()
        else:
            content = structure_file
        
        file_format = SequenceReader.determine_file_format(content).lower()
        if file_format == 'pdb':
            view.addModel(content, 'pdb')
        elif file_format == 'mmcif':
            view.addModel(content, 'cif')
        else:
            raise ValueError("Unsupported file format. Please provide a PDB or mmCIF file.")

        view.setStyle({f'{style}': {'color': 'grey'}}) # spectrum
        for cluster, ranges in self.cluster_intervals.items():
            col = color_range[cluster % len(color_range)]
            for start, end in ranges:
                view.setStyle({'resi': f'{start+1}-{end+1}'}, {f'{style}': {'color': col}})
        if add_surface:
            view.addSurface(py3Dmol.VDW,{'opacity':surface_opacity,'color':'white'})

        view.zoomTo()
        view.show()


    def print_fasta(self, sequence_file: FilePath, prefix: str = "", width: int = 60) -> None:
        """
        Print the sequences corresponding to each cluster in FASTA format.

        Parameters:
        - sequence_file (FilePath): The path to the sequence file.
        - prefix (str, optional): The prefix to add to the sequence headers. Defaults to self.sequence_reader.name.
        - width (int, optional): The width of the sequence lines. Defaults to 60.

        Raises:
        - ValueError: If the cluster intervals are not defined.
        """
        if self.sequence_reader is None:
            self.sequence_reader = SequenceReader(sequence_file, seq_length=self.pae_matrix.shape[0])
        parsed_sequences = self.sequence_reader.clusters_to_fasta(self.cluster_intervals)
        
        if not prefix:
            prefix = self.sequence_reader.name

        for i, seq in parsed_sequences.items():
            interval = self.cluster_intervals[i]
            interval_str = "_".join([f"{start+1}-{end+1}" for start, end in interval])
            print(f">{prefix}_{i+1} {interval_str}")
            for i in range(0, len(seq), width):
                print(seq[i:i+width])


    def save_fasta(self, sequence_file: FilePath, output_file: FilePath, prefix: str = "", width: int = 60) -> None:
        """
        Save the sequences corresponding to each cluster in FASTA format to a file.

        Parameters:
        - sequence_file (FilePath): The path to the sequence file.
        - output_file (FilePath): The path to save the output file.
        - prefix (str, optional): The prefix to add to the sequence headers. Defaults to self.sequence_reader.name.
        - width (int, optional): The width of the sequence lines. Defaults to 60.

        Raises:
        - ValueError: If the cluster intervals are not defined.
        """
        if self.sequence_reader is None:
            self.sequence_reader = SequenceReader(sequence_file, seq_length=self.pae_matrix.shape[0])
        parsed_sequences = self.sequence_reader.clusters_to_fasta(self.cluster_intervals)
        
        if not prefix:
            prefix = self.sequence_reader.name

        output_file = Path(output_file)
        
        with open(output_file, 'w') as f:
            for i, seq in parsed_sequences.items():
                interval = self.cluster_intervals[i]
                interval_str = "_".join([f"{start+1}-{end+1}" for start, end in interval])
                f.write(f">{prefix}_{i+1} {interval_str}\n")
                for i in range(0, len(seq), width):
                    f.write(seq[i:i+width] + "\n")
