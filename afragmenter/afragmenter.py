from pathlib import Path
from typing import Union, Tuple
import os
import sys
from io import StringIO
import csv

from rich.console import Console
from rich.table import Table
import numpy as np
from matplotlib import image, axes

from .pae_handler import PAEHandler
from .sequence_reader import SequenceReader
from .graph import create_graph, cluster_graph, default_resolutions
from .intervals import find_cluster_intervals, filter_cluster_intervals
from afragmenter import format_result
from afragmenter import plotting


FilePath = Union[str, Path]


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
    

    def run(self, resolution: Union[float, None] = None, 
            objective_function: str = "modularity", 
            n_iterations: int = -1, 
            min_size: int = 10,
            attempt_merge: bool = True) -> 'AFragmenter':
        """Alias for the cluster method."""
        return self.cluster(resolution=resolution, 
                            objective_function=objective_function, 
                            n_iterations=n_iterations, 
                            min_size=min_size, 
                            attempt_merge=attempt_merge)


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


    def __repr__(self) -> str:
        return (f"AFragmenter(pae_matrix=np.ndarray(shape={self.pae_matrix.shape}), "
                f"params={self.params} "
                f"cluster_intervals={self.cluster_intervals if hasattr(self, 'cluster_intervals') else None})")


    def __str__(self) -> str:
        return (f"AFragmenter(pae_matrix=np.ndarray(shape={self.pae_matrix.shape}), "
                f"params={self.params} "
                f"cluster_intervals={self.cluster_intervals if hasattr(self, 'cluster_intervals') else None})")
    

    def print_result(self, format: str = 'auto', delimiter: str = ',') -> None:
        """
        Print the clustering results in a table format. Will use rich if the output is a terminal or a jupyter notebook, else will use csv.

        Parameters:
        - format (str, optional): The format to use. [auto, csv, rich]
        - delimiter (str, optional): The delimiter to use for the csv format.

        Returns:
        - None

        Raises:
        - ValueError: If the cluster intervals are not defined.
        """
        if not hasattr(self, "cluster_intervals"):
            raise ValueError("No clustering results found, please run the cluster method first")
        format_result.print_result(self.cluster_intervals, format=format, delimiter=delimiter)

    
    def save_result(self, output_file: FilePath, format: str = 'csv', delimiter: str = ',') -> None:
        """
        Save the clustering results in a table format to a file.

        Parameters:
        - output_file (FilePath): The path to save the output file.
        - format (str, optional): The format to use. [auto, csv, rich]
        - delimiter (str, optional): The delimiter to use for the csv format.

        Returns:
        - None

        Raises:
        - ValueError: If the cluster intervals are not defined.
        """
        if not hasattr(self, "cluster_intervals"):
            raise ValueError("No clustering results found, please run the cluster method first")
        format_result.save_result(self.cluster_intervals, output_file, format=format, delimiter=delimiter)
        

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
            import py3Dmol # type: ignore
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


    def format_fasta_sequences(self, parsed_sequences, prefix, width):
        """
        Generator function to format sequences in FASTA format.

        Parameters:
        - parsed_sequences (dict): The parsed sequences.
        - prefix (str): The prefix to add to the sequence headers.
        - width (int): The width of the sequence lines.

        Yields:
        - str: The formatted FASTA sequence.
        """
        for i, seq in parsed_sequences.items():
            interval = self.cluster_intervals[i]
            interval_str = "_".join([f"{start+1}-{end+1}" for start, end in interval])
            yield f">{prefix}_{i+1} {interval_str}"
            for j in range(0, len(seq), width):
                yield seq[j:j+width]


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

        for line in self.format_fasta_sequences(parsed_sequences, prefix, width):
            print(line)


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

        with output_file.open('w') as f:
            for line in self.format_fasta_sequences(parsed_sequences, prefix, width):
                f.write(line + "\n")