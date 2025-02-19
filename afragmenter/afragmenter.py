from pathlib import Path
from typing import Optional, Union, Tuple
import os
from io import StringIO

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
    - params (dict): The parameters used for clustering.
    - cluster_intervals (dict): The cluster intervals obtained from the clustering results. Each key represents a cluster number,
                                and the value is a list of tuples containing the start and end indices of the cluster intervals.
    - sequence_reader (SequenceReader): The SequenceReader object used to read the sequence file. This is automatically created
                                        when the print_fasta or save_fasta methods are called.

    Methods:
    - cluster: Cluster the graph using the Leiden algorithm.
    - run: Alias for the cluster method.
    - plot_pae: Plot the Predicted Aligned Error matrix as a heatmap.
    - plot_result: Plot the clustering results on top of the Predicted Aligned Error matrix.
    - print_result: Print the clustering results in a csv or table format.
    - save_result: Save the clustering results in a csv or table format to a file.
    - py3Dmol: Visualize the 3D structure of the protein using py3Dmol. (Requires the py3Dmol library to be installed)
    - print_fasta: Print the sequences corresponding to each cluster in FASTA format.
    - save_fasta: Save the sequences corresponding to each cluster in FASTA format to a file.
    """

    def __init__(self, pae_matrix: Union[np.ndarray, FilePath, list, dict, StringIO], threshold: float = 5.0):
        if isinstance(pae_matrix, (list, dict)):
            pae_matrix = PAEHandler.process_pae_data(pae_matrix)
        elif isinstance(pae_matrix, (str, Path, StringIO)):
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
    

    def cluster(self, resolution: Optional[float] = None, 
                objective_function: str = "modularity", 
                n_iterations: int = -1, 
                min_size: int = 10,
                attempt_merge: bool = True,
                **kwargs) -> 'AFragmenter':
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
        - **kwargs: Additional keyword arguments to be passed to the community_leiden function from igraph.

        Returns:
        - AFragmenter: The AFragmenter object with the cluster intervals stored in the cluster_intervals attribute.

        Raises:
        - ValueError: If the minimum cluster size is less than 0 or greater than the number of residues.
        """

        if min_size < 0 or min_size > self.pae_matrix.shape[0]:
            raise ValueError(f"Minimum cluster size must be between 0 and the number of residues in the protein ({self.pae_matrix.shape[0]})")
    
        # Update the parameters with the input values in case of error during clustering
        # Overwritten by default values, if applicable, after clustering
        self.params.update({ 
            "resolution": resolution,
            "objective_function": objective_function,
            "n_iterations": n_iterations,
            "min_size": min_size,
            "attempt_merge": attempt_merge
        })
        
        clusters, params = cluster_graph(graph=self.graph, 
                                 resolution=resolution, 
                                 n_iterations=n_iterations,
                                 objective_function=objective_function,
                                 return_params=True,
                                 **kwargs)
        self.params.update(params) # Update the parameters with the actual values used for clustering
        cluster_intervals = find_cluster_intervals(clusters)
        self.cluster_intervals = filter_cluster_intervals(intervals=cluster_intervals, 
                                                          min_size=min_size, 
                                                          attempt_merge=attempt_merge)
        return self
    

    def run(self, resolution: Optional[float] = None, 
            objective_function: str = "modularity", 
            n_iterations: int = -1, 
            min_size: int = 10,
            attempt_merge: bool = True,
            **kwargs) -> 'AFragmenter':
        """Alias for the cluster method."""
        return self.cluster(resolution=resolution, 
                            objective_function=objective_function, 
                            n_iterations=n_iterations, 
                            min_size=min_size, 
                            attempt_merge=attempt_merge,
                            **kwargs)
    

    def check_has_cluster_intervals(self) -> None:
        """
        Check if the cluster_intervals attribute is defined.
        If no results are found, raise a ValueError.

        Returns:
        - None

        Raises:
        - ValueError: If the cluster_intervals attribute is not defined.
        """
        if not hasattr(self, "cluster_intervals"):
            raise ValueError("No clustering results found, please run the cluster method first")


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
        self.check_has_cluster_intervals()
        
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
    

    def print_result(self, base_name: Optional[str] = None, format: str = 'auto', delimiter: str = ',') -> None:
        """
        Print the clustering results in a table format. Will use rich if the output is a terminal or a jupyter notebook, else will use csv.

        Parameters:
        - base_name (str, optional): the base name to use for the cluster intervals.
        - format (str, optional): The format to use. [auto, csv, rich]
        - delimiter (str, optional): The delimiter to use for the csv format.

        Returns:
        - None

        Raises:
        - ValueError: If the cluster intervals are not defined.
        """
        self.check_has_cluster_intervals()
        
        # not 'if not base_name' because base_name can be an empty string such that only the cluster number is used as name
        if base_name is None:
            base_name = self.sequence_reader.name if self.sequence_reader else None

        format_result.print_result(self.cluster_intervals, format=format, delimiter=delimiter, base_name=base_name)

    
    def save_result(self, output_file: FilePath, base_name: Optional[str] = None, format: str = 'csv', delimiter: str = ',') -> None:
        """
        Save the clustering results in a table format to a file.

        Parameters:
        - output_file (FilePath): The path to save the output file.
        - base_name (str, optional): the base name to use for the cluster intervals.
        - format (str, optional): The format to use. [auto, csv, rich]
        - delimiter (str, optional): The delimiter to use for the csv format.

        Returns:
        - None

        Raises:
        - ValueError: If the cluster intervals are not defined.
        """
        self.check_has_cluster_intervals()
        
        if base_name is None:
            base_name = self.sequence_reader.name if self.sequence_reader else None

        format_result.save_result(self.cluster_intervals, output_file, base_name=base_name, format=format, delimiter=delimiter)
    
    
    def py3Dmol(self, structure_file: str,
                color_palette: Optional[list] = None,
                style: str = 'cartoon',
                displace_domains: bool = False,
                tether_strength: float = 0.1,
                repulse_strength: float = 100,
                steps: int = 1000,
                dt: float = 0.01,
                **kwargs) -> 'py3Dmol.view':
        """
        Visualize the 3D structure of the protein using py3Dmol. Color the residues based on the clusters.

        Parameters:
        - structure_file (str): The path to the PDB or mmcif file.
        - color_palette (list, optional): A list of colors to use for the clusters, expects color names or hex-codes.
        - style (str, optional): The style to use. Defaults to 'cartoon'.
        - displace_domains (bool, optional): Whether to displace the domains based on the clusters. Defaults to False.
        - tether_strength (float, optional): The spring constant pulling each group toward its original center. Defaults 
                                             to 0.1. Only used if displace_domains is True.
        - repulse_strength (float, optional): The constant for the inverse-square repulsion between groups. Defaults to 100.
                                              Only used if displace_domains is True.
        - steps (int, optional): The number of simulation iterations. Defaults to 1000. Only used if displace_domains is True.
        - dt (float, optional): The time step for the simulation. Defaults to 0.01. Only used if displace_domains is True.
        - **kwargs: Additional keyword arguments to be passed to the py3Dmol.view function.

        Returns:
        - py3Dmol.view: The py3Dmol viewer object.

        Raises:
        - ImportError: If the py3Dmol library is not installed.
        - ValueError: If the cluster intervals are not defined.
        - ValueError: If the structure file is not a PDB or mmCIF file.
        """
        
        try:
            # Check if py3Dmol is installed for view_py3Dmol function
            import py3Dmol # type: ignore 
            from .structure_viewer import view_py3Dmol
        except ImportError:
            raise ImportError(
                "The py3Dmol library is required for the py3Dmol function. "
                "Please install it using 'pip install py3Dmol'."
            )
        
        self.check_has_cluster_intervals()

        if color_palette is None:
            from .structure_viewer import COLOR_PALETTE
            color_palette = COLOR_PALETTE

        view = view_py3Dmol(structure_file,
                            self.cluster_intervals,
                            displace_domains=displace_domains,
                            color_palette=color_palette,
                            style=style,
                            tether_strength=tether_strength,
                            repulse_strength=repulse_strength,
                            steps=steps,
                            dt=dt,
                            **kwargs)
        return view


    def _format_fasta_sequences(self, parsed_sequences: dict, header_name: str, width: int) -> str: # type: ignore (generator)
        """
        Generator function to format sequences in FASTA format.

        Parameters:
        - parsed_sequences (dict): The parsed sequences.
        - prefix (str): The name to add to the sequence headers.
        - width (int): The width of the sequence lines.

        Yields:
        - str: The formatted FASTA sequence.
        """
        for i, seq in parsed_sequences.items():
            interval = self.cluster_intervals[i]
            interval_str = "_".join([f"{start+1}-{end+1}" for start, end in interval])
            name = f"{header_name}_{i+1}" if header_name else f"{i+1}"
            yield f">{name} {interval_str}"
            for j in range(0, len(seq), width):
                yield seq[j:j+width]


    def _get_parsed_sequence(self, sequence_file: FilePath) -> SequenceReader:
        """
        Get the parsed sequences from the sequence file.
        
        Parameters:
        - sequence_file (FilePath): The path to the sequence file.

        Returns:
        - dict: A dictionary where the keys are the cluster indices and the values are the corresponding protein sequences.
        """
        self.check_has_cluster_intervals()

        if self.sequence_reader is None:
            self.sequence_reader = SequenceReader(sequence_file, seq_length=self.pae_matrix.shape[0])
        parsed_sequences = self.sequence_reader.clusters_to_fasta(self.cluster_intervals)

        return parsed_sequences


    def print_fasta(self, sequence_file: FilePath, header_name: Optional[str] = None, width: int = 60) -> None:
        """
        Print the sequences corresponding to each cluster in FASTA format.

        Parameters:
        - sequence_file (FilePath): The path to the sequence file.
        - header_name (str, optional): The name to add to the sequence headers. Defaults to self.sequence_reader.name.
        - width (int, optional): The width of the sequence lines. Defaults to 60.

        Raises:
        - ValueError: If the cluster intervals are not defined.
        """
        parsed_sequences = self._get_parsed_sequence(sequence_file)
        
        if header_name is None:
            header_name = self.sequence_reader.name

        for line in self._format_fasta_sequences(parsed_sequences, header_name, width):
            print(line)


    def save_fasta(self, sequence_file: FilePath, output_file: FilePath, header_name: Optional[str] = None, width: int = 60) -> None:
        """
        Save the sequences corresponding to each cluster in FASTA format to a file.

        Parameters:
        - sequence_file (FilePath): The path to the sequence file.
        - output_file (FilePath): The path to save the output file.
        - header_name (str, optional): The name to add to the sequence headers. Defaults to self.sequence_reader.name.
        - width (int, optional): The width of the sequence lines. Defaults to 60.

        Raises:
        - ValueError: If the cluster intervals are not defined.
        """
        parsed_sequences = self._get_parsed_sequence(sequence_file)
        
        if header_name is None:
            header_name = self.sequence_reader.name

        output_file = Path(output_file)

        with output_file.open('w') as f:
            for line in self._format_fasta_sequences(parsed_sequences, header_name, width):
                f.write(line + "\n")
