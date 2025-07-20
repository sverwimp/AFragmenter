from pathlib import Path
from typing import Optional, Union, Tuple
from io import StringIO

import numpy as np
from matplotlib import image, axes

from .pae_handler import process_pae_data, load_pae
from .sequence_reader import SequenceReader
from .graph import create_graph, cluster_graph
from .intervals import find_cluster_intervals, filter_cluster_intervals
from afragmenter import plotting
from .result import ClusteringResult


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
    - sequence_reader (SequenceReader): The SequenceReader object used to read the sequence file. This is automatically created
                                        when the print_fasta or save_fasta methods are called.

    Methods:
    - cluster: Cluster the graph using the Leiden algorithm.
    - run: Alias for the cluster method.
    - plot_pae: Plot the Predicted Aligned Error matrix as a heatmap.
    """

    def __init__(self, pae_matrix: Union[np.ndarray, FilePath, list, dict, StringIO], threshold: float = 5.0, sequence_file: Optional[FilePath] = None):
        if isinstance(pae_matrix, (list, dict)):
            pae_matrix = process_pae_data(pae_matrix)
        elif isinstance(pae_matrix, (str, Path, StringIO)):
            pae_matrix = load_pae(pae_matrix)
        elif not isinstance(pae_matrix, np.ndarray):
            raise TypeError("pae_matrix must be a numpy array, a file path, or a dictionary containing PAE data")
        
        self.pae_matrix = pae_matrix
        self.edge_weights_matrix = AFragmenter._pae_transform(pae_matrix, threshold)
        self.graph = create_graph(self.edge_weights_matrix)
        self.params = {"threshold": threshold}
        self.sequence_reader = None
        if sequence_file:
            self.sequence_reader = SequenceReader(sequence_file, seq_length=self.pae_matrix.shape[0])
        

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
                **kwargs) -> 'ClusteringResult':
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
        params = self.params.copy()
        params.update({ 
            "resolution": resolution,
            "objective_function": objective_function,
            "n_iterations": n_iterations,
            "min_size": min_size,
            "attempt_merge": attempt_merge
        })
        
        clusters, cluster_params = cluster_graph(graph=self.graph, 
                                 resolution=resolution, 
                                 n_iterations=n_iterations,
                                 objective_function=objective_function,
                                 return_params=True,
                                 **kwargs)
        params.update(cluster_params) # Update the parameters with the actual values used for clustering
        cluster_intervals = find_cluster_intervals(clusters)
        cluster_intervals = filter_cluster_intervals(intervals=cluster_intervals, 
                                                          min_size=min_size, 
                                                          attempt_merge=attempt_merge)
        return ClusteringResult(self.pae_matrix, cluster_intervals, params, self.sequence_reader)
    

    def run(self, resolution: Optional[float] = None, 
            objective_function: str = "modularity", 
            n_iterations: int = -1, 
            min_size: int = 10,
            attempt_merge: bool = True,
            **kwargs) -> 'ClusteringResult':
        """Alias for the cluster method."""
        return self.cluster(resolution=resolution, 
                            objective_function=objective_function, 
                            n_iterations=n_iterations, 
                            min_size=min_size, 
                            attempt_merge=attempt_merge,
                            **kwargs)

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


    def __repr__(self) -> str:
        return (f"AFragmenter(pae_matrix=np.ndarray(shape={self.pae_matrix.shape}), "
                f"params={self.params})")


    def __str__(self) -> str:
        return (f"AFragmenter(pae_matrix=np.ndarray(shape={self.pae_matrix.shape}), "
                f"params={self.params})")
