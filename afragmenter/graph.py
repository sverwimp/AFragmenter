import igraph
from typing import Union, Optional
import numpy as np


default_resolutions = {"modularity": 0.7, "cpm": 0.3}


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

    return graph


def cluster_graph(graph: igraph.Graph,
            resolution: Optional[float] = None, 
            n_iterations: int = -1, 
            objective_function: str = "modularity",
            return_params: bool = False,
            **kwargs) -> Union[igraph.VertexClustering, tuple[igraph.VertexClustering, dict]]:
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
    - kwargs: Additional keyword arguments for the community_leiden function from igraph.

    Returns:
    - igraph.VertexClustering: The resulting vertex clustering.
    - dict: A dictionary containing the parameters used for clustering. Only returned if return_params is True.
    """

    objective_funtions = list(default_resolutions.keys())
    if objective_function.lower() not in objective_funtions:
        raise ValueError(f"Objective function must be one of {objective_funtions}, got {objective_function}")

    if resolution is None:
        resolution = default_resolutions[objective_function.lower()]

    partition = graph.community_leiden(
        weights="weight",
        resolution=resolution, 
        n_iterations=n_iterations, 
        objective_function=objective_function,
        **kwargs
    )

    if return_params:
        return partition, {"resolution": resolution, "n_iterations": n_iterations, "objective_function": objective_function}
    return partition
