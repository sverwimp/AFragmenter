import pytest

import numpy as np
from afragmenter.graph import create_graph, cluster_graph, default_resolutions
from afragmenter import AFragmenter


def test_create_graph():
    """Test creating a simple graph with 3 nodes and 3 edges."""
    weights_matrix = np.array([[0, 1, 2], [1, 0, 3], [2, 3, 0]])

    graph = create_graph(weights_matrix)

    assert graph.is_weighted() == True
    assert graph.is_directed() == False
    assert graph.vcount() == 3
    assert graph.ecount() == 3
    assert graph.get_edgelist() == [(0, 1), (0, 2), (1, 2)]
    assert graph.es['weight'] == [1, 2, 3]


def test_cluster_graph_basic():
    """Test clustering a simple graph with 3 nodes and 2 clusters."""
    weights = np.array([
        [0, 20, 0],
        [20, 0, 0],
        [0, 0, 0]
    ])
    graph = create_graph(weights)
    partition = cluster_graph(graph, resolution=1.0, objective_function='modularity')

    memberships = partition.membership
    assert len(memberships) == 3
    assert memberships[0] == memberships[1]
    assert memberships[2] != memberships[0]


def test_cluster_graph_default_resolution():
    """Test clustering a simple graph with 2 nodes and 2 clusters using default resolution."""
    data = np.array([[0, 1], [1, 0]])
    a = AFragmenter(data)   

    _ = a.cluster(objective_function='modularity', min_size=0) # min_size must not be larger than the number of nodes, default is 10
    assert a.params.get('resolution') == default_resolutions['modularity']
    _ = a.cluster(objective_function='cpm', min_size=0)
    assert a.params.get('resolution') == default_resolutions['cpm']


def test_cluster_graph_invalid_objective_function():
    """Test that an invalid objective function raises a ValueError."""
    graph = create_graph(np.array([[0, 1], [1, 0]]))
    
    with pytest.raises(ValueError):
        cluster_graph(graph, objective_function='invalid')


def test_cluster_graph_no_edges():
    """Test clustering a graph with no edges."""
    weights = np.zeros((3,3))
    graph = create_graph(weights)
    
    partition = cluster_graph(graph, objective_function='modularity', n_iterations=1000)
    
    memberships = partition.membership
    assert len(set(memberships)) == 3


def test_cluster_graph_no_nodes():
    """Test clustering an empty graph."""
    weights = np.zeros((0,0))
    graph = create_graph(weights)
    
    partition = cluster_graph(graph, objective_function='modularity', n_iterations=1000)
    
    memberships = partition.membership
    assert len(set(memberships)) == 0


def test_cluster_graph_kwargs_valid():
    """Test passing additional keyword arguments to the clustering function."""
    graph = create_graph(np.array([[0, 1], [1, 0]]))
    
    _ = cluster_graph(graph, objective_function='modularity', resolution=0.5, n_iterations=10, beta=0.5)


def test_cluster_graph_kwargs_invalid():
    """Test passing an invalid keyword argument to the clustering function."""
    graph = create_graph(np.array([[0, 1], [1, 0]]))
    
    with pytest.raises(TypeError):
        _ = cluster_graph(graph, objective_function='modularity', seed=42) # seed is not a valid keyword argument
