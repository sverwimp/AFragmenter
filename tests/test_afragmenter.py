import pytest
from pathlib import Path
from io import StringIO
import numpy as np

from afragmenter.afragmenter import AFragmenter
from afragmenter.graph import default_resolutions


@pytest.fixture
def json_file():
    return Path(__file__).parent / "data" / "B1LFV6" / "B1LFV6.json"

@pytest.fixture
def cif_file():
    return Path(__file__).parent / "data" / "B1LFV6" / "B1LFV6.cif"

@pytest.fixture
def pdb_file():
    return Path(__file__).parent / "data" / "B1LFV6" / "B1LFV6.pdb"

@pytest.fixture
def pae_matrix():
    return np.array([
        [1, 2, 3],
        [2, 1, 2],
        [3, 2, 1]
    ])


def test_init_str_path(json_file):
    str_path = str(json_file)
    a = AFragmenter(str_path, threshold=1.0)
    assert hasattr(a, 'pae_matrix')
    assert hasattr(a, 'edge_weights_matrix')
    assert hasattr(a, 'graph')
    assert hasattr(a, 'params')
    assert a.params.get('threshold') == 1.0
    assert a.sequence_reader == None # Assigmed None in __init__, later reassigned


def test_init_pathlib_path(json_file):
    a = AFragmenter(json_file, threshold=1.0)
    assert a.params.get('threshold') == 1.0
    assert hasattr(a, 'pae_matrix')
    assert hasattr(a, 'edge_weights_matrix')
    assert hasattr(a, 'graph')
    assert hasattr(a, 'params')
    assert a.sequence_reader == None


def test_init_string_content(json_file):
    """Passing a string with the content of the file. Should fail"""
    with pytest.raises(FileNotFoundError):
        with open(json_file, 'r') as f:
            content = f.read()
        a = AFragmenter(content, threshold=1.0)


def test_init_dict_content():
    dict_content = {"predicted_aligned_error": [
        [0.1, 0.2, 0.3],
        [0.2, 0.1, 0.4],
        [0.3, 0.4, 0.1]
    ]}
    a = AFragmenter(dict_content, threshold=1.0)
    assert a.params.get('threshold') == 1.0
    assert hasattr(a, 'pae_matrix')
    assert hasattr(a, 'edge_weights_matrix')
    assert hasattr(a, 'graph')
    assert hasattr(a, 'params')
    assert a.sequence_reader == None
    

def test_init_stringio(json_file):
    """Passing a StringIO object with the content of the file"""
    with open(json_file, 'r') as f:
        content = f.read()
    a = AFragmenter(StringIO(content), threshold=1.0)
    assert hasattr(a, 'pae_matrix')
    assert hasattr(a, 'edge_weights_matrix')
    assert hasattr(a, 'graph')
    assert hasattr(a, 'params')
    assert a.params.get('threshold') == 1.0
    assert a.sequence_reader == None


def test_init_invalid_file():
    with pytest.raises(FileNotFoundError):
        _ = AFragmenter("invalid.json", threshold=1.0)


def test_init_invalid_args(json_file):
    with pytest.raises(TypeError):
        _ = AFragmenter(json_file, threshold="1.0", invalid="arg")


def test_threshold_range(pae_matrix):
    _ = AFragmenter._pae_transform(pae_matrix, threshold=0)
    _ = AFragmenter._pae_transform(pae_matrix, threshold=0.001)
    _ = AFragmenter._pae_transform(pae_matrix, threshold=5)
    _ = AFragmenter._pae_transform(pae_matrix, threshold=31.74)
    _ = AFragmenter._pae_transform(pae_matrix, threshold=31.75)

    with pytest.raises(ValueError, match="Threshold must be between 0 and 31.75"):
        _ = AFragmenter._pae_transform(pae_matrix, threshold=-0.001)

    with pytest.raises(ValueError, match="Threshold must be between 0 and 31.75"):
        _ = AFragmenter._pae_transform(pae_matrix, threshold=31.76)


def test_cluster_params_update(json_file):
    a = AFragmenter(json_file)
    assert a.params.get('resolution') == None
    a = a.cluster(resolution=0.8, n_iterations=2, objective_function="CPM", min_size=42, attempt_merge=False)
    assert a.params.get('resolution') == 0.8
    assert a.params.get('n_iterations') == 2
    assert a.params.get('objective_function') == "CPM"
    assert a.params.get('min_size') == 42
    assert a.params.get('attempt_merge') == False

    b = AFragmenter(json_file)
    assert b.params.get('resolution') == None
    b = b.cluster(objective_function='modularity')
    # Check if the default resolution is set
    assert b.params.get('resolution') == default_resolutions.get('modularity')


def test_cluster_intervals(json_file):
    """Test if the cluster method assigns the cluster_intervals attribute"""
    a = AFragmenter(json_file)
    assert not hasattr(a, 'cluster_intervals')
    a = a.cluster()
    assert hasattr(a, 'cluster_intervals')


def test_run(json_file):
    """Test alias for cluster method"""
    a = AFragmenter(json_file)
    assert not hasattr(a, 'cluster_intervals')
    a = a.run()
    assert hasattr(a, 'cluster_intervals')


def test_no_result_output(json_file):
    a = AFragmenter(json_file)
    
    with pytest.raises(ValueError, match="No clustering results found, please run the cluster method first"):
        a.print_result()

    with pytest.raises(ValueError, match="No clustering results found, please run the cluster method first"):
        a.plot_result()
    
    with pytest.raises(ValueError, match="No clustering results found, please run the cluster method first"):
        a.save_result("results.csv")

    with pytest.raises(ValueError, match="No clustering results found, please run the cluster method first"):
        a.save_fasta("sequence_file.cif", "output.fasta")

    with pytest.raises(ValueError, match="No clustering results found, please run the cluster method first"):
        a.print_fasta("sequence_file.cif")
    




