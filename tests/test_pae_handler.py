import json

import numpy as np
import pytest

from afragmenter.pae_handler import load_pae, process_pae_data, _validate_pae


def test_validate_pae_valid_input():
    """Test that a valid PAE matrix does not raise any exceptions."""
    pae = np.array([[1.0, 2.0, 3.0],
                   [4.0, 5.0, 6.0],
                   [7.0, 8.0, 9.0]])
    
    _validate_pae(pae)


def test_validate_pae_invalid_type():
    """Test that a non-numpy array input raises TypeError."""
    # Test with a list instead of numpy array
    pae = [[1.0, 2.0], [3.0, 4.0]]
    
    with pytest.raises(TypeError, match="pae must be a numpy array"):
        _validate_pae(pae)


def test_validate_pae_invalid_ndim():
    """Test that non-2D matrices raise ValueError."""
    # Test with 1D array
    pae = np.array([1.0, 2.0, 3.0])
    
    with pytest.raises(ValueError, match="PAE matrix must be 2D"):
        _validate_pae(pae)


def test_validate_pae_invalid_shape():
    """Test that non-square matrices raise ValueError."""
    # Test with a rectangular (non-square) matrix
    pae = np.array([[1.0, 2.0],
                   [3.0, 4.0],
                   [5.0, 6.0]])
    
    with pytest.raises(ValueError, match="PAE matrix must be square"):
        _validate_pae(pae)


def test_validate_pae_negative_values():
    """Test that matrices with negative values raise ValueError."""
    # Test with a valid shape but contains negative values
    pae = np.array([[1.0, -2.0],
                   [-3.0, 4.0]])
    
    with pytest.raises(ValueError, match="PAE values must be non-negative"):
        _validate_pae(pae)


def test_validate_pae_empty_array():
    """Test that an empty array raises appropriate exceptions."""
    # Test with an empty 2D array
    pae = np.array([])
    
    with pytest.raises(ValueError, match="PAE matrix is empty"):
        _validate_pae(pae)


def test_validate_pae_zero_matrix():
    """Test that a zero matrix (all zeros) does not raise exceptions."""
    pae = np.zeros((2, 2))
    _validate_pae(pae)


def test_process_pae_data_afdb_v1_v2_format():
    """Test PAE processing of AFDB v1 and v2 formats."""
    pae = [{
        'residue1': [0, 1, 2],
        'residue2': [0, 1, 2],
        'distance': [0.5, 0.6, 0.7]
    }]
    expected_matrix = np.array([
        [0.5, 0.0, 0.0],
        [0.0, 0.6, 0.0],
        [0.0, 0.0, 0.7]
    ])
    result_matrix = process_pae_data(pae)
    np.testing.assert_array_equal(result_matrix, expected_matrix)


def test_process_pae_data_afdb_format():
    """Test PAE processing of AFDB v3 and v4 format."""
    pae = {
        'predicted_aligned_error': [
            [0.5, 0.6, 0.7],
            [0.6, 0.5, 0.8],
            [0.7, 0.8, 0.5]
        ]
    }
    expected_matrix = np.array([
        [0.5, 0.6, 0.7],
        [0.6, 0.5, 0.8],
        [0.7, 0.8, 0.5]
    ])
    result_matrix = process_pae_data(pae)
    np.testing.assert_array_almost_equal(result_matrix, expected_matrix)


def test_process_pae_data_list():
    """Test PAE processing of a list of PAE matrices."""
    pae = [{'predicted_aligned_error': 
            [[0.5, 0.6, 0.7], 
             [0.6, 0.5, 0.8], 
             [0.7, 0.8, 0.5]]}]
    expected_matrix = np.array([
        [0.5, 0.6, 0.7],
        [0.6, 0.5, 0.8],
        [0.7, 0.8, 0.5]
    ])
    result_matrix = process_pae_data(pae)
    np.testing.assert_array_almost_equal(result_matrix, expected_matrix)


def test_process_pae_data_missing_pae():
    """Test that missing PAE data raises ValueError."""
    pae_data = {}
    with pytest.raises(ValueError, match="PAE data not found in JSON file"):
        process_pae_data(pae_data)


def test_process_pae_data_invalid_format():
    """Test that invalid PAE data format raises TypeError."""
    pae_data = "invalid_format"
    with pytest.raises(TypeError, match="Invalid PAE data format, expected a dictionary"):
        process_pae_data(pae_data)


def test_load_pae_valid_input(tmp_path):
    """Test that a valid PAE JSON file."""
    pae = {
        'predicted_aligned_error': [
            [0.5, 0.6, 0.7],
            [0.6, 0.5, 0.8],
            [0.7, 0.8, 0.5]
        ]
    }
    expected_matrix = np.array([
        [0.5, 0.6, 0.7],
        [0.6, 0.5, 0.8],
        [0.7, 0.8, 0.5]
    ])
    json_file = tmp_path / "pae.json"
    with open(json_file, 'w') as f:
        json.dump(pae, f)
    
    result_matrix = load_pae(json_file)
    np.testing.assert_array_almost_equal(result_matrix, expected_matrix)


def test_load_pae_file_not_found():
    """Test loading a non-existent PAE JSON file."""
    with pytest.raises(FileNotFoundError, match="not found"):
        load_pae("non_existent_file.json")

def test_load_pae_invalid_json(tmp_path):
    """Test loading an invalid JSON file."""
    invalid_json_file = tmp_path / "invalid.json"
    with open(invalid_json_file, 'w') as f:
        f.write("invalid json content")
    
    with pytest.raises(json.JSONDecodeError):
        load_pae(invalid_json_file)

def test_load_pae_missing_pae_data(tmp_path):
    """Test loading a JSON file with missing PAE data."""
    invalid_pae_data = {
        'some_other_key': [
            [0.5, 0.6, 0.7],
            [0.6, 0.5, 0.8],
            [0.7, 0.8, 0.5]
        ]
    }
    
    # Create a temporary JSON file
    json_file = tmp_path / "invalid_pae.json"
    with open(json_file, 'w') as f:
        json.dump(invalid_pae_data, f)
    
    with pytest.raises(ValueError, match="PAE data not found in JSON file"):
        load_pae(json_file)


def test_load_pae_invalid_format(tmp_path):
    """Test loading a JSON file with invalid format."""
    invalid_pae_data = "invalid_format"
    
    # Create a temporary JSON file
    json_file = tmp_path / "invalid_pae.json"
    with open(json_file, 'w') as f:
        f.write(invalid_pae_data)
    
    with pytest.raises(json.JSONDecodeError):
        load_pae(json_file)