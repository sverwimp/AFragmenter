import json
import os
from pathlib import Path
from typing import Union
from io import StringIO

import numpy as np


FilePath = Union[str, Path]


def _validate_pae(pae: np.ndarray) -> None:
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
    if pae.size == 0:
        raise ValueError("PAE matrix is empty")
    if pae.ndim != 2:
        raise ValueError("PAE matrix must be 2D")
    if pae.shape[0] != pae.shape[1]:
        raise ValueError("PAE matrix must be square")
    if np.min(pae) < 0:
        raise ValueError("PAE values must be non-negative")


def process_pae_data(pae_data: Union[list, dict]) -> np.ndarray:
    """
    Parameters:
    - pae_data (list, dict): json-like format containing the PAE data

    Returns:
    - np.ndarray: The PAE matrix
    
    Raises:
    - ValueError: If the PAE data is not found in the JSON file.
    """
    # AF2 format loads as a list containing a dictionary, AF3 and colabfold directly load the dictionary
    if isinstance(pae_data, list):
        pae_data = pae_data[0]

    if not isinstance(pae_data, dict):
        raise TypeError("Invalid PAE data format, expected a dictionary")

    # AFDB v1 and v2 have different keys for the PAE data
    if "distance" in pae_data:
        nrows = max(pae_data.get('residue1'))
        pae_matrix = np.zeros((nrows + 1, nrows + 1))
        for r, c, v in zip(pae_data.get('residue1'), pae_data.get('residue2'), pae_data.get('distance')):
            pae_matrix[r, c] = v
    else:
        pae = pae_data.get("predicted_aligned_error") or pae_data.get("pae")
        if pae is None:
            raise ValueError("PAE data not found in JSON file")
        pae_matrix = np.stack(pae, axis=0)

    _validate_pae(pae_matrix)
    return pae_matrix


def load_pae(json_source: Union[FilePath, StringIO]) -> np.ndarray:
    """
    Reads a json file and calls the _process_pae_data function to
    load the Predicted Aligned Error (PAE) data as a numpy array

    Parameters:
    - json_source (FilePath, StringIO): The path to the JSON file. (str or Path)

    Returns:
    - np.ndarray: The PAE matrix.

    Raises:
    - FileNotFoundError: If the JSON file does not exist.
    """

    if isinstance(json_source, StringIO):
        pae_data = json.load(json_source)

    elif isinstance(json_source, (str, Path)):
        if not os.path.exists(json_source):
            raise FileNotFoundError(f"{json_source} not found")
        with open(json_source, "r") as f:
            pae_data = json.load(f)

    else:
        raise TypeError("Invalid json_source type, expected str, Path or StringIO")

    return process_pae_data(pae_data)