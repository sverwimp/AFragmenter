from collections import namedtuple
from typing import Literal

import requests


def fetch_url_content(url: str) -> requests.Response:
    """
    Fetch the content of a URL and return the response object.
    
    Parameters
    - url (str): The URL to fetch the content from.
    
    Returns
    - response (requests.Response): The response object.
    
    Raises
    - requests.HTTPError: If the response status code is not ok.
    """
    response = requests.get(url)
    
    if not response.ok:
        response.raise_for_status()
    
    return response

def fetch_afdb_data(
    uniprot_id: str, 
    structure_format: Literal['cif', 'mmcif', 'pdb'] = 'mmcif'
) -> namedtuple:
    """
    Fetch the PAE and structure data for a given UniProt ID from the AlphaFold Database.
    
    Parameters
    - uniprot_id (str): The UniProt ID to fetch the data for.
    - structure_format (Literal['cif', 'mmcif', 'pdb']): The desired format for the structure data.
                                                 Defaults to 'mmcif'.
    
    Returns
    - data (namedtuple): A named tuple containing the PAE (pae_data) and structure data (structure_data).
    
    Raises
    - ValueError: If no PAE data is available for the given UniProt ID.
    - ValueError: If no CIF or PDB data is available for the given UniProt ID based on the requested format.
    - requests.HTTPError: If there's an issue fetching data from the URLs.
    """
    afdb_prediction_url = f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id.upper()}"
    afdb_prediction_json = fetch_url_content(afdb_prediction_url).json()
    
    # Handle cases where the API might return a single object directly or a list with one object
    afdb_prediction = afdb_prediction_json[0] if isinstance(afdb_prediction_json, list) and afdb_prediction_json else afdb_prediction_json
    
    if not afdb_prediction:
        raise ValueError(f"No prediction data found for UniProt ID: {uniprot_id}")

    pae_url = afdb_prediction.get('paeDocUrl', '')
    if not pae_url:
        raise ValueError(f"No PAE data URL available for {uniprot_id}")
    
    structure_download_url = None

    if structure_format == 'mmcif' or structure_format == 'cif':
        cif_url = afdb_prediction.get('cifUrl', '')
        if cif_url:
            structure_download_url = cif_url
        else:
            raise ValueError(f"No CIF data available for {uniprot_id} as requested.")
    elif structure_format == 'pdb':
        pdb_url = afdb_prediction.get('pdbUrl', '')
        if pdb_url:
            structure_download_url = pdb_url
        else:
            raise ValueError(f"No PDB data available for {uniprot_id} as requested.")
    else:
        # This case should ideally not be reached due to Literal type hint,
        # but good for robustness.
        raise ValueError(f"Invalid structure_format: {structure_format}. Must be 'mmcif' or 'pdb'.")

    pae_data_json = fetch_url_content(pae_url).json()
    # Handle cases where PAE data might also be a list or single object
    pae_data = pae_data_json[0] if isinstance(pae_data_json, list) and pae_data_json else pae_data_json

    structure_data = fetch_url_content(structure_download_url).text
    
    Data = namedtuple('AFDB_Data', ['pae_data', 'structure_data'])
    return Data(pae_data, structure_data)

    
    
    
    
    
    
    
    
    
    