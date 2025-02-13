from collections import namedtuple

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
    

def fetch_afdb_data(uniprot_id: str) -> namedtuple:
    """
    Fetch the PAE and structure data for a given UniProt ID from the AlphaFold Database.
    
    Parameters
    - uniprot_id (str): The UniProt ID to fetch the data for.
    
    Returns
    - data (namedtuple): A named tuple containing the PAE (pae_data) and structure data (structure_data).
    
    Raises
    - ValueError: If no PAE data is available for the given UniProt ID.
    - ValueError: If no CIF or PDB data is available for the given UniProt ID.
    """
    afdb_prediction_url = f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id.upper()}"
    afdb_prediction_json = fetch_url_content(afdb_prediction_url).json()
    afdb_prediction = afdb_prediction_json[0] if isinstance(afdb_prediction_json, list) else afdb_prediction_json
    
    pae_url = afdb_prediction.get('paeDocUrl', '')
    if not pae_url:
        raise ValueError(f"No PAE data available for {uniprot_id}")
    
    cif_url = afdb_prediction.get('cifUrl', '')
    pdb_url = afdb_prediction.get('pdbUrl', '')
    if not cif_url or not pdb_url:
        raise ValueError(f"No CIF or PDB data available for {uniprot_id}")
    
    pae_data_json = fetch_url_content(pae_url).json()
    pae_data = pae_data_json[0] if isinstance(pae_data_json, list) else pae_data_json
    structure_data = fetch_url_content(cif_url).text if cif_url else fetch_url_content(pdb_url).text
    
    Data = namedtuple('AFDB_Data', ['pae_data', 'structure_data'])
    return Data(pae_data, structure_data)
    
    
    
    
    
    
    
    
    
    
    