import pytest
from unittest.mock import patch, Mock

from afragmenter.afdb_client import fetch_url_content, fetch_afdb_data


def test_url_content():
    url = "https://example.com"
    mock_response = Mock()
    mock_response.ok = True
    mock_response.text = "content"
    
    with patch('requests.get', return_value=mock_response):
        response = fetch_url_content(url)
        assert response.text == "content"


def test_fetch_afdb_data_success():
    uniprot_id = "P12345"
    mock_afdb_response = Mock()
    mock_afdb_response.ok = True
    mock_afdb_response.json.return_value = [{
        'paeDocUrl': 'https://example.com/pae',
        'cifUrl': 'https://example.com/cif',
        'pdbUrl': 'https://example.com/pdb'
    }]
    
    mock_pae_response = Mock()
    mock_pae_response.ok = True
    mock_pae_response.json.return_value = {"pae": "data"}
    
    mock_structure_response = Mock()
    mock_structure_response.ok = True
    mock_structure_response.text = "structure data"
    
    with patch('requests.get', side_effect=[mock_afdb_response, mock_pae_response, mock_structure_response]):
        data = fetch_afdb_data(uniprot_id)
        assert data.pae_data == {"pae": "data"}
        assert data.structure_data == "structure data"


def test_fetch_afdb_data_no_pae():
    uniprot_id = "P12345"
    mock_afdb_response = Mock()
    mock_afdb_response.ok = True
    mock_afdb_response.json.return_value = [{
        'paeDocUrl': '',
        'cifUrl': 'https://example.com/cif',
        'pdbUrl': 'https://example.com/pdb'
    }]
    
    with patch('requests.get', return_value=mock_afdb_response):
        with pytest.raises(ValueError, match="No PAE data available for P12345"):
            fetch_afdb_data(uniprot_id)


def test_fetch_afdb_data_no_structure():
    uniprot_id = "P12345"
    mock_afdb_response = Mock()
    mock_afdb_response.ok = True
    mock_afdb_response.json.return_value = [{
        'paeDocUrl': 'https://example.com/pae',
        'cifUrl': '',
        'pdbUrl': ''
    }]
    
    with patch('requests.get', return_value=mock_afdb_response):
        with pytest.raises(ValueError, match="No CIF or PDB data available for P12345"):
            fetch_afdb_data(uniprot_id)