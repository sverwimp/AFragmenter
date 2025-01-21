from pathlib import Path
from typing import Union
from collections.abc import Iterable

from Bio.PDB import MMCIF2Dict, PDBParser
from Bio import SeqIO

FilePath = Union[str, Path]

# This dictionary was obtained from the hh-suite codebase: https://github.com/soedinglab/hh-suite/blob/master/scripts/cif2fasta.py
THREE2ONE = {
    'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K', 'ILE': 'I', 'PRO': 'P', 
    'THR': 'T', 'PHE': 'F', 'ASN': 'N', 'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 
    'TRP': 'W', 'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M', 'MSE': 'M',
    'HYP': 'P', 'MLY': 'K', 'SEP': 'S', 'TPO': 'T', 'CSO': 'C', 'PTR': 'Y', 'KCX': 'K',
    'CME': 'C', 'CSD': 'A', 'CAS': 'C', 'MLE': 'L', 'DAL': 'A', 'CGU': 'E', 'DLE': 'L',
    'FME': 'M', 'DVA': 'V', 'OCS': 'C', 'DPR': 'P', 'MVA': 'V', 'TYS': 'Y', 'M3L': 'K',
    'SMC': 'C', 'ALY': 'K', 'CSX': 'C', 'DCY': 'C', 'NLE': 'L', 'DGL': 'E', 'DSN': 'S',
    'CSS': 'C', 'DLY': 'K', 'MLZ': 'K', 'DPN': 'F', 'DAR': 'R', 'PHI': 'F', 'IAS': 'D',
    'DAS': 'D', 'HIC': 'H', 'MP8': 'P', 'DTH': 'T', 'DIL': 'I', 'MEN': 'N', 'DTY': 'Y',
    'CXM': 'M', 'DGN': 'G', 'DTR': 'W', 'SAC': 'S', 'DSG': 'N', 'MME': 'M', 'MAA': 'A',
    'YOF': 'Y', 'FP9': 'P', 'FVA': 'V', 'MLU': 'L', 'OMY': 'Y', 'FGA': 'E', 'MEA': 'F',
    'CMH': 'C', 'DHI': 'H', 'SEC': 'C', 'OMZ': 'Y', 'SCY': 'C', 'MHO': 'M', 'MED': 'M',
    'CAF': 'C', 'NIY': 'Y', 'OAS': 'S', 'SCH': 'C', 'MK8': 'L', 'SME': 'M', 'LYZ': 'K'
}


def translate_three2one_iter(residue: Iterable) -> str:
    """
    Translates a sequence of three-letter amino acid codes to one-letter codes.
    
    Parameters:
    - residue (Iterable): A sequence of three-letter amino acid codes.

    Returns:
    - str: A string containing the translated one-letter amino acid codes.
    """
    return ''.join([THREE2ONE[res] for res in residue])


def _read_mmcif_sequence(mmcif_file: FilePath):
    """
    Reads in a protin sequence from a mmCIF file.
    The function first tries to read the sequence from the '_entity_poly.pdbx_seq_one_letter_code' field.
    If the field is not found, it tries to read the sequence from the '_entity_poly_seq.mon_id' field.

    Parameters:
    - mmcif_file (FilePath): The path to the mmCIF file.

    Returns:
    - str: The protein sequence.

    Raises:
    - ValueError: If the sequence is not found in the mmCIF file.
    """
    mmcif_dict = MMCIF2Dict.MMCIF2Dict(str(mmcif_file))
    
    seq = mmcif_dict.get('_entity_poly.pdbx_seq_one_letter_code')
    if seq:
        return seq[0].replace('\n', '')
    
    seq = mmcif_dict.get('_entity_poly_seq.mon_id')
    if seq:
        return translate_three2one_iter(seq)
    
    raise ValueError("Could not find the sequence in the mmCIF file.")


def read_mmcif_sequence(mmcif_file: FilePath, seq_length: int) -> str:
    """
    Intermediary function that calls _read_mmcif_sequence to read in the sequence and checks the sequence length.

    Parameters:
    - mmcif_file (FilePath): The path to the PDB file.
    - seq_length (int): The expected length of the protein sequence.

    Returns:
    - str: The protein sequence.
    
    Raises:
    - ValueError: If the sequence length does not match the expected length.
    """
    seq = _read_mmcif_sequence(mmcif_file)
    if len(seq) != seq_length:
        raise ValueError(f"Sequence length mismatch: expected {seq_length} residues, got {len(seq)}")
    return seq


def _read_pdb_sequence(pdb_file: FilePath, query_chain: str) -> str:
    """
    Reads a protein sequence from a PDB file.
    
    Parameters:
    - pdb_file (FilePath): The path to the PDB file.
    - query_chain (str): The chain identifier to use when reading the PDB file.

    Returns:
    - str: The protein sequence.

    Raises:
    - ValueError: If the chain identifier is not found in the PDB file.
    """
    parser = PDBParser()
    structure = parser.get_structure('structure', pdb_file)
    for model in structure:
        for chain in model:
            if chain.id == query_chain:
                return ''.join([THREE2ONE[residue.resname] for residue in chain.get_residues()])
    else:
        raise ValueError(f"Chain {query_chain} not found in the PDB file.")


def read_pdb_sequence(pdb_file: FilePath, query_chain: str, seq_length: int) -> str:
    """
    Intermediary function that calls _read_pdb_sequence to read in the sequence and checks the sequence length.

    Parameters:
    - pdb_file (FilePath): The path to the PDB file.
    - query_chain (str): The chain identifier to use when reading the PDB file.
    - seq_length (int): The expected length of the protein sequence.

    Returns:
    - str: The protein sequence.

    Raises:
    - ValueError: If the sequence length does not match the expected length.
    """
    seq = _read_pdb_sequence(pdb_file, query_chain)
    if len(seq) != seq_length:
        raise ValueError(f"Sequence length mismatch: expected {seq_length} residues, got {len(seq)}")
    return seq


def read_fasta_sequence(fasta_file: FilePath, seq_length: int) -> str:
    """
    Reads a protein sequence from a FASTA file.
    Iterates over the records in the FASTA file and returns the first sequence that matches the expected length.

    Parameters:
    - fasta_file (FilePath): The path to the FASTA file.
    - seq_length (int): The expected length of the protein sequence.

    Returns:
    - str: The protein sequence.

    Raises:
    - ValueError: If no sequence of the expected length is found in the FASTA file.
    """
    for record in SeqIO.parse(fasta_file, 'fasta'):
        if len(record.seq) == seq_length:
            return str(record.seq)
    raise ValueError(f"No sequence of length {seq_length} found in the FASTA file.")
    

def read_first_valid_line(file_path: FilePath) -> Union[str, None]:
    """
    Reads the first valid line from a file. 
    A valid line is a line that is not empty and does not start with a comment character ('#').

    Parameters:
    - file_path (FilePath): The path to the file.
    
    Returns:
    - Union[str, None]: The first valid line from the file, or None if no valid line is found.
    """
    with open(file_path, 'r') as file:
        for line in file:
            stripped_line = line.strip()
            if stripped_line and not stripped_line.startswith('#'):
                return stripped_line
    return None


def determine_file_format(file_path: FilePath) -> str:
    """
    Determines the format of a given file based on its first line.
    Supported formats are FASTA, PDB, and mmCIF.

    Parameters:
    - file_path (FilePath): The path to the file.
    
    Returns:
    - str: The format of the file. Possible values are 'fasta', 'PDB', or 'mmCIF'
    
    Raises:
    - ValueError: If the file is empty or contains only comments and the format cannot be inferred.
    - ValueError: If the file format is not supported.
    """
    first_line = read_first_valid_line(file_path)
    if first_line is None:
        raise ValueError("Unable to infer file format, file is empty or contains only comments")

    if first_line.startswith('>'):
        return 'FASTA'
    elif first_line.startswith('HEADER') or first_line.startswith('ATOM'):
        return 'PDB'
    elif first_line.startswith('data_'):
        return 'mmCIF'
    else:
        raise ValueError("Unsupported file format, please provide a FASTA, PDB or mmCIF file")


def read_sequence(file_path: FilePath, seq_length: int, query_chain: str = "A") -> str:
    """
    Reads a protein sequence from a file.
    
    Parameters:
    - file_path (FilePath): The path to the file.
    - seq_length (int): The expected length of the protein sequence.
    - query_chain (str): The chain identifier to use when reading a PDB file. Default is 'A'.
    
    Returns:
    - str: The protein sequence.
    
    Raises:
    - ValueError: If the file format is not supported.
    """
    file_format = determine_file_format(file_path)
    
    if file_format == 'FASTA':
        return read_fasta_sequence(file_path, seq_length)
    elif file_format == 'PDB':
        return read_pdb_sequence(file_path, query_chain=query_chain.upper(), seq_length=seq_length)
    elif file_format == 'mmCIF':
        return read_mmcif_sequence(file_path, seq_length)
    else:
        raise ValueError("Unsupported file format, please provide a FASTA, PDB or mmCIF file")
