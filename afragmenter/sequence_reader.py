import os
from io import StringIO
import re
from typing import Union, Optional, Iterable
from pathlib import Path

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


def _get_content(file_input: Union[FilePath, StringIO]) -> str:
    """
    Reads the content of a file or returns the input string.
    """
    if not isinstance(file_input, (str, Path, StringIO)):
        raise TypeError("file_input must be a string, Path, or StringIO object")

    if isinstance(file_input, Path) or (isinstance(file_input, str) and os.path.exists(file_input)):
        with open(file_input, 'r') as file:
            return file.read()
    elif isinstance(file_input, StringIO):
        return file_input.getvalue()
    elif isinstance(file_input, str):
        return file_input
    
    raise TypeError(f"Unsupported input type: {type(file_input)}")


def _determine_file_format(content: str) -> str:
    """
    Determines the format of a given file based on its first line.
    """
    first_line = _read_first_valid_line(content)
    if first_line is None:
        raise ValueError("Unable to infer file format, file is empty or contains only comments")

    pdb_firstlines = ('HEADER', 'ATOM', 'MODEL')

    if first_line.startswith('>'):
        return 'fasta'
    elif first_line.startswith(pdb_firstlines):
        return 'pdb'
    elif first_line.startswith('data_'):
        return 'mmcif'
    else:
        raise ValueError("Unsupported file format, please provide a FASTA, PDB or mmCIF file")


def _read_first_valid_line(content: str) -> Optional[str]:
    """
    Reads the first valid line from a file content.
    """
    for line in content.splitlines():
        stripped_line = line.strip()
        if stripped_line and not stripped_line.startswith('#'):
            return stripped_line
    return None


def _clean_name(name: str) -> str:
    """
    Cleans the name of a protein.
    """
    name = re.sub(r'^AF-', '', name)
    name = re.sub(r'-F\d+$', '', name)
    return '_'.join(name.split())


def _get_name_cif(cif_dict: dict) -> str:
    """
    Extracts the name of the protein from a CIF dictionary.
    """
    name = cif_dict.get('_entry.id', [''])[0]
    return _clean_name(name)


def _get_name_pdb(structure) -> str:
    """
    Extracts the name of the protein from a PDB structure.
    """
    name = structure.header.get('compound', {}).get('1', {}).get('molecule', '')
    if not name:
        name = structure.header.get('name', '')
    return _clean_name(name)


def _get_name_fasta(record: SeqIO.SeqRecord) -> str:
    """
    Extracts the name of the protein from a FASTA record.
    """
    return _clean_name(record.id)


def _translate_three2one_iter(residue: Iterable) -> str:
    """
    Translates a sequence of three-letter amino acid codes to one-letter codes.
    """
    return ''.join([THREE2ONE[res.upper()] for res in residue])


def _read_fasta_sequence(content: str, seq_length: Optional[int] = None) -> tuple[str, str]:
    """
    Reads the first protein sequence from a FASTA file or content.
    """
    handle = StringIO(content)
    
    if seq_length:
        for record in SeqIO.parse(handle, "fasta"):
            if len(record.seq) == seq_length:
                name = _get_name_fasta(record)
                return name, str(record.seq)
        raise ValueError(f"No sequence of length {seq_length} found in the FASTA content.")
    
    for record in SeqIO.parse(handle, "fasta"):
        name = _get_name_fasta(record)
        return name, str(record.seq)
    raise ValueError("No sequence found in FASTA content.")


def _read_pdb_sequence(content: str, query_chain: str = "A") -> tuple[str, str]:
    """
    Reads a protein sequence from a PDB file or content for a specific chain.
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('structure', StringIO(content))
    name = _get_name_pdb(structure)
    for model in structure:
        for chain in model:
            if chain.id == query_chain:
                sequence = ''.join([THREE2ONE[residue.resname] for residue in chain if residue.resname in THREE2ONE])
                return name, sequence
    raise ValueError(f"Chain {query_chain} not found in PDB file.")


def _read_mmcif_sequence(content: str) -> tuple[str, str]:
    """
    Reads a protein sequence from a mmCIF file or content.
    """
    mmcif_dict = MMCIF2Dict.MMCIF2Dict(StringIO(content))
    name = _get_name_cif(mmcif_dict)
    seq = mmcif_dict.get('_entity_poly.pdbx_seq_one_letter_code')
    if seq:
        return name, seq[0].replace('\n', '')
    
    seq = mmcif_dict.get('_entity_poly_seq.mon_id')
    if seq:
        return name, _translate_three2one_iter(seq)
    
    raise ValueError("Could not find the sequence in the mmCIF dictionary.")


class SequenceReader:
    def __init__(self, file_input: FilePath, format: str = 'auto', seq_length: Optional[int] = None, query_chain: str = "A"):
        self.name = ''
        self.query_chain = query_chain

        content = _get_content(file_input)

        if format.lower() == 'auto':
            format = _determine_file_format(content)
        self.format = format.lower()

        if self.format == 'fasta':
            self.name, self.sequence = _read_fasta_sequence(content, seq_length=seq_length)
        elif self.format == 'pdb':
            self.name, self.sequence = _read_pdb_sequence(content, query_chain=query_chain.upper())
        elif self.format in ['mmcif', 'cif']:
            self.name, self.sequence = _read_mmcif_sequence(content)
        else:
            raise ValueError("Unsupported file format, please provide a FASTA, PDB or mmCIF file")

        if seq_length and len(self.sequence) != seq_length:
            raise ValueError(f"Sequence length mismatch: expected {seq_length} residues, got {len(self.sequence)}")
        self.seq_length = len(self.sequence)

        if not self.name:
            if isinstance(file_input, (str, Path)) and os.path.isfile(file_input):
                self.name = Path(file_input).stem
            else:
                self.name = 'name'

    def clusters_to_fasta(self, cluster_intervals: dict) -> dict:
        """
        Parse the sequence file and return the sequences corresponding to each cluster.
        """
        parsed_sequences = {}
        for cluster_id, interval in cluster_intervals.items():
            cluster_sequence = []
            for i, (start, end) in enumerate(interval):
                if i > 0:
                    previous_end = interval[i-1][1]
                    if start > previous_end:
                        gap_size = start - previous_end - 1
                        cluster_sequence.extend(['-'] * gap_size)
                
                cluster_sequence.extend(self.sequence[start:end+1])
            parsed_sequences[cluster_id] = ''.join(cluster_sequence)
        return parsed_sequences
