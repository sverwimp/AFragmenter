import os
from io import StringIO
import re
from typing import Union, Optional, Iterable
from pathlib import Path

from Bio.PDB import MMCIF2Dict, PDBParser
from Bio import SeqIO


FilePath = Union[str, Path]


class SequenceReader:
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


    def __init__(self, file_input: FilePath, format: str = 'auto', seq_length: Optional[int] = None, query_chain: str = "A"):
        self.name = ''
        self.query_chain = query_chain

        content = SequenceReader.get_content(file_input)

        if format.lower() == 'auto':
            format = SequenceReader.determine_file_format(content)
        self.format = format

        self.seq_length = seq_length # Used to check if the sequence length matches the expected length when reading a FASTA file
        # The sequence is read during initialization
        # And the name of the protein is set during the reading process
        self.sequence = self.read_sequence(content, self.format, query_chain=query_chain)

        if seq_length and len(self.sequence) != seq_length:
            raise ValueError(f"Sequence length mismatch: expected {seq_length} residues, got {len(self.sequence)}")
        self.seq_length = len(self.sequence)

        if self.name == '':
            if os.path.isfile(file_input):
                self.name = Path(file_input).stem
            else:
                self.name = 'name'
        

    @staticmethod
    def get_content(file_input: Union[FilePath, StringIO]) -> str:
        """
        Reads the content of a file or returns the input string.

        Parameters:
        - file_input (FilePath or StringIO): The path to the file or the content of the file.

        Returns:
        - str: The content of the file.
        """
        if not isinstance(file_input, (str, Path, StringIO)):
            raise TypeError("file_input must be a string or a Path object")

        if os.path.exists(file_input):
            with open(file_input, 'r') as file:
                return file.read()
        return file_input


    @staticmethod
    def determine_file_format(content: str) -> str:
        """
        Determines the format of a given file based on its first line.
        Supported formats are FASTA, PDB, and mmCIF.

        Parameters:
        - content (str): The content of a file.
        
        Returns:
        - str: The format of the file. Possible values are 'fasta', 'pdb', or 'mmcif'
        
        Raises:
        - ValueError: If the file is empty or contains only comments and the format cannot be inferred.
        - ValueError: If the file format is not supported.
        """

        first_line = SequenceReader.read_first_valid_line(content)
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


    @staticmethod
    def read_first_valid_line(content: str) -> Optional[str]:
        """
        Reads the first valid line from a file.
        A valid line is a line that is not empty and does not start with a comment character ('#').
        
        Parameters:
        - content (str): The content of a file.
        
        Returns:
        - Union[str, None]: The first valid line from the file, or None if no valid line is found.
        """
        for line in content.splitlines():
            stripped_line = line.strip()
            if stripped_line and not stripped_line.startswith('#'):
                return stripped_line
        return None


    def read_sequence(self, content: str, format: str, query_chain: str = "A"):
        """
        Reads a protein sequence from a file.
        
        Parameters:
        - content (str): The content of the file.
        - format (str): The format of the file. Possible values are 'fasta', 'pdb', or 'mmcif'.
        - query_chain (str): The chain identifier to use when reading a PDB file. Default is 'A'.
        
        Returns:
        - str: The protein sequence.
        
        Raises:
        - ValueError: If the file format is not supported.
        """
        format = format.lower()
        if format == 'fasta':
            return self.read_fasta_sequence(content, seq_length=self.seq_length)
        elif format == 'pdb':
            return self.read_pdb_sequence(content, query_chain=query_chain.upper())
        elif format == 'mmcif' or format == 'cif':
            return self.read_mmcif_sequence(content)
        else:
            raise ValueError("Unsupported file format, please provide a FASTA, PDB or mmCIF file")


    def read_fasta_sequence(self, content: str, seq_length: Optional[int] = None) -> str:
        """
        Reads the first protein sequence from a FASTA file or content.
        If seq_length is provided, only sequences with a length equal to seq_length are considered.

        Parameters:
        - content (str): The content of the FASTA file.
        - seq_length (int) [Optional]: The expected length of the protein sequence.

        Returns:
        - str: The protein sequence.

        Raises:
        - ValueError: If no sequence is found in the FASTA content.
        """
        handle = StringIO(content)
        
        # If seq_length is provided, only consider sequences with the specified length
        if seq_length:
            for record in SeqIO.parse(handle, "fasta"):
                if len(record.seq) == seq_length:
                    self.name = SequenceReader.get_name_fasta(record)
                    return str(record.seq)
            raise ValueError(f"No sequence of length {seq_length} found in the FASTA content.")
        
        # Else, return the first sequence found
        for record in SeqIO.parse(handle, "fasta"):
            self.name = SequenceReader.get_name_fasta(record)
            return str(record.seq)
        raise ValueError("No sequence found in FASTA content.")
    

    def read_pdb_sequence(self, content: str, query_chain: str = "A") -> str:
        """
        Reads a protein sequence from a PDB file or content for a specific chain.

        Parameters:
        - content (str): The content of the PDB file.
        - query_chain (str): The chain identifier to use when reading the PDB file.

        Returns:
        - str: The protein sequence.

        Raises:
        - ValueError: If the chain is not found in the PDB structure.
        """
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure('structure', StringIO(content))
        self.name = SequenceReader.get_name_pdb(structure)
        for model in structure:
            for chain in model:
                if chain.id == query_chain:
                    #self.name = SequenceReader.get_name_pdb_chain(chain)
                    return ''.join([SequenceReader.THREE2ONE[residue.resname] for residue in chain if residue.resname in SequenceReader.THREE2ONE])
        else:
            raise ValueError(f"Chain {query_chain} not found in PDB file.")


    def read_mmcif_sequence(self, content: str) -> str:
        """
        Reads a protein sequence from a mmCIF file or content.

        Parameters:
        - content (str): The content of the mmCIF file.

        Returns:
        - str: The protein sequence.

        Raises:
        - ValueError: If the sequence is not found in the mmCIF dictionary.
        """
        mmcif_dict = MMCIF2Dict.MMCIF2Dict(StringIO(content))
        seq = mmcif_dict.get('_entity_poly.pdbx_seq_one_letter_code')
        if seq:
            self.name = SequenceReader.get_name_cif(mmcif_dict)
            return seq[0].replace('\n', '')
        
        seq = mmcif_dict.get('_entity_poly_seq.mon_id')
        if seq:
            self.name = SequenceReader.get_name_cif(mmcif_dict)
            return SequenceReader.translate_three2one_iter(seq)
        
        raise ValueError("Could not find the sequence in the mmCIF dictionary.")


    @staticmethod
    def translate_three2one_iter(residue: Iterable) -> str:
        """
        Translates a sequence of three-letter amino acid codes to one-letter codes.

        Parameters:
        - residue (Iterable): A sequence of three-letter amino acid codes.

        Returns:
        - str: A string containing the translated one-letter amino acid codes.
        """
        return ''.join([SequenceReader.THREE2ONE[res.upper()] for res in residue])
    

    def clusters_to_fasta(self, cluster_intervals: dict) -> dict:
        """
        Parse the sequence file and return the sequences corresponding to each cluster.

        Parameters:
        - cluster_intervals (dict): A dictionary where the keys are the cluster indices and the 
                                    values are the corresponding intervals. (from AFragmenter object)

        Returns:
        - dict: A dictionary where the keys are the cluster indices and the values are the corresponding protein sequences.

        Raises:
        - ValueError: If the cluster intervals are not defined.
        """
        parsed_sequences = {}
        for cluster_id, interval in cluster_intervals.items():
            cluster_sequence = []
            for i, (start, end) in enumerate(interval):
                if i > 0: # If it isn't the first interval, check for gaps
                    previous_end = interval[i-1][1]
                    if start > previous_end:
                        gap_size = start - previous_end - 1
                        cluster_sequence.extend(['-'] * gap_size)
                
                cluster_sequence.extend(self.sequence[start:end+1])
            parsed_sequences[cluster_id] = ''.join(cluster_sequence)
        return parsed_sequences


    @staticmethod
    def get_name_cif(cif_dict: dict) -> str:
        """
        Extracts the name of the protein from a CIF dictionary.

        Parameters:
        - cif_dict (dict): The CIF dictionary from Bio.PDB.MMCIF2Dict.

        Returns:
        - str: The name of the protein. If the entry ID is not found, an empty string is returned.
        """
        name = cif_dict.get('_entry.id', [''])[0]
        return SequenceReader.clean_name(name)
    

    @staticmethod
    def get_name_pdb_chain(chain) -> str:
        """
        Extracts the name of the protein from a chain object.

        Parameters:
        - chain: The chain object from Bio.PDB.Chain.

        Returns:
        - str: The name of the protein. If the chain ID is not found, an empty string is returned.
        """
        try:
            name = chain.full_id[0]
        except AttributeError:
            name = ''
        return SequenceReader.clean_name(name)
    

    @staticmethod
    def get_name_pdb(structure) -> str:
        """
        Extracts the name of the protein from a PDB structure.

        Parameters:
        - structure: The PDB structure from Bio.PDB.PDBParser.

        Returns:
        - str: The name of the protein. If molecule and name are not found, an empty string is returned.
        """
        name = structure.header.get('compound', {}).get('1', {}).get('molecule', '')
        if not name:
            name = structure.header.get('name', '')
        print(name)
        return SequenceReader.clean_name(name)
    

    @staticmethod
    def get_name_fasta(record: SeqIO.SeqRecord) -> str:
        """
        Extracts the name of the protein from a FASTA record.

        Parameters:
        - record: The FASTA record from Bio.SeqIO.

        Returns:
        - str: The record ID cleaned using the clean_name method.
        """
        return SequenceReader.clean_name(record.id)


    @staticmethod
    def clean_name(name: str) -> str:
        """
        Cleans the name of a protein.
        Removes the 'AF-' prefix and the '-F1' suffix from the name.
        Converts spaces to underscores.

        Parameters:
        - name (str): The name of the protein.

        Returns:
        - str: The cleaned name.
        """
        name = re.sub(r'^AF-', '', name)
        name = re.sub(r'-F\d+$', '', name)
        return '_'.join(name.split())