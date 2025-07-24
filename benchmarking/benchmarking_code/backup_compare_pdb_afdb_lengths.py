import requests
from Bio import PDB
from io import StringIO
from afragmenter import fetch_afdb_data
from functools import lru_cache
import time
import logging
from typing import Optional
import re
import warnings
warnings.filterwarnings('ignore')

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class ProteinCoverageAnalyzer:
    """Analyzes PDB coverage of AlphaFold structures with robust error handling."""
    
    def __init__(self, retry_attempts: int = 3, delay_between_requests: float = 0.2):
        self.retry_attempts = retry_attempts
        self.delay_between_requests = delay_between_requests
        self.session = requests.Session()
        self.session.headers.update({'User-Agent': 'ProteinCoverageAnalyzer/1.0'})
        
    def get_pdb_chain_length(self, pdb_id: str, chain_id: str) -> Optional[int]:
        """Get the length of a specific chain from a PDB structure with retry logic."""
        for attempt in range(self.retry_attempts):
            try:
                # Download PDB file
                url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
                response = self.session.get(url, timeout=30)
                response.raise_for_status()
                
                # Check if we got a valid PDB file
                if len(response.text) < 100:
                    logger.warning(f"PDB file for {pdb_id} seems too small")
                    return None
                
                # Parse PDB structure
                parser = PDB.PDBParser(QUIET=True)
                structure = parser.get_structure(pdb_id, StringIO(response.text))
                
                # Get the specific chain
                for model in structure:
                    for chain in model:
                        if chain.id == chain_id:
                            # Count residues (excluding water and other non-protein atoms)
                            residues = [res for res in chain if res.id[0] == ' ']
                            if residues:
                                return len(residues)
                
                logger.warning(f"Chain {chain_id} not found in {pdb_id}")
                return None
                
            except requests.exceptions.RequestException as e:
                logger.warning(f"Request error for {pdb_id} (attempt {attempt + 1}): {e}")
                if attempt < self.retry_attempts - 1:
                    time.sleep(2 ** attempt)  # Exponential backoff
                    continue
            except Exception as e:
                logger.error(f"Unexpected error processing {pdb_id}:{chain_id} - {e}")
                return None
        
        return None
    
    def get_alphafold_length(self, uniprot_id: str) -> Optional[int]:
        """Get the length of the AlphaFold structure with better error handling for both mmCIF and PDB formats."""
        for attempt in range(self.retry_attempts):
            try:
                json_data, structure_data = fetch_afdb_data(uniprot_id)
                
                if not structure_data:
                    logger.warning(f"No structure data returned for {uniprot_id}")
                    return None
                
                # Determine file format based on content
                lines = structure_data.split('\n')
                is_mmcif = any(line.startswith('data_') or line.startswith('_') for line in lines[:10])
                
                if is_mmcif:
                    logger.debug(f"Processing mmCIF format for {uniprot_id}")
                    return self._parse_mmcif_length(structure_data, uniprot_id)
                else:
                    logger.debug(f"Processing PDB format for {uniprot_id}")
                    return self._parse_pdb_length(structure_data, uniprot_id)
                    
            except requests.exceptions.HTTPError as e:
                if e.response.status_code == 404:
                    logger.info(f"AlphaFold structure not available for {uniprot_id} (404)")
                    return None
                else:
                    logger.warning(f"HTTP error for {uniprot_id} (attempt {attempt + 1}): {e}")
                    if attempt < self.retry_attempts - 1:
                        time.sleep(2 ** attempt)
                        continue
            except Exception as e:
                logger.error(f"Unexpected error fetching AlphaFold for {uniprot_id} - {e}")
                if attempt < self.retry_attempts - 1:
                    time.sleep(1)
                    continue
        
        return None

    def _parse_mmcif_length(self, mmcif_data: str, uniprot_id: str) -> Optional[int]:
        """Parse mmCIF format to get structure length."""
        try:
            # Method 1: Use BioPython MMCIF parser (most reliable)
            from Bio.PDB import MMCIFParser
            parser = MMCIFParser(QUIET=True)
            structure = parser.get_structure(uniprot_id, StringIO(mmcif_data))
            
            max_res_num = 0
            residue_count = 0
            for model in structure:
                for chain in model:
                    for residue in chain:
                        if residue.id[0] == ' ':  # Standard residue (not heteroatom)
                            residue_count += 1
                            res_num = residue.id[1]
                            if isinstance(res_num, int):
                                max_res_num = max(max_res_num, res_num)
            
            if residue_count > 0:
                logger.debug(f"mmCIF BioPython: {residue_count} residues, max number {max_res_num}")
                return max(residue_count, max_res_num)
                
        except Exception as e:
            logger.debug(f"mmCIF BioPython parsing failed: {e}")
        
        try:
            # Method 2: Manual parsing of mmCIF atom_site records
            lines = mmcif_data.split('\n')
            in_atom_site = False
            residue_numbers = set()
            
            for line in lines:
                line = line.strip()
                if line.startswith('_atom_site.'):
                    in_atom_site = True
                    continue
                elif line.startswith('_') and in_atom_site:
                    break
                elif in_atom_site and line and not line.startswith('#'):
                    # Parse atom site data
                    fields = line.split()
                    if len(fields) >= 8:  # Typical mmCIF atom_site has many fields
                        try:
                            # In mmCIF, residue number is typically in auth_seq_id field
                            # Field positions may vary, so we look for numeric values that could be residue numbers
                            for field in fields:
                                if field.isdigit():
                                    res_num = int(field)
                                    if 1 <= res_num <= 10000:  # Reasonable range for residue numbers
                                        residue_numbers.add(res_num)
                        except (ValueError, IndexError):
                            continue
            
            if residue_numbers:
                max_res = max(residue_numbers)
                logger.debug(f"mmCIF manual parsing: residues from {min(residue_numbers)} to {max_res}")
                return max_res
                
        except Exception as e:
            logger.debug(f"mmCIF manual parsing failed: {e}")
        
        return None

    def _parse_pdb_length(self, pdb_data: str, uniprot_id: str) -> Optional[int]:
        """Parse PDB format to get structure length."""
        lines = pdb_data.split('\n')
        
        try:
            # Method 1: Use BioPython PDB parser (most reliable)
            from Bio.PDB import PDBParser
            parser = PDBParser(QUIET=True)
            structure = parser.get_structure(uniprot_id, StringIO(pdb_data))
            
            max_res_num = 0
            residue_count = 0
            for model in structure:
                for chain in model:
                    for residue in chain:
                        if residue.id[0] == ' ':  # Standard residue (not heteroatom)
                            residue_count += 1
                            res_num = residue.id[1]
                            if isinstance(res_num, int):
                                max_res_num = max(max_res_num, res_num)
            
            if residue_count > 0:
                logger.debug(f"PDB BioPython: {residue_count} residues, max number {max_res_num}")
                return max(residue_count, max_res_num)
                
        except Exception as e:
            logger.debug(f"PDB BioPython parsing failed: {e}")
        
        try:
            # Method 2: Manual PDB parsing
            atom_lines = [line for line in lines if line.startswith('ATOM')]
            
            if not atom_lines:
                logger.warning(f"No ATOM lines found for {uniprot_id}")
                return None
            
            # Parse residue numbers from PDB format
            residue_numbers = set()
            for line in atom_lines:
                if len(line) >= 26:
                    # PDB format: residue number is in columns 23-26 (0-indexed: 22-25)
                    res_num_str = line[22:26].strip()
                    if res_num_str and res_num_str.isdigit():
                        res_num = int(res_num_str)
                        residue_numbers.add(res_num)
            
            if residue_numbers:
                max_res = max(residue_numbers)
                logger.debug(f"PDB manual parsing: residues from {min(residue_numbers)} to {max_res}")
                return max_res
            
            # Method 3: Count CA atoms as fallback
            ca_atoms = [line for line in atom_lines if len(line) >= 16 and line[12:16].strip() == 'CA']
            if ca_atoms:
                logger.debug(f"PDB CA count: {len(ca_atoms)} CA atoms")
                return len(ca_atoms)
                
        except Exception as e:
            logger.debug(f"PDB manual parsing failed: {e}")
        
        return None

    def backup_get_alphafold_length(self, uniprot_id: str) -> Optional[int]:
        """Get the length of the AlphaFold structure with better error handling."""
        for attempt in range(self.retry_attempts):
            try:
                json_data, pdb_data = fetch_afdb_data(uniprot_id)
                
                if not pdb_data:
                    logger.warning(f"No PDB data returned for {uniprot_id}")
                    return None
                
                
                # Debug: Print first few lines to understand format
                lines = pdb_data.split('\n')
                logger.debug(f"First 3 ATOM lines for {uniprot_id}:")
                atom_count = 0
                for line in lines:
                    if line.startswith('ATOM') and atom_count < 3:
                        logger.debug(f"  '{line}'")
                        atom_count += 1
                
                # Count lines starting with ATOM (excluding HETATM)
                atom_lines = [line for line in lines if line.startswith('ATOM')]
                
                if not atom_lines:
                    logger.warning(f"No ATOM lines found for {uniprot_id}")
                    return None
                
                logger.debug(f"Found {len(atom_lines)} ATOM lines for {uniprot_id}")
                
                # Method 1: Parse residue numbers from PDB format
                residue_numbers = set()
                for line in atom_lines:
                    try:
                        # Ensure line is long enough
                        if len(line) < 26:
                            continue
                        
                        # PDB format: residue number is in columns 23-26 (0-indexed: 22-25)
                        res_num_str = line[22:26].strip()
                        if res_num_str and res_num_str.isdigit():
                            res_num = int(res_num_str)
                            residue_numbers.add(res_num)
                    except (ValueError, IndexError):
                        continue
                
                if residue_numbers:
                    max_res = max(residue_numbers)
                    logger.debug(f"Method 1: Found residue numbers from {min(residue_numbers)} to {max_res}")
                    return max_res
                
                # Method 2: Count CA atoms (one per residue)
                ca_atoms = []
                for line in atom_lines:
                    if len(line) >= 16 and line[12:16].strip() == 'CA':
                        ca_atoms.append(line)
                
                if ca_atoms:
                    logger.debug(f"Method 2: Found {len(ca_atoms)} CA atoms for {uniprot_id}")
                    return len(ca_atoms)
                
                # Method 3: Use BioPython parser as fallback
                try:
                    parser = PDB.PDBParser(QUIET=True)
                    structure = parser.get_structure(uniprot_id, StringIO(pdb_data))
                    
                    max_res_num = 0
                    residue_count = 0
                    for model in structure:
                        for chain in model:
                            for residue in chain:
                                if residue.id[0] == ' ':  # Standard residue
                                    residue_count += 1
                                    res_num = residue.id[1]
                                    max_res_num = max(max_res_num, res_num)
                    
                    if residue_count > 0:
                        logger.debug(f"Method 3: BioPython found {residue_count} residues, max number {max_res_num}")
                        return max(residue_count, max_res_num)
                        
                except Exception as parser_error:
                    logger.debug(f"BioPython parsing failed for {uniprot_id}: {parser_error}")
                
                # Method 4: Use JSON data if available
                if json_data and isinstance(json_data, dict):
                    if 'uniprotEnd' in json_data:
                        length = json_data['uniprotEnd']
                        logger.debug(f"Method 4: Using JSON uniprotEnd = {length}")
                        return length
                    elif 'uniprotStart' in json_data and 'uniprotEnd' in json_data:
                        length = json_data['uniprotEnd'] - json_data['uniprotStart'] + 1
                        logger.debug(f"Method 4: Calculated length from JSON = {length}")
                        return length
                
                logger.warning(f"All parsing methods failed for {uniprot_id}")
                return None
                    
            except requests.exceptions.HTTPError as e:
                if e.response.status_code == 404:
                    logger.info(f"AlphaFold structure not available for {uniprot_id} (404)")
                    return None
                else:
                    logger.warning(f"HTTP error for {uniprot_id} (attempt {attempt + 1}): {e}")
                    if attempt < self.retry_attempts - 1:
                        time.sleep(2 ** attempt)
                        continue
            except Exception as e:
                logger.error(f"Unexpected error fetching AlphaFold for {uniprot_id} - {e}")
                if attempt < self.retry_attempts - 1:
                    time.sleep(1)
                    continue
    
        return None

        
    def calculate_coverage(self, pdb_id: str, chain_id: str) -> Optional[float]:
        """Calculate PDB coverage of AlphaFold structures."""
        uniprot_id = cached_pdb_to_uniprot(pdb_id, chain_id)

        if not uniprot_id:
            logger.error(f"Failed to find UniProt ID for {pdb_id}:{chain_id}")
            return None
        
        if not is_valid_uniprot_pattern(uniprot_id):
            logger.error(f"Invalid UniProt ID format for {pdb_id}:{chain_id} -> {uniprot_id}")
            return None

        time.sleep(self.delay_between_requests)  # Respect API rate limits
        pdb_length = self.get_pdb_chain_length(pdb_id, chain_id)
        alphafold_length = self.get_alphafold_length(uniprot_id)
        
        coverage = None
        #pass_filter = False
        
        if pdb_length is not None and alphafold_length is not None:
            coverage = pdb_length / alphafold_length
            #passes_filter = coverage >= min_coverage
            logger.debug(f"  PDB: {pdb_length}, AlphaFold: {alphafold_length}, Coverage: {coverage:.3f}")
        else:
            logger.warning(f"  Failed to get lengths (PDB: {pdb_length}, AlphaFold: {alphafold_length})")
            
        if coverage == None:
            logger.error(f"Coverage calculation failed for {pdb_id}:{chain_id}")
            return None, None, None # Or should we return None?
            
        return coverage, pdb_length, alphafold_length

    '''
    def preprocess_dataframe(self, df: pd.DataFrame) -> pd.DataFrame:
        """Preprocess the dataframe with better validation."""
        logger.info(f"Starting preprocessing with {len(df)} rows")
        
        # Filter out rows with NaN uniprot_id
        initial_count = len(df)
        df = df.dropna(subset=['uniprot_id'])
        logger.info(f"Removed {initial_count - len(df)} rows with missing uniprot_id")
        
        # Filter out empty string uniprot_id
        df = df[df['uniprot_id'].astype(str).str.len() > 0]
        logger.info(f"Remaining after filtering empty uniprot_id: {len(df)}")
        
        # Validate UniProt ID format (basic validation)
        valid_uniprot_pattern = df['uniprot_id'].astype(str).str.match(r'^[A-Z0-9]{6,10}$')
        df = df[valid_uniprot_pattern]
        logger.info(f"Remaining after UniProt ID validation: {len(df)}")
        
        # Deduplicate
        before_dedup = len(df)
        df = df.drop_duplicates(subset=['uniprot_id'], keep='first').copy()
        logger.info(f"Removed {before_dedup - len(df)} duplicates, remaining: {len(df)}")
        
        return df
    '''

@lru_cache(maxsize=128)
def cached_pdb_to_uniprot(pdb_id: str, chain_id: str) -> Optional[str]:
    """Cached version of PDB to UniProt ID conversion to speed up repeated calls."""
    return pdb_to_uniprot(pdb_id, chain_id)

# Convert PDB id to UniProt ID
def pdb_to_uniprot(pdb_id: str, chain_id: str) -> Optional[str]:
    """Convert PDB ID to UniProt ID using PDB API."""
    
    pdb_id = pdb_id.lower()  # Ensure PDB ID is in lowercase
    
    def fetch_pdb_info(pdb_id):
        time.sleep(0.1)  # To avoid hitting the API too hard
        url = f"https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/{pdb_id}"
        response = requests.get(url)
        if response.status_code == 200:
            return response.json()
        else:
            print(f"Error fetching data for {pdb_id}: {response.status_code}")
            return None
        
    fetch_result  = fetch_pdb_info(pdb_id)
    if not fetch_result:
        return None
    for id in fetch_result.get(pdb_id).get('UniProt').keys():
        mapping = fetch_result.get(pdb_id).get('UniProt').get(id).get('mappings')
        for m in mapping:
            if m.get('chain_id') == chain_id:
                return id # = UniProt ID
    

def is_valid_uniprot_pattern(uniprot_id: str) -> bool:
    """Quick check if UniProt ID matches the expected pattern."""
    return bool(re.match(r'^[A-Z0-9]{6,10}$', uniprot_id))