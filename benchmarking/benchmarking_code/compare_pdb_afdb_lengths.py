import asyncio
import aiohttp
import json
import time
import logging
import re
import warnings
from typing import Optional, Dict, Tuple
from Bio import PDB
from io import StringIO
from functools import lru_cache

warnings.filterwarnings('ignore')

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class ProteinCoverageAnalyzer:
    """Optimized protein processor with concurrent requests and caching."""
    
    def __init__(self, max_concurrent: int = 20, request_delay: float = 0.05, min_plddt: float = 0.0):
        self.max_concurrent = max_concurrent
        self.request_delay = request_delay
        self.min_plddt = min_plddt  # Minimum pLDDT threshold
        self.session = None
        self.semaphore = None
        self.uniprot_cache = {}
        self.pdb_cache = {}
        self.alphafold_cache = {}
    
    async def __aenter__(self):
        connector = aiohttp.TCPConnector(limit=self.max_concurrent, limit_per_host=10)
        timeout = aiohttp.ClientTimeout(total=30)
        self.session = aiohttp.ClientSession(
            connector=connector,
            timeout=timeout,
            headers={'User-Agent': 'ProteinCoverageAnalyzer/1.0'}
        )
        self.semaphore = asyncio.Semaphore(self.max_concurrent)
        return self
    
    async def __aexit__(self, exc_type, exc_val, exc_tb):
        if self.session:
            await self.session.close()
    
    async def fetch_with_retry(self, url: str, max_retries: int = 3) -> Optional[str]:
        """Fetch URL with retry logic and rate limiting."""
        async with self.semaphore:
            for attempt in range(max_retries):
                try:
                    await asyncio.sleep(self.request_delay)  # Rate limiting
                    async with self.session.get(url) as response:
                        if response.status == 200:
                            return await response.text()
                        elif response.status == 404:
                            return None
                        else:
                            logger.warning(f"HTTP {response.status} for {url} (attempt {attempt + 1})")
                except Exception as e:
                    logger.warning(f"Request failed for {url} (attempt {attempt + 1}): {e}")
                    if attempt < max_retries - 1:
                        await asyncio.sleep(2 ** attempt)
            return None
    
    async def get_uniprot_id(self, pdb_id: str, chain_id: str) -> Optional[str]:
        """Get UniProt ID with caching."""
        cache_key = f"{pdb_id}:{chain_id}"
        if cache_key in self.uniprot_cache:
            return self.uniprot_cache[cache_key]
        
        pdb_id_lower = pdb_id.lower()
        url = f"https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/{pdb_id_lower}"
        
        try:
            response_text = await self.fetch_with_retry(url)
            if not response_text:
                self.uniprot_cache[cache_key] = None
                return None
            
            data = json.loads(response_text)
            uniprot_mappings = data.get(pdb_id_lower, {}).get('UniProt', {})
            
            for uniprot_id, mapping_data in uniprot_mappings.items():
                mappings = mapping_data.get('mappings', [])
                for mapping in mappings:
                    if mapping.get('chain_id') == chain_id:
                        self.uniprot_cache[cache_key] = uniprot_id
                        return uniprot_id
            
            self.uniprot_cache[cache_key] = None
            return None
            
        except Exception as e:
            logger.error(f"Error fetching UniProt ID for {pdb_id}:{chain_id}: {e}")
            self.uniprot_cache[cache_key] = None
            return None
    
    async def get_pdb_chain_length(self, pdb_id: str, chain_id: str) -> Optional[int]:
        """Get PDB chain length with caching."""
        cache_key = f"{pdb_id}:{chain_id}"
        if cache_key in self.pdb_cache:
            return self.pdb_cache[cache_key]
        
        url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
        
        try:
            pdb_data = await self.fetch_with_retry(url)
            if not pdb_data or len(pdb_data) < 100:
                self.pdb_cache[cache_key] = None
                return None
            
            # Parse PDB structure
            parser = PDB.PDBParser(QUIET=True)
            structure = parser.get_structure(pdb_id, StringIO(pdb_data))
            
            for model in structure:
                for chain in model:
                    if chain.id == chain_id:
                        residues = [res for res in chain if res.id[0] == ' ']
                        if residues:
                            length = len(residues)
                            self.pdb_cache[cache_key] = length
                            return length
            
            self.pdb_cache[cache_key] = None
            return None
            
        except Exception as e:
            logger.error(f"Error fetching PDB chain length for {pdb_id}:{chain_id}: {e}")
            self.pdb_cache[cache_key] = None
            return None
    
    async def get_alphafold_length(self, uniprot_id: str) -> Optional[int]:
        """Get AlphaFold length with caching and pLDDT check."""
        if uniprot_id in self.alphafold_cache:
            return self.alphafold_cache[uniprot_id]
        
        # Try AlphaFold v4 first, then v3
        urls = [
            f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v4.pdb",
            f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v3.pdb"
        ]
        
        for url in urls:
            try:
                pdb_data = await self.fetch_with_retry(url)
                if pdb_data:
                    length, avg_plddt = self.parse_alphafold_length_and_plddt(pdb_data)
                    if length and avg_plddt and avg_plddt >= self.min_plddt:
                        logger.debug(f"AlphaFold structure for {uniprot_id} passed pLDDT check: {avg_plddt:.2f}")
                        self.alphafold_cache[uniprot_id] = length
                        return length
                    elif length and avg_plddt:
                        logger.debug(f"AlphaFold structure for {uniprot_id} failed pLDDT check: {avg_plddt:.2f} < {self.min_plddt}")
            except Exception as e:
                logger.debug(f"Error fetching AlphaFold for {uniprot_id}: {e}")
                continue
        
        self.alphafold_cache[uniprot_id] = None
        return None
    
    def parse_alphafold_length_and_plddt(self, pdb_data: str) -> Tuple[Optional[int], Optional[float]]:
        """Parse AlphaFold PDB data to get length and average pLDDT."""
        try:
            lines = pdb_data.split('\n')
            atom_lines = [line for line in lines if line.startswith('ATOM')]
            
            if not atom_lines:
                return None, None
            
            # Get residue numbers and pLDDT scores from CA atoms
            residue_data = {}
            for line in atom_lines:
                if len(line) >= 80:  # Ensure line is long enough
                    atom_type = line[12:16].strip()
                    if atom_type == 'CA':  # Only process CA atoms
                        res_num_str = line[22:26].strip()
                        if res_num_str and res_num_str.isdigit():
                            res_num = int(res_num_str)
                            # pLDDT score is stored in the B-factor column (positions 60-66)
                            try:
                                plddt_score = float(line[60:66].strip())
                                residue_data[res_num] = plddt_score
                            except ValueError:
                                logger.debug(f"Could not parse pLDDT score from line: {line[60:66]}")
                                continue
            
            if not residue_data:
                return None, None
            
            # Calculate length and average pLDDT
            length = max(residue_data.keys()) if residue_data else None
            avg_plddt = sum(residue_data.values()) / len(residue_data) if residue_data else None
            
            return length, avg_plddt
            
        except Exception as e:
            logger.error(f"Error parsing AlphaFold PDB data: {e}")
            return None, None
    
    def parse_alphafold_length(self, pdb_data: str) -> Optional[int]:
        """Parse AlphaFold PDB data to get length (legacy method for backward compatibility)."""
        length, _ = self.parse_alphafold_length_and_plddt(pdb_data)
        return length
    
    async def get_alphafold_info(self, uniprot_id: str) -> Optional[Dict]:
        """Get AlphaFold structure info including length and pLDDT statistics."""
        if uniprot_id in self.alphafold_cache:
            # If cached, return basic info
            return {"length": self.alphafold_cache[uniprot_id]}
        
        # Try AlphaFold v4 first, then v3
        urls = [
            f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v4.pdb",
            f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v3.pdb"
        ]
        
        for url in urls:
            try:
                pdb_data = await self.fetch_with_retry(url)
                if pdb_data:
                    length, avg_plddt = self.parse_alphafold_length_and_plddt(pdb_data)
                    if length and avg_plddt:
                        info = {
                            "length": length,
                            "avg_plddt": avg_plddt,
                            "passes_plddt_check": avg_plddt >= self.min_plddt,
                            "plddt_threshold": self.min_plddt
                        }
                        
                        # Only cache if it passes the pLDDT check
                        if avg_plddt >= self.min_plddt:
                            self.alphafold_cache[uniprot_id] = length
                        
                        return info
            except Exception as e:
                logger.debug(f"Error fetching AlphaFold for {uniprot_id}: {e}")
                continue
        
        return None

'''
class ProteinCoverageAnalyzer:
    """Original analyzer class for backward compatibility."""
    
    def __init__(self, retry_attempts: int = 3, delay_between_requests: float = 0.2):
        self.retry_attempts = retry_attempts
        self.delay_between_requests = delay_between_requests
        # Keep the original implementation for compatibility
        
    def calculate_coverage(self, pdb_id: str, chain_id: str):
        """Original method - kept for backward compatibility."""
        # Implementation would go here if needed
        pass
'''

@lru_cache(maxsize=128)
def cached_pdb_to_uniprot(pdb_id: str, chain_id: str) -> Optional[str]:
    """Cached version of PDB to UniProt ID conversion."""
    return pdb_to_uniprot(pdb_id, chain_id)

def pdb_to_uniprot(pdb_id: str, chain_id: str) -> Optional[str]:
    """Convert PDB ID to UniProt ID using PDB API."""
    import requests
    
    pdb_id = pdb_id.lower()
    
    def fetch_pdb_info(pdb_id):
        time.sleep(0.1)
        url = f"https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/{pdb_id}"
        response = requests.get(url)
        if response.status_code == 200:
            return response.json()
        else:
            print(f"Error fetching data for {pdb_id}: {response.status_code}")
            return None
        
    fetch_result = fetch_pdb_info(pdb_id)
    if not fetch_result:
        return None
        
    for id in fetch_result.get(pdb_id).get('UniProt').keys():
        mapping = fetch_result.get(pdb_id).get('UniProt').get(id).get('mappings')
        for m in mapping:
            if m.get('chain_id') == chain_id:
                return id

def is_valid_uniprot_pattern(uniprot_id: str) -> bool:
    """Quick check if UniProt ID matches the expected pattern."""
    return bool(re.match(r'^[A-Z0-9]{6,10}$', uniprot_id))