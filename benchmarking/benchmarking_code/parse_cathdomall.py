"""
CATH Domall File Parser

Parses CATH Domall files to extract PDB ID, chain ID, and domain strings
in chainsaw format (e.g., '0-9_21-30,10-20').

Example usage:
    parser = CathDomallParser()
    domains = parser.parse_file('CathDomall.v4.4.0')
    
Example output:
    domains = {
        '101m': {'A': '0-153'},
        '102l': {'A': '1-162'},
        '102m': {'A': '0-153'},
        '108l': {'A': '1-162'},
        '108m': {'A': '0-153'},
        '109l': {'A': '1-162'},
        '109m': {'A': '0-153'},
        '10gs': {'A': '2-78_187-208,79-186', 'B': '2-78_187-208,79-186'},
        '10mh': {'A': '1-186_285-327,187-284'},
        '1amg': {'A': '1-360,362-417'}
        }
"""

from typing import List, Dict, Optional


class CathDomallParser:
    def __init__(self):
        self.domains = []
    
    def parse_line(self, line: str) -> Optional[Dict]:
        """Parse a single line from the CATH Domall file."""
        if line.startswith('#') or not line.strip():
            return None

        parts = line.strip().split()

        # Find Dxx and Fxx tokens
        d_idx = next((i for i, p in enumerate(parts) if p.startswith('D') and len(p) == 3), None)
        f_idx = next((i for i, p in enumerate(parts) if p.startswith('F') and len(p) == 3), None)

        if d_idx is None or f_idx is None or f_idx != d_idx + 1:
            return None  # Malformed line

        chain_name = parts[0].strip()
        num_domains = int(parts[d_idx][1:])  # e.g., 'D02' -> 2
        num_fragments = int(parts[f_idx][1:])  # e.g., 'F00' -> 0

        # Extract PDB ID and chain ID from chain name
        if len(chain_name) == 4:
            pdb_id = chain_name
            chain_id = "A"
        elif len(chain_name) == 5:
            pdb_id = chain_name[:4]
            chain_id = chain_name[4:]
        else:
            pdb_id = chain_name[:4]
            chain_id = chain_name[4:] if len(chain_name) > 4 else "A"

        # Everything after the Fxx field is domain data
        domain_data_parts = parts[f_idx + 1:]

        domains = self._parse_domains(domain_data_parts, num_domains, num_fragments)

        return {
            'pdb_id': pdb_id,
            'chain_id': chain_id,
            'chain_name': chain_name,
            'num_domains': num_domains,
            'num_fragments': num_fragments,
            'domains': domains
        }
    
    def _parse_domains(self, parts: List[str], num_domains: int, num_fragments: int) -> List[Dict]:
        """Parse domain segment information from the line parts."""
        domains = []
        i = 0
        
        # Remove fragment information first (everything from first '(' onwards)
        clean_parts = []
        for part in parts:
            if '(' in part:
                break
            clean_parts.append(part)
        
        #print(f"DEBUG: Parsing {num_domains} domains from clean parts: {clean_parts}")
        
        for domain_idx in range(num_domains):
            if i >= len(clean_parts):
                #print(f"DEBUG: Ran out of parts at domain {domain_idx + 1}")
                break
                
            # Get number of segments for this domain
            try:
                num_segments = int(clean_parts[i])
                #print(f"DEBUG: Domain {domain_idx + 1} has {num_segments} segments, starting at index {i}")
                i += 1
            except (ValueError, IndexError):
                #print(f"DEBUG: Could not parse segment count at index {i}")
                break
            
            segments = []
            for seg_idx in range(num_segments):
                if i + 5 < len(clean_parts):
                    try:
                        # Parse segment: chain_id start_pos insert_code chain_id end_pos insert_code
                        chain_start = clean_parts[i]
                        start_pos = int(clean_parts[i + 1])
                        start_insert = clean_parts[i + 2] if clean_parts[i + 2] != '-' else ''
                        chain_end = clean_parts[i + 3]
                        end_pos = int(clean_parts[i + 4])
                        end_insert = clean_parts[i + 5] if clean_parts[i + 5] != '-' else ''
                        
                        #print(f"DEBUG: Segment {seg_idx + 1}: {chain_start} {start_pos}{start_insert} to {chain_end} {end_pos}{end_insert}")
                        
                        segments.append({
                            'start_pos': start_pos,
                            'start_insert': start_insert,
                            'end_pos': end_pos,
                            'end_insert': end_insert
                        })
                        i += 6
                    except (ValueError, IndexError) as e:
                        #print(f"DEBUG: Error parsing segment {seg_idx + 1}: {e}")
                        break
                else:
                    #print(f"DEBUG: Not enough parts for segment {seg_idx + 1} at index {i}")
                    break
            
            domains.append({
                'domain_idx': domain_idx + 1,
                'segments': segments
            })
            #print(f"DEBUG: Completed domain {domain_idx + 1} with {len(segments)} segments")
        
        return domains
    
    def _format_domain_string(self, segments: List[Dict]) -> str:
        """Convert domain segments to chainsaw format string."""
        segment_strings = []
        
        for segment in segments:
            start = segment['start_pos']
            end = segment['end_pos']
            
            # Add insert codes if present
            start_str = str(start)
            if segment['start_insert']:
                start_str += segment['start_insert']
            
            end_str = str(end)
            if segment['end_insert']:
                end_str += segment['end_insert']
            
            segment_strings.append(f"{start_str}-{end_str}")
        
        return '_'.join(segment_strings)
    
    def parse_file(self, filename: str) -> List[Dict]:
        """Parse entire CATH Domall file."""
        results = []
        
        with open(filename, 'r') as f:
            for line_num, line in enumerate(f, 1):
                try:
                    parsed = self.parse_line(line)
                    if parsed:
                        # Generate domain strings for each domain
                        for domain in parsed['domains']:
                            domain_string = self._format_domain_string(domain['segments'])
                            
                            results.append({
                                'pdb_id': parsed['pdb_id'],
                                'chain_id': parsed['chain_id'],
                                'domain_string': domain_string,
                                'domain_idx': domain['domain_idx'],
                                'chain_name': parsed['chain_name'],
                                'num_domains': parsed['num_domains']
                            })
                
                except Exception as e:
                    print(f"Error parsing line {line_num}: {e}")
                    print(f"Line content: {line.strip()}")
        
        return results
    
    def parse_string(self, content: str) -> Dict[str, Dict[str, str]]:
        """Parse CATH Domall content from string and return nested dictionary."""
        result_dict = {}

        for line_num, line in enumerate(content.split('\n'), 1):
            try:
                parsed = self.parse_line(line)
                if parsed:
                    pdb_id = parsed['pdb_id']
                    chain_id = parsed['chain_id']

                    for domain in parsed['domains']:
                        domain_string = self._format_domain_string(domain['segments'])

                        # Initialize nested dictionary
                        if pdb_id not in result_dict:
                            result_dict[pdb_id] = {}
                        if chain_id not in result_dict[pdb_id]:
                            result_dict[pdb_id][chain_id] = []

                        result_dict[pdb_id][chain_id].append(domain_string)
            except Exception as e:
                print(f"Error parsing line {line_num}: {e}")
                print(f"Line content: {line.strip()}")

        # Join domain strings with ','
        for pdb_id in result_dict:
            for chain_id in result_dict[pdb_id]:
                result_dict[pdb_id][chain_id] = ','.join(result_dict[pdb_id][chain_id])

        return result_dict
    
    def parse_entry(self, entry: str) -> Optional[Dict[str, str]]:
        """Parse a single CATH entry string."""
        try:
            parsed = self.parse_line(entry)
            if parsed:
                pdb_id = parsed['pdb_id']
                chain_id = parsed['chain_id']
                domain_string = ','.join(
                    self._format_domain_string(domain['segments']) for domain in parsed['domains']
                )
                return {
                    'pdb_id': pdb_id,
                    'chain_id': chain_id,
                    'domain_string': domain_string
                }
        except Exception as e:
            print(f"Error parsing entry: {e}")
            return None
        return None
    