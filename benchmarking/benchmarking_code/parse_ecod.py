import re
from collections import defaultdict


class ECODParser:
    """
    Parser for ECOD (Evolutionary Classification of Domains) data files.
    
    Efficiently processes large ECOD files by streaming through the data
    and building domain dictionaries with metadata.
    """
    
    def __init__(self, file_path=None):
        """
        Initialize ECOD parser.
        
        Args:
            file_path (str, optional): Path to ECOD data file
        """
        self.file_path = file_path
        self.domains = {}
        self.metadata = {}
        
        if file_path:
            self.parse_file(file_path)
    
    def parse_file(self, file_path):
        """
        Parse ECOD data file efficiently using streaming.
        
        Args:
            file_path (str): Path to the ECOD data file
        """
        self.file_path = file_path
        self.domains = defaultdict(lambda: defaultdict(list))
        self.metadata = {}

        # Read file and remove duplicates
        with open(file_path, 'r') as file:
            lines = file.readlines()  # Read all lines into a list
        unique_lines = list(set(lines))  # Remove duplicates by converting to a set and back to a list

        # Parse header and process data lines
        for line in unique_lines:
            line = line.strip()
            if line.startswith('#'):
                self._parse_header_line(line)
            elif line:  # Process non-comment lines
                self._process_data_line(line)

        # Convert defaultdict to regular dict and merge domains
        self._finalize_domains()
    
    def parse_string(self, data_string):
        """
        Parse ECOD data from string.
        
        Args:
            data_string (str): ECOD data as string
        """
        self.domains = defaultdict(lambda: defaultdict(list))
        self.metadata = {}
        
        lines = data_string.strip().split('\n')
        
        # Parse header and data
        for line in lines:
            line = line.strip()
            if line.startswith('#'):
                self._parse_header_line(line)
            elif line:
                self._process_data_line(line)
                #print("Current domains state:")
                #print(self.domains)  # For some reason, nothing is printed here

        # Convert defaultdict to regular dict and merge domains
        self._finalize_domains()
    
    def _parse_header_line(self, line):
        """Parse header line for metadata."""
        if line.startswith('#ECOD version'):
            self.metadata['ecod_version'] = line.split('version', 1)[1].strip()
        elif line.startswith('#Domain list version'):
            self.metadata['domain_list_version'] = line.split('version', 1)[1].strip()
        elif 'Grishin Lab' in line:
            self.metadata['source'] = line.replace('#', '').strip()
        elif line.startswith('#/data/ecod/'):
            self.metadata['database_path'] = line.replace('#', '').strip()
    
    def _process_data_line(self, line):
        """Process a single data line."""
        columns = line.split()
        if len(columns) < 7:
            print(f"Skipping line due to insufficient columns: {line}")
            return
        
        source_id = self._clean_source_id(columns[0])
        d_type = columns[2]
        range_info = columns[6]
        
        # Parse chain and range information
        if d_type == 'PDB':
            if ':' in range_info:
                chain_id = range_info.split(':', 1)[0]
                range_part = range_info.split(':', 1)[1]
            else:
                chain_id = 'A'
                range_part = range_info
        else:  # AFDB
            chain_id = 'A'
            range_part = range_info
        
        # Parse range information
        domain_ranges = []
        for range_segment in range_part.split(','):
            range_segment = range_segment.strip()
            if range_segment:
                if ':' in range_segment:
                    range_segment = range_segment.split(':', 1)[1]
                domain_ranges.append(range_segment)
        
        # Join ranges with underscore for multiple segments
        domain_string = '_'.join(domain_ranges)
        #print(f"Source ID: {source_id}, Chain ID: {chain_id}, Domain String: {domain_string}")
        self.domains[source_id][chain_id].append(domain_string)
    
    def _clean_source_id(self, source_id):
        """Remove _Fx suffix from source_id where x is any number."""
        return re.sub(r'_F\d+$', '', source_id)
    
    def _finalize_domains(self):
        """Convert defaultdict to regular dict and merge domains with pipe separator."""
        final_domains = {}
        for source_id, chains in self.domains.items():
            final_domains[source_id] = {}
            for chain_id, domain_list in chains.items():
                final_domains[source_id][chain_id] = ','.join(domain_list)
        
        self.domains = final_domains
    
    def get_domain_count(self):
        """Get total number of domain entries."""
        count = 0
        for chains in self.domains.values():
            for domain_string in chains.values():
                count += len(domain_string.split(','))
        return count
    
    def get_protein_count(self):
        """Get total number of unique proteins."""
        return len(self.domains)
    
    def get_chain_count(self):
        """Get total number of chains across all proteins."""
        count = 0
        for chains in self.domains.values():
            count += len(chains)
        return count


def order_domains_by_start(parser, by_size=False):
    """
    Order domains in the parser results by start coordinate or size.
    
    Args:
        parser (ECODParser): Parsed ECOD data
        by_size (bool): If True, order by size (largest first), 
                       otherwise by start coordinate (smallest first)
    
    Returns:
        dict: Reordered domain dictionary
    """
    def get_domain_start(domain_string):
        """Get the start coordinate of the first range in a domain."""
        first_range = domain_string.split('_')[0]
        if '-' in first_range:
            try:
                return int(first_range.split('-')[0])
            except ValueError:
                return float('inf')
        return float('inf')
    
    def get_domain_size(domain_string):
        """Calculate total size of a domain."""
        total_size = 0
        for segment in domain_string.split('_'):
            if '-' in segment:
                try:
                    start, end = segment.split('-')
                    total_size += int(end) - int(start) + 1
                except ValueError:
                    continue
        return total_size
    
    ordered_domains = {}
    
    for source_id, chains in parser.domains.items():
        ordered_domains[source_id] = {}
        for chain_id, domain_string in chains.items():
            domains = domain_string.split(',')
            
            if by_size:
                domains.sort(key=get_domain_size, reverse=True)
            else:
                domains.sort(key=get_domain_start)
            
            ordered_domains[source_id][chain_id] = ','.join(domains)
    
    return ordered_domains


# Example usage
if __name__ == "__main__":
    sample_data = """#/data/ecod/database_versions/v292;
#ECOD version develop292
#Domain list version 1.7
#Grishin Lab (http://prodata.swmed.edu/ecod)
#source_id	domain_id	d_type	t_id	f_id	ecod_uid	range
12as	e12asA1	PDB	314.1.1	314.1.1.7	000007022	A:4-330	MANUAL_REP
16vp	e16vpA1	PDB	816.1.1	816.1.1.1	000007584	A:1-303,A:349-356	MANUAL_REP
1a0a	e1a0aA1	PDB	105.1.1	105.1.1.1	000004281	A:1-63	MANUAL_REP
1a0i	e1a0iA2	PDB	206.1.3	206.1.3.1	000007114	A:1-239	MANUAL_REP
1a0p	e1a0pA1	PDB	101.1.8	101.1.8.1	000002807	A:109-237	MANUAL_REP
1a0p	e1a0pA3	PDB	186.1.1	186.1.1.1	000003058	A:1-98	MANUAL_REP
1a3q	e1a3qA2	PDB	11.1.5	11.1.5.3	000001869	A:1-184	MANUAL_REP
1a3x	e1a3xA3	PDB	7518.1.1	7518.1.1.1	000010312	A:367-500	MANUAL_REP
1a5t	e1a5tA3	PDB	138.1.1	138.1.1.3	000003384	A:208-330	MANUAL_REP
1a5t	e1a5tA4	PDB	148.1.3	148.1.3.39	001002401	A:166-207	MANUAL_REP
1a5t	e1a5tA5	PDB	2004.1.1	2004.1.1.187	001002400	A:1-165	MANUAL_REP
1a5z	e1a5zA1	PDB	279.1.1	279.1.1.1	000007483	A:141-312	MANUAL_REP
X8FP97_F1	nD1	AFDB	2003.1.1	2003.1.1.45	003961586	1-175	AUTOMATIC_NONREP
X8FP97_F1	nD2	AFDB	129.1.1	129.1.1.3	003961578	176-260	AUTOMATIC_NONREP
X8FP97_F1	nD3	AFDB	2003.1.12	2003.1.12.1	003961575	261-390	AUTOMATIC_NONREP
X8FPB7_F1	nD1	AFDB	850.1.1	850.1.1.1	003959490	26-75	AUTOMATIC_NONREP
X8FPB7_F1	nD2	AFDB	304.48.1	304.48.1.10	003959489	6-25,76-175	AUTOMATIC_NONREP
X8FPB7_F1	nD3	AFDB	102.1.1	102.1.1.23	003962409	176-235	AUTOMATIC_NONREP
X8FPB7_F1	nD4	AFDB	302.1.1	302.1.1.1	003959474	236-375	AUTOMATIC_NONREP
X8FPC3_F1	nD1	AFDB	2011.1.1	2011.1.1.8	003959498	21-220	AUTOMATIC_NONREP
X8FPC8_F1	nD1	AFDB	101.1.1	101.1.1.5	003962396	11-80	AUTOMATIC_NONREP
X8FPG6_F1	nD1	AFDB	304.11.1	304.11.1.1	003961871	26-95	AUTOMATIC_NONREP
X8FPG6_F1	nD3	AFDB	222.1.1	222.1.1.12	003961872	251-390	AUTOMATIC_NONREP
X8FPG9_F1	nD1	AFDB	298.3.1	298.3.1.2	003961177	1-115	AUTOMATIC_NONREP
X8FPH6_F1	nD2	AFDB	2003.1.1	2003.1.1.37	003961931	51-185	AUTOMATIC_NONREP
X8FPM3_F1	nD3	AFDB	267.1.1	267.1.1.2	003960063	231-365	AUTOMATIC_NONREP
X8FPQ0_F1	nD2	AFDB	633.11.1	633.11.1.1	003958054	46-150	AUTOMATIC_NONREP
X8FPQ7_F1	nD1	AFDB	7581.1.1	7581.1.1.2	003957633	1-170	AUTOMATIC_NONREP
X8FPT3_F1	nD2	AFDB	191.1.1	191.1.1.41	003960409	106-220	AUTOMATIC_NONREP
X8FPX5_F1	nD1	AFDB	2003.1.2	2003.1.2.5	003961503	6-165	AUTOMATIC_NONREP
X8FPZ0_F1	nD2	AFDB	2003.1.5	2003.1.5.21	003962080	156-320	AUTOMATIC_NONREP
X8FQ03_F1	nD1	AFDB	7581.1.1	7581.1.1.2	003962911	41-315	AUTOMATIC_NONREP
Z4YHZ5_F1	nD1	AFDB	3777.1.1	3777.1.1.1	003899303	21-270	AUTOMATIC_NONREP
Z4YHZ5_F1	nD2	AFDB	10.12.1	10.12.1.53	003899304	271-485	AUTOMATIC_NONREP"""
    
    # Parse the data
    parser = ECODParser()
    parser.parse_string(sample_data)
    
    print("Metadata:")
    for key, value in parser.metadata.items():
        print(f"  {key}: {value}")
    
    print(f"\nStats:")
    print(f"  Proteins: {parser.get_protein_count()}")
    print(f"  Chains: {parser.get_chain_count()}")
    print(f"  Domains: {parser.get_domain_count()}")
    
    print(f"\nSample domains (original order):")
    for source_id, chains in list(parser.domains.items())[:5]:
        print(f"  {source_id}: {chains}")
    
    # Order by start coordinate
    ordered_by_start = order_domains_by_start(parser)
    print(f"\nSample domains (ordered by start coordinate):")
    for source_id, chains in list(ordered_by_start.items())[:5]:
        print(f"  {source_id}: {chains}")
