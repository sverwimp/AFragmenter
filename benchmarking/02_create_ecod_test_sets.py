import random
import json
import asyncio
import time
from typing import List, Dict, Optional, Set
from benchmarking_code import ProteinCoverageAnalyzer
from benchmarking_code.parse_ecod import ECODParser, order_domains_by_start
import os

# Configuration
TARGET_PDB_DATAPOINTS = 500
TARGET_AFDB_DATAPOINTS = 500
BATCH_SIZE = 200
ALPHAFOLD_BATCH_SIZE = 100
MIN_COVERAGE = 0.9
MAX_COVERAGE = 1.1
MAX_INDEX_RATIO = 1.1
RANDOM_SEED = 42
MAX_CONCURRENT_REQUESTS = 20
REQUEST_DELAY = 0.05


def validate_chopping_indices(true_chopping: str, pdb_length: int, alphafold_length: int, 
                             pdb_id: str, chain_id: str) -> bool:
    """Validate that chopping indices are reasonable compared to protein lengths."""
    try:
        indexes = []
        for part in true_chopping.split(','):
            part = part.strip()
            for p in part.split('_'):
                p = p.strip()
                if p.isdigit():
                    indexes.append(int(p))
                else:
                    # Handle ranges like '1-10'
                    if '-' in p:
                        end = p.split('-')[-1]  # Take the end part
                        end = int(end.strip())
                        indexes.append(end)
                    else:
                        print(f"Unexpected part in true_chopping: {p}")
        
        if not indexes:
            print(f"No valid indexes found in true_chopping for {pdb_id}:{chain_id}")
            return False
            
        max_index = max(indexes)
        
        # Check if max index is too high compared to alphafold length
        if max_index > MAX_INDEX_RATIO * alphafold_length:
            print(f"Max index {max_index} for {pdb_id}:{chain_id} is too high compared to lengths: {pdb_length}, {alphafold_length}")
            return False
            
        return True
        
    except Exception as e:
        print(f"Error validating chopping indices for {pdb_id}:{chain_id}: {e}")
        return False

def create_pdb_chain_list(pdb_domains: Dict) -> List[tuple]:
    """Create a list of (pdb_id, chain_id) tuples from pdb_domains."""
    pdb_chain_list = []
    for pdb_id, chains in pdb_domains.items():
        for chain_id in chains.keys():
            pdb_chain_list.append((pdb_id, chain_id))
    return pdb_chain_list

def sample_unique_pdb_chains(pdb_domains: Dict, target_count: int) -> List[tuple]:
    """Sample unique PDB IDs and randomly select one chain from each."""
    random.seed(RANDOM_SEED)
    
    # Get all unique PDB IDs
    unique_pdb_ids = list(pdb_domains.keys())
    
    # Sample the required number of unique PDB IDs
    sampled_pdb_ids = random.sample(unique_pdb_ids, min(target_count, len(unique_pdb_ids)))
    
    # For each PDB ID, randomly select one chain
    sampled_pairs = []
    for pdb_id in sampled_pdb_ids:
        available_chains = list(pdb_domains[pdb_id].keys())
        selected_chain = random.choice(available_chains)
        sampled_pairs.append((pdb_id, selected_chain))
    
    return sampled_pairs

async def process_single_pdb_protein(processor: ProteinCoverageAnalyzer, 
                                    pdb_id: str, 
                                    chain_id: str,
                                    pdb_domains: Dict,
                                    seen_uniprot_ids: Set[str]) -> Optional[Dict]:
    """Process a single PDB protein asynchronously."""
    try:
        # Get UniProt ID
        uniprot_id = await processor.get_uniprot_id(pdb_id, chain_id)
        if not uniprot_id:
            return None
        
        # Check if UniProt ID is already in the dataset
        if uniprot_id in seen_uniprot_ids:
            print(f"Skipping PDB {pdb_id}:{chain_id} - UniProt ID {uniprot_id} already in dataset")
            return None
        
        # Get lengths concurrently
        pdb_length_task = processor.get_pdb_chain_length(pdb_id, chain_id)
        alphafold_length_task = processor.get_alphafold_length(uniprot_id)
        
        pdb_length, alphafold_length = await asyncio.gather(
            pdb_length_task, alphafold_length_task, return_exceptions=True
        )
        
        # Handle exceptions
        if isinstance(pdb_length, Exception):
            print(f"Error: PDB length error for {pdb_id}:{chain_id}: {pdb_length}")
            return None
        if isinstance(alphafold_length, Exception):
            print(f"Error: AlphaFold length error for {uniprot_id}: {alphafold_length}")
            return None
        
        if pdb_length is None or alphafold_length is None:
            return None
        
        # Calculate coverage
        coverage = pdb_length / alphafold_length
        
        # Check coverage bounds
        if coverage < MIN_COVERAGE or coverage > MAX_COVERAGE:
            return None
        
        # Get true chopping
        true_chopping = pdb_domains[pdb_id][chain_id]
        
        # Validate chopping indices
        if not validate_chopping_indices(true_chopping, pdb_length, alphafold_length, pdb_id, chain_id):
            return None
        
        return {
            'pdb_id': pdb_id,
            'pdb_length': pdb_length,
            'alphafold_length': alphafold_length,
            'chain_id': chain_id,
            'uniprot_id': uniprot_id,
            'coverage': round(coverage, 5),
            'true_chopping': true_chopping,
        }
        
    except Exception as e:
        print(f"Error processing PDB {pdb_id}:{chain_id}: {e}")
        return None

async def process_single_afdb_protein(processor: ProteinCoverageAnalyzer, 
                                     uniprot_id: str,
                                     chain_id: str,
                                     afdb_domains: Dict,
                                     seen_uniprot_ids: Set[str]) -> Optional[Dict]:
    """Process a single AlphaFold protein asynchronously."""
    try:
        # Check if UniProt ID is already in the dataset
        if uniprot_id in seen_uniprot_ids:
            print(f"Skipping AFDB {uniprot_id}:{chain_id} - UniProt ID {uniprot_id} already in dataset")
            return None
        
        # Get AlphaFold length
        alphafold_length = await processor.get_alphafold_length(uniprot_id)
        
        if alphafold_length is None:
            return None
        
        # Get true chopping
        true_chopping = afdb_domains[uniprot_id][chain_id]
        
        return {
            'uniprot_id': uniprot_id,
            'alphafold_length': alphafold_length,
            'chain_id': chain_id,
            'true_chopping': true_chopping,
        }
        
    except Exception as e:
        print(f"Error processing AlphaFold {uniprot_id}: {e}")
        return None

async def process_pdb_batch(processor: ProteinCoverageAnalyzer,
                           pdb_chain_pairs: List[tuple], 
                           pdb_domains: Dict,
                           seen_uniprot_ids: Set[str]) -> List[Dict]:
    """Process a batch of PDB proteins concurrently."""
    print(f"Processing PDB batch of {len(pdb_chain_pairs)} proteins...")
    
    tasks = []
    for pdb_id, chain_id in pdb_chain_pairs:
        task = process_single_pdb_protein(processor, pdb_id, chain_id, pdb_domains, seen_uniprot_ids)
        tasks.append(task)
    
    results = await asyncio.gather(*tasks, return_exceptions=True)
    
    # Filter out None results and exceptions, and update seen_uniprot_ids
    successful_results = []
    for result in results:
        if isinstance(result, dict):
            # Add UniProt ID to seen set to prevent duplicates in future batches
            seen_uniprot_ids.add(result['uniprot_id'])
            successful_results.append(result)
        elif isinstance(result, Exception):
            print(f"Error: PDB batch processing error: {result}")
    
    return successful_results

async def process_afdb_batch(processor: ProteinCoverageAnalyzer,
                            afdb_pairs: List[tuple], 
                            afdb_domains: Dict,
                            seen_uniprot_ids: Set[str]) -> List[Dict]:
    """Process a batch of AlphaFold proteins concurrently."""
    print(f"Processing AlphaFold batch of {len(afdb_pairs)} proteins...")
    
    tasks = []
    for uniprot_id, chain_id in afdb_pairs:
        task = process_single_afdb_protein(processor, uniprot_id, chain_id, afdb_domains, seen_uniprot_ids)
        tasks.append(task)
    
    results = await asyncio.gather(*tasks, return_exceptions=True)
    
    # Filter out None results and exceptions, and update seen_uniprot_ids
    successful_results = []
    for result in results:
        if isinstance(result, dict):
            # Add UniProt ID to seen set to prevent duplicates in future batches
            seen_uniprot_ids.add(result['uniprot_id'])
            successful_results.append(result)
        elif isinstance(result, Exception):
            print(f"Error: AlphaFold batch processing error: {result}")
    
    return successful_results

async def sample_pdb_until_target(target_count: int, 
                                 pdb_domains: Dict,
                                 processor: ProteinCoverageAnalyzer,
                                 seen_uniprot_ids: Set[str]) -> List[Dict]:
    """Sample PDB proteins until target is reached."""
    print(f"Starting PDB sampling to get {target_count} successful datapoints...")
    
    data = []
    processed_pairs = set()
    
    # Get all possible PDB chain pairs
    all_pdb_pairs = create_pdb_chain_list(pdb_domains)
    
    iteration = 0
    while len(data) < target_count:
        iteration += 1
        
        # Check remaining pairs
        remaining_pairs = [pair for pair in all_pdb_pairs if pair not in processed_pairs]
        if not remaining_pairs:
            print(f"Warning: No more PDB pairs to process! Got {len(data)} datapoints.")
            break
        
        # Calculate batch size
        needed = target_count - len(data)
        estimated_needed = min(needed * 3, len(remaining_pairs))
        actual_batch_size = min(BATCH_SIZE, estimated_needed)
        
        print(f"\n--- PDB Iteration {iteration} ---")
        print(f"Current successful datapoints: {len(data)}")
        print(f"Unique UniProt IDs seen: {len(seen_uniprot_ids)}")
        print(f"Target: {target_count}")
        print(f"Remaining pairs: {len(remaining_pairs)}")
        print(f"Sampling batch size: {actual_batch_size}")
        
        # Sample batch ensuring unique PDB IDs
        sampled_pairs = sample_unique_pdb_chains(
            {pdb_id: pdb_domains[pdb_id] for pdb_id, _ in remaining_pairs 
             if pdb_id not in [pair[0] for pair in processed_pairs]}, 
            actual_batch_size
        )
        
        # Process batch
        start_time = time.time()
        batch_results = await process_pdb_batch(processor, sampled_pairs, pdb_domains, seen_uniprot_ids)
        batch_time = time.time() - start_time
        
        # Update data
        data.extend(batch_results)
        processed_pairs.update(sampled_pairs)
        
        success_rate = len(batch_results) / len(sampled_pairs) if sampled_pairs else 0
        print(f"PDB Batch {iteration} completed in {batch_time:.1f}s: {len(batch_results)}/{len(sampled_pairs)} successful ({success_rate:.1%})")
        print(f"Total successful datapoints: {len(data)}")
        
        # Safety check
        if iteration > 50:
            print(f"Warning: Reached maximum iterations. Stopping.")
            break
    
    return data

async def sample_afdb_until_target(target_count: int, 
                                  afdb_domains: Dict,
                                  processor: ProteinCoverageAnalyzer,
                                  seen_uniprot_ids: Set[str]) -> List[Dict]:
    """Sample AlphaFold proteins until target is reached."""
    print(f"Starting AlphaFold sampling to get {target_count} successful datapoints...")
    
    data = []
    processed_pairs = set()
    
    # Get all possible AlphaFold pairs
    all_afdb_pairs = []
    for uniprot_id, chains in afdb_domains.items():
        for chain_id in chains.keys():
            all_afdb_pairs.append((uniprot_id, chain_id))
    
    random.seed(RANDOM_SEED)
    
    iteration = 0
    while len(data) < target_count:
        iteration += 1
        
        # Check remaining pairs
        remaining_pairs = [pair for pair in all_afdb_pairs if pair not in processed_pairs]
        if not remaining_pairs:
            print(f"Warning: No more AlphaFold pairs to process! Got {len(data)} datapoints.")
            break
        
        # Calculate batch size
        needed = target_count - len(data)
        estimated_needed = min(needed * 3, len(remaining_pairs))
        actual_batch_size = min(ALPHAFOLD_BATCH_SIZE, estimated_needed)
        
        print(f"\n--- AlphaFold Iteration {iteration} ---")
        print(f"Current successful datapoints: {len(data)}")
        print(f"Unique UniProt IDs seen: {len(seen_uniprot_ids)}")
        print(f"Target: {target_count}")
        print(f"Remaining pairs: {len(remaining_pairs)}")
        print(f"Sampling batch size: {actual_batch_size}")
        
        # Sample batch
        sampled_pairs = random.sample(remaining_pairs, actual_batch_size)
        
        # Process batch
        start_time = time.time()
        batch_results = await process_afdb_batch(processor, sampled_pairs, afdb_domains, seen_uniprot_ids)
        batch_time = time.time() - start_time
        
        # Update data
        data.extend(batch_results)
        processed_pairs.update(sampled_pairs)
        
        success_rate = len(batch_results) / len(sampled_pairs) if sampled_pairs else 0
        print(f"AlphaFold Batch {iteration} completed in {batch_time:.1f}s: {len(batch_results)}/{len(sampled_pairs)} successful ({success_rate:.1%})")
        print(f"Total successful datapoints: {len(data)}")
        
        # Safety check
        if iteration > 50:
            print(f"Warning: Reached maximum iterations. Stopping.")
            break
    
    return data

# Main execution
async def main():
    # Set random seed for reproducibility
    random.seed(RANDOM_SEED)
    
    parser = ECODParser("ecod-F40-domains.txt")
    ecod_domain_dict = order_domains_by_start(parser)

    def seperate_pdb_and_afdb(ecod_domain_dict):
        pdb_domains = {}
        afdb_domains = {}
        for id, domains in ecod_domain_dict.items():
            if len(id) == 4:  # PDB ID
                pdb_domains[id] = domains
            elif len(id) >= 6:  # AFDB ID
                afdb_domains[id] = domains
            else:
                print(f"Unknown ID format: {id}")
        return pdb_domains, afdb_domains

    pdb_domains, afdb_domains = seperate_pdb_and_afdb(ecod_domain_dict)
    
    print(f"Number of PDB domains: {len(pdb_domains)}")
    print(f"Number of AFDB domains: {len(afdb_domains)}")
    
    # Initialize shared UniProt ID tracking set
    seen_uniprot_ids = set()
    
    async with ProteinCoverageAnalyzer(max_concurrent=MAX_CONCURRENT_REQUESTS, 
                                       request_delay=REQUEST_DELAY) as processor:
        
        # Sample PDB proteins
        print("\n" + "="*50)
        print("SAMPLING PDB PROTEINS")
        print("="*50)
        start_time = time.time()
        pdb_data = await sample_pdb_until_target(TARGET_PDB_DATAPOINTS, pdb_domains, processor, seen_uniprot_ids)
        pdb_time = time.time() - start_time
        
        # Sample AlphaFold proteins (continuing with the same seen_uniprot_ids set)
        print("\n" + "="*50)
        print("SAMPLING ALPHAFOLD PROTEINS")
        print("="*50)
        start_time = time.time()
        afdb_data = await sample_afdb_until_target(TARGET_AFDB_DATAPOINTS, afdb_domains, processor, seen_uniprot_ids)
        afdb_time = time.time() - start_time

    os.makedirs('test_data', exist_ok=True)
    
    pdb_data = pdb_data[:TARGET_PDB_DATAPOINTS]
    # Save PDB results
    pdb_output_file = f'test_data/ecod_test_set_{TARGET_PDB_DATAPOINTS}_pdb.json'
    with open(pdb_output_file, 'w') as f:
        json.dump(pdb_data, f, indent=4)
    
    afdb_data = afdb_data[:TARGET_AFDB_DATAPOINTS]
    # Save AlphaFold results
    afdb_output_file = f'test_data/ecod_test_set_{TARGET_AFDB_DATAPOINTS}_afdb.json'
    with open(afdb_output_file, 'w') as f:
        json.dump(afdb_data, f, indent=4)
    
    print(f"\n" + "="*50)
    print("FINAL RESULTS")
    print("="*50)
    print(f"Saved {len(pdb_data)} PDB datapoints to {pdb_output_file}")
    print(f"Saved {len(afdb_data)} AlphaFold datapoints to {afdb_output_file}")
    print(f"Total unique UniProt IDs across both datasets: {len(seen_uniprot_ids)}")
    print(f"PDB processing time: {pdb_time/60:.1f} minutes")
    print(f"AlphaFold processing time: {afdb_time/60:.1f} minutes")
    print(f"Total processing time: {(pdb_time + afdb_time)/60:.1f} minutes")
    
    # Print some statistics
    if pdb_data:
        coverages = [d['coverage'] for d in pdb_data]
        pdb_lengths = [d['pdb_length'] for d in pdb_data]
        af_lengths = [d['alphafold_length'] for d in pdb_data]
        pdb_uniprot_ids = [d['uniprot_id'] for d in pdb_data]
        
        print(f"\nPDB Dataset statistics:")
        print(f"Coverage range: {min(coverages):.3f} - {max(coverages):.3f}")
        print(f"PDB length range: {min(pdb_lengths)} - {max(pdb_lengths)}")
        print(f"AlphaFold length range: {min(af_lengths)} - {max(af_lengths)}")
        print(f"Average coverage: {sum(coverages)/len(coverages):.3f}")
        #print(f"Unique UniProt IDs in PDB dataset: {len(set(pdb_uniprot_ids))}")
        #print(f"PDB UniProt ID uniqueness: {len(set(pdb_uniprot_ids)) == len(pdb_uniprot_ids)}")
    
    if afdb_data:
        afdb_lengths = [d['alphafold_length'] for d in afdb_data]
        afdb_uniprot_ids = [d['uniprot_id'] for d in afdb_data]
        
        print(f"\nAlphaFold Dataset statistics:")
        print(f"AlphaFold length range: {min(afdb_lengths)} - {max(afdb_lengths)}")
        print(f"Average AlphaFold length: {sum(afdb_lengths)/len(afdb_lengths):.1f}")
        #print(f"Unique UniProt IDs in AFDB dataset: {len(set(afdb_uniprot_ids))}")
        #print(f"AFDB UniProt ID uniqueness: {len(set(afdb_uniprot_ids)) == len(afdb_uniprot_ids)}")
    
    # Check for overlap between PDB and AFDB datasets
    if pdb_data and afdb_data:
        pdb_uniprot_set = set(d['uniprot_id'] for d in pdb_data)
        afdb_uniprot_set = set(d['uniprot_id'] for d in afdb_data)
        overlap = pdb_uniprot_set.intersection(afdb_uniprot_set)
        print(f"\nCross-dataset UniProt ID overlap: {len(overlap)} UniProt IDs")
        if overlap:
            print(f"Note: These UniProt IDs appear in both PDB and AFDB datasets: {list(overlap)[:10]}{'...' if len(overlap) > 10 else ''}")

if __name__ == "__main__":
    asyncio.run(main())