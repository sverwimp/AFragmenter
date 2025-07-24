import random
import json
import asyncio
import time
from typing import List, Dict, Optional, Set
from benchmarking_code import ProteinCoverageAnalyzer, CathDomallParser
import os

# Configuration
TARGET_DATAPOINTS = 1000
BATCH_SIZE = 200
MIN_COVERAGE = 0.9
MAX_COVERAGE = 1.1
MAX_INDEX_RATIO = 1.1
RANDOM_SEED = 42
MAX_CONCURRENT_REQUESTS = 20
REQUEST_DELAY = 0.05


def read_cath_dataset_unique_pdb_chain(file_path: str) -> list:
    """Read CATH dataset and return unique PDB+chain IDs, dropping domain info."""
    unique_pdb_chain_ids = set()
    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()
            if not line:
                continue
            # Extract pdb_id + chain (first 5 chars)
            pdb_chain = line[:5]  # e.g. '12asA' → '12asa'
            unique_pdb_chain_ids.add(pdb_chain)
    return list(unique_pdb_chain_ids)

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

async def process_single_protein(processor: ProteinCoverageAnalyzer, 
                                pdb_chain_id: str, 
                                cath_domain_dict: Dict,
                                seen_uniprot_ids: Set[str]) -> Optional[Dict]:  # Removed seen_pdb_ids parameter
    """Process a single protein asynchronously."""
    pdb_id = pdb_chain_id[:4]
    chain_id = pdb_chain_id[4]
    
    try:
        # Get UniProt ID
        uniprot_id = await processor.get_uniprot_id(pdb_id, chain_id)
        if not uniprot_id:
            return None
        
        # Check if UniProt ID is already in the dataset
        if uniprot_id in seen_uniprot_ids:
            #print(f"Skipping {pdb_id}:{chain_id} - UniProt ID {uniprot_id} already in dataset")
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
        true_chopping = cath_domain_dict.get(pdb_id, {}).get(chain_id, None)
        if not true_chopping:
            return None
        
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
        print(f"Error processing {pdb_id}:{chain_id}: {e}")
        return None

async def process_batch(processor: ProteinCoverageAnalyzer,
                       pdb_chain_ids: List[str], 
                       cath_domain_dict: Dict,
                       seen_uniprot_ids: Set[str],
                       seen_pdb_ids: Set[str]) -> List[Dict]: 
    """Process a batch of proteins concurrently."""
    print(f"Processing batch of {len(pdb_chain_ids)} proteins...")
    
    tasks = []
    for pdb_chain_id in pdb_chain_ids:
        task = process_single_protein(processor, pdb_chain_id, cath_domain_dict, seen_uniprot_ids)  # Removed seen_pdb_ids
        tasks.append(task)
    
    results = await asyncio.gather(*tasks, return_exceptions=True)
    
    # Filter out None results and exceptions, and update seen sets
    successful_results = []
    for result in results:
        if isinstance(result, dict):
            seen_uniprot_ids.add(result['uniprot_id'])
            seen_pdb_ids.add(result['pdb_id'])
            successful_results.append(result)
        elif isinstance(result, Exception):
            print(f"Error: Batch processing error: {result}")
    
    return successful_results

async def sample_until_target_async(target_count: int, 
                                   unique_pdb_chain_ids: List[str],
                                   cath_domain_dict: Dict) -> List[Dict]:
    """Sample proteins until target is reached using async processing."""
    
    random.seed(RANDOM_SEED)
    data = []
    processed_ids = set()
    seen_uniprot_ids = set() 
    seen_pdb_ids = set() 
    
    print(f"Starting async sampling to get {target_count} successful datapoints...")
    print(f"Available unique PDB chain IDs: {len(unique_pdb_chain_ids)}")
    
    async with ProteinCoverageAnalyzer(max_concurrent=MAX_CONCURRENT_REQUESTS, 
                                       request_delay=REQUEST_DELAY) as processor:
        iteration = 0
        
        while len(data) < target_count:
            iteration += 1
            
            # Check remaining proteins
            remaining_proteins = len(unique_pdb_chain_ids) - len(processed_ids)
            if remaining_proteins == 0:
                print(f"Warning: No more proteins to process! Got {len(data)} datapoints.")
                break
            
            # Calculate batch size
            needed = target_count - len(data)
            estimated_needed = min(needed * 3, remaining_proteins)
            actual_batch_size = min(BATCH_SIZE, estimated_needed, remaining_proteins)
            
            print(f"\n--- Iteration {iteration} ---")
            print(f"Current successful datapoints: {len(data)}")
            print(f"Unique UniProt IDs seen: {len(seen_uniprot_ids)}")
            print(f"Unique PDB IDs seen: {len(seen_pdb_ids)}")
            print(f"Target: {target_count}")
            print(f"Remaining proteins: {remaining_proteins}")
            print(f"Sampling batch size: {actual_batch_size}")
            
            if actual_batch_size == 0:
                print("Warning: No more valid proteins to process (all remaining have duplicate PDB IDs)")
                break
            
            # Sample batch ensuring no duplicate PDB IDs within the batch
            available_proteins = [pid for pid in unique_pdb_chain_ids if pid not in processed_ids]
            batch_ids = []
            batch_pdb_ids = set()
            
            # Shuffle available proteins for random sampling
            random.shuffle(available_proteins)
            
            # Select proteins ensuring no duplicate PDB IDs within this batch
            for pdb_chain_id in available_proteins:
                pdb_id = pdb_chain_id[:4]
                
                # Skip if this PDB ID is already in the current batch or already seen
                if pdb_id in batch_pdb_ids or pdb_id in seen_pdb_ids:
                    continue
                
                batch_ids.append(pdb_chain_id)
                batch_pdb_ids.add(pdb_id)
                
                # Stop when we have enough for this batch
                if len(batch_ids) >= actual_batch_size:
                    break
            
            # Update actual batch size based on what we could actually select
            actual_batch_size = len(batch_ids)
            
            # Process batch
            start_time = time.time()
            batch_results = await process_batch(processor, batch_ids, cath_domain_dict, seen_uniprot_ids, seen_pdb_ids)
            batch_time = time.time() - start_time
            
            # Update data
            data.extend(batch_results)
            processed_ids.update(batch_ids)
            
            success_rate = len(batch_results) / len(batch_ids) if batch_ids else 0
            duplicate_skips = len(batch_ids) - len(batch_results) - sum(1 for r in batch_results if r is None)
            
            print(f"Batch {iteration} completed in {batch_time:.1f}s: {len(batch_results)}/{len(batch_ids)} successful ({success_rate:.1%})")
            print(f"Skipped {duplicate_skips} duplicates or invalid entries in this batch.")
            print(f"Total successful datapoints: {len(data)}")
            
            # Safety check
            if iteration > 50:
                print(f"Warning: Reached maximum iterations. Stopping.")
                break
    
    print(f"\nAsync sampling complete!")
    print(f"Final successful datapoints: {len(data)}")
    print(f"Unique UniProt IDs in dataset: {len(seen_uniprot_ids)}")
    print(f"Unique PDB IDs in dataset: {len(seen_pdb_ids)}")
    print(f"Total proteins processed: {len(processed_ids)}")
    print(f"Overall success rate: {len(data)/len(processed_ids):.1%}")
    
    return data

# Main execution
async def main():
    # Set random seed for reproducibility
    random.seed(RANDOM_SEED)
    
    # Load CATH data
    with open('cath-domain-boundaries.txt', 'r') as f:
        content = f.read()
    parser = CathDomallParser()
    cath_domain_dict = parser.parse_string(content)
    
    # Read CATH dataset
    print("Reading CATH dataset...")
    unique_pdb_chain_ids = read_cath_dataset_unique_pdb_chain('cath-dataset-nonredundant-S40.list')
    print(f"Found {len(unique_pdb_chain_ids)} unique PDB chain IDs")
    
    # Sample until we get target number of successful datapoints
    start_time = time.time()
    data = await sample_until_target_async(TARGET_DATAPOINTS, unique_pdb_chain_ids, cath_domain_dict)
    total_time = time.time() - start_time
    data = data[:TARGET_DATAPOINTS] # Ensure we only keep the target number of datapoints
    
    os.makedirs('test_data', exist_ok=True)
    # Save results
    output_file = f'test_data/cath_test_set.json'
    with open(output_file, 'w') as f:
        json.dump(data, f, indent=4)
    
    print(f"\nSaved {len(data)} successful datapoints to {output_file}")
    print(f"Total processing time: {total_time/60:.1f} minutes")
    
    # Print some statistics
    if data:
        coverages = [d['coverage'] for d in data]
        pdb_lengths = [d['pdb_length'] for d in data]
        af_lengths = [d['alphafold_length'] for d in data]
        
        print(f"\nDataset statistics:")
        print(f"Coverage range: {min(coverages):.3f} - {max(coverages):.3f}")
        print(f"PDB length range: {min(pdb_lengths)} - {max(pdb_lengths)}")
        print(f"AlphaFold length range: {min(af_lengths)} - {max(af_lengths)}")
        print(f"Average coverage: {sum(coverages)/len(coverages):.3f}")
        
        # Verify uniqueness
        unique_pdb_ids = set(d['pdb_id'] for d in data)
        unique_uniprot_ids = set(d['uniprot_id'] for d in data)
        print(f"\nUniqueness verification:")
        print(f"Total datapoints: {len(data)}")
        print(f"Unique PDB IDs: {len(unique_pdb_ids)}")
        print(f"Unique UniProt IDs: {len(unique_uniprot_ids)}")
        
        if len(unique_pdb_ids) != len(data):
            print(f"⚠️  WARNING: Found duplicate PDB IDs!")
        if len(unique_uniprot_ids) != len(data):
            print(f"⚠️  WARNING: Found duplicate UniProt IDs!")

if __name__ == "__main__":
    asyncio.run(main())