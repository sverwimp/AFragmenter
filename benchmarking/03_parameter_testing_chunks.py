#!/usr/bin/env python3
"""
Chunked version of your original code for multiple datasets
Usage: python 03_paramter_testing_chunks.py <start_index> <end_index>
"""

import sys
import multiprocessing
import pickle
import os
import signal
import json
import re

from afragmenter import AFragmenter, fetch_afdb_data
from benchmarking_code import iou_chain

timeout_sec = 2
objective_function = "modularity"

class TimeoutError(Exception):
    pass

def timeout_handler(signum, frame):
    raise TimeoutError("Clustering timed out")

def run_cluster_with_timeout(af, resolution, timeout_sec=2):
    try:
        signal.signal(signal.SIGALRM, timeout_handler)
        signal.alarm(timeout_sec)
        af.cluster(objective_function=objective_function, resolution=resolution)
        signal.alarm(0)
        return True
    except TimeoutError:
        af.cluster(objective_function=objective_function, resolution=resolution, n_iterations=20)
        print(f"-----Timed out result: {af.cluster_intervals}-----")
        return True
    except Exception as e:
        af.cluster_intervals = af.cluster_intervals = {'0': [(0, len(af.pae_matrix)-1)]}  # Full protein as one domain  #{'0': [(0, 0)]}
        return False
    finally:
        signal.alarm(0)

def convert_to_1based_indexing(cluster_intervals):
    """Convert 0-based cluster intervals to 1-based indexing"""
    converted = {}
    for domain_id, intervals in cluster_intervals.items():
        converted[domain_id] = [(start + 1, end + 1) for start, end in intervals]
    return converted

def cluster_intervals_to_domain_string(af):
    # Convert to 1-based indexing first
    converted_intervals = convert_to_1based_indexing(af.cluster_intervals)

    domain_strings = []
    for domain_id, domain in converted_intervals.items():
        if not domain:  # Skip empty domains
            continue
        single_string = "_".join([f"{d[0]}-{d[1]}" for d in domain if len(d) == 2])
        if single_string:  # Only add non-empty strings
            domain_strings.append(single_string)

    result = ",".join(domain_strings)

    # Debug: Check for problematic patterns
    if not result or result.startswith(',') or result.endswith(',') or ',,' in result:
        print(f"Warning: Malformed domain string: '{result}'")
        print(f"Original cluster_intervals: {af.cluster_intervals}")
        print(f"Converted cluster_intervals: {converted_intervals}")

    return result

def fix_negative_chopping(true_chopping):
    """Fix negative starting indices in true_chopping strings"""
    if true_chopping.startswith('-'):
        # Replace the first negative number with 1 (for 1-based indexing)
        fixed = re.sub(r'^-\d+', '1', true_chopping)
        return fixed
    return true_chopping

def create_checkpoint_dir(checkpoint_dir):
    if not os.path.exists(checkpoint_dir):
        os.makedirs(checkpoint_dir, exist_ok=True)
    return checkpoint_dir

def determine_dataset_info(item, cath_data=None, pdb_data=None):
    """Determine which dataset an item belongs to and return appropriate checkpoint directory and identifier"""
    # First check if we have the dataset source marked (preferred method)
    dataset_source = item.get('_dataset_source')

    if dataset_source == 'cath':
        return 'checkpoints_cath', item['pdb_id']
    elif dataset_source == 'pdb':
        return 'checkpoints_ecod_pdb', item['pdb_id']
    elif dataset_source == 'afdb':
        return 'checkpoints_ecod_afdb', item['uniprot_id']

    # Fallback to original logic if dataset source is not marked
    # Check if it's AFDB data (no pdb_id or pdb_length)
    if 'pdb_id' not in item or 'pdb_length' not in item:
        return 'checkpoints_ecod_afdb', item['uniprot_id']

    # For CATH and ECOD-PDB data (both have pdb_id and pdb_length), we need to check which dataset it belongs to
    pdb_id = item['pdb_id']

    # Check if this pdb_id exists in CATH data
    if cath_data is not None:
        cath_pdb_ids = {entry['pdb_id'] for entry in cath_data}
        if pdb_id in cath_pdb_ids:
            return 'checkpoints_cath', pdb_id

    # Check if this pdb_id exists in PDB data
    if pdb_data is not None:
        pdb_pdb_ids = {entry['pdb_id'] for entry in pdb_data}
        if pdb_id in pdb_pdb_ids:
            return 'checkpoints_ecod_pdb', pdb_id

    # If we can't determine from the provided data, default to PDB
    return 'checkpoints_ecod_pdb', pdb_id

def get_checkpoint_filename(identifier, checkpoint_dir):
    return os.path.join(checkpoint_dir, f"{identifier}.pkl")

def save_protein_checkpoint(identifier, protein_results, checkpoint_dir):
    checkpoint_file = get_checkpoint_filename(identifier, checkpoint_dir)
    with open(checkpoint_file, "wb") as f:
        pickle.dump(protein_results, f)

def load_protein_checkpoint(identifier, checkpoint_dir):
    checkpoint_file = get_checkpoint_filename(identifier, checkpoint_dir)
    if os.path.exists(checkpoint_file):
        with open(checkpoint_file, "rb") as f:
            return pickle.load(f)
    return None

def get_completed_proteins(checkpoint_dir):
    if not os.path.exists(checkpoint_dir):
        return set()

    completed = set()
    for filename in os.listdir(checkpoint_dir):
        if filename.endswith(".pkl"):
            identifier = filename.replace(".pkl", "")
            completed.add(identifier)
    return completed

def evaluate_protein_all_params(args):
    item, thresholds, resolutions, cath_data, pdb_data = args
    uniprot_id = item["uniprot_id"]
    true_chopping = fix_negative_chopping(item["true_chopping"])

    # Determine checkpoint directory and identifier based on dataset
    checkpoint_dir, identifier = determine_dataset_info(item, cath_data, pdb_data)
    create_checkpoint_dir(checkpoint_dir)

    # Check if this protein has already been processed
    existing_results = load_protein_checkpoint(identifier, checkpoint_dir)
    if existing_results is not None:
        print(f"Loading cached results for {identifier} from {checkpoint_dir}")
        return existing_results

    # Single fetch per protein
    try:
        json_data, pdb = fetch_afdb_data(uniprot_id)
    except Exception as e:
        print(f"Skipping {uniprot_id} (identifier: {identifier}) due to error: {e}")
        return None

    # Test all parameter combinations on this protein
    protein_results = {}

    for threshold in thresholds:
        for resolution in resolutions:
            try:
                # Create new AFragmenter instance for each parameter combination
                af = AFragmenter(json_data, threshold=threshold)
                success = run_cluster_with_timeout(af, resolution, timeout_sec=timeout_sec)

                print(f"Clustering {uniprot_id} (identifier: {identifier}) with threshold={threshold}, resolution={resolution} (success={success})")

                if not success or af.cluster_intervals == {} or len(af.cluster_intervals) == 0:
                    # Handle empty clusters case
                    af.cluster_intervals = {'0': [(0, len(af.pae_matrix)-1)]}  # Full protein as one domain

                prediction = cluster_intervals_to_domain_string(af)
                print(f"Processing {uniprot_id} (identifier: {identifier}) with threshold={threshold}, resolution={resolution}, k={prediction}")
                if prediction == '':
                    print(f"Converted cluster intervals to domain string for {uniprot_id}: {af.cluster_intervals}")

                # Handle iou_chain return value
                iou_result = iou_chain(true_chopping, prediction)
                if isinstance(iou_result, tuple):
                    iou, iou_per_domain = iou_result
                else:
                    iou = iou_result
                    iou_per_domain = [iou]

                print(f"IOU for {uniprot_id} (identifier: {identifier}): {iou}")
                print(f"    with true chopping: {true_chopping}")
                print(f"    with predicted chopping: {prediction}")

                protein_results[(threshold, resolution)] = {
                    "uniprot_id": uniprot_id,
                    "identifier": identifier,
                    "threshold": threshold,
                    "resolution": resolution,
                    "iou": iou,
                    "iou_per_domain": iou_per_domain,
                    "prediction": prediction,
                    "true_chopping": true_chopping,
                    "success": success,
                    "num_domains_predicted": len(af.cluster_intervals) if hasattr(af, 'cluster_intervals') and af.cluster_intervals else 0
                }
                print()

            except Exception as e:
                print(f"Error processing {uniprot_id} (identifier: {identifier}) with threshold={threshold}, resolution={resolution}: {e}")
                # Store a default result for this parameter combination
                protein_results[(threshold, resolution)] = {
                    "uniprot_id": uniprot_id,
                    "identifier": identifier,
                    "threshold": threshold,
                    "resolution": resolution,
                    "iou": 0.0,
                    "iou_per_domain": [0.0],
                    "prediction": "1-1",
                    "true_chopping": true_chopping,
                    "success": False,
                    "num_domains_predicted": 0
                }

    # Save checkpoint for this protein using the appropriate identifier
    save_protein_checkpoint(identifier, protein_results, checkpoint_dir)
    print(f"Finished processing {uniprot_id} (identifier: {identifier}) with {len(protein_results)} parameter combinations (saved to {checkpoint_dir})")
    return protein_results

def get_item_identifier(item, cath_data=None, pdb_data=None):
    """Get the identifier for an item based on its dataset type"""
    checkpoint_dir, identifier = determine_dataset_info(item, cath_data, pdb_data)
    return identifier

def main():
    if len(sys.argv) != 3:
        print("Usage: python chunked_protein_analysis.py <start_index> <end_index>")
        sys.exit(1)

    start_index = int(sys.argv[1]) - 1  # Convert to 0-based indexing
    end_index = int(sys.argv[2]) - 1

    # Load all datasets and track their source
    datasets = []
    cath_data = []
    pdb_data = []
    afdb_data = []

    # Load CATH data
    try:
        with open('test_data/cath_test_set.json', 'r') as f:
            cath_data = json.load(f)
        # Add dataset source to each item
        for item in cath_data:
            item['_dataset_source'] = 'cath'
        datasets.extend(cath_data)
        print(f"Loaded {len(cath_data)} CATH proteins")
    except FileNotFoundError:
        print("Warning: cath_test_set.json not found")

    try:
        with open('test_data/ecod_test_set_500_pdb.json', 'r') as f:
            pdb_data = json.load(f)[:500]
        # Add dataset source to each item
        for item in pdb_data:
            item['_dataset_source'] = 'pdb'
        datasets.extend(pdb_data)
        print(f"Loaded {len(pdb_data)} PDB proteins")
    except FileNotFoundError:
        print("Warning: ecod_test_set_500_pdb.json not found")

    try:
        with open('test_data/ecod_test_set_500_afdb.json', 'r') as f:
            afdb_data = json.load(f)[:500]
        # Add dataset source to each item
        for item in afdb_data:
            item['_dataset_source'] = 'afdb'
        datasets.extend(afdb_data)
        print(f"Loaded {len(afdb_data)} AFDB proteins")
    except FileNotFoundError:
        print("Warning: ecod_test_set_500_afdb.json not found")

    # Get subset of data for this chunk
    data_chunk = datasets[start_index:end_index + 1]
    print(f"Processing proteins {start_index + 1} to {end_index + 1} ({len(data_chunk)} proteins)")

    # Create checkpoint directories
    create_checkpoint_dir("checkpoints_cath")
    create_checkpoint_dir("checkpoints_ecod_pdb")
    create_checkpoint_dir("checkpoints_ecod_afdb")

    # Check which proteins have already been processed PER DATASET
    cath_completed = get_completed_proteins("checkpoints_cath")
    pdb_completed = get_completed_proteins("checkpoints_ecod_pdb")
    afdb_completed = get_completed_proteins("checkpoints_ecod_afdb")

    print(f"Found {len(cath_completed)} completed CATH proteins")
    print(f"Found {len(pdb_completed)} completed PDB proteins")
    print(f"Found {len(afdb_completed)} completed AFDB proteins")

    # Filter out already completed proteins based on their dataset-specific completion status
    df_remaining = []
    for item in data_chunk:
        # Use the dataset source we added to determine checkpoint directory and identifier
        dataset_source = item.get('_dataset_source', 'unknown')

        if dataset_source == 'cath':
            checkpoint_dir = 'checkpoints_cath'
            identifier = item['pdb_id']
            completed_set = cath_completed
        elif dataset_source == 'pdb':
            checkpoint_dir = 'checkpoints_ecod_pdb'
            identifier = item['pdb_id']
            completed_set = pdb_completed
        elif dataset_source == 'afdb':
            checkpoint_dir = 'checkpoints_ecod_afdb'
            identifier = item['uniprot_id']
            completed_set = afdb_completed
        else:
            print(f"Warning: Unknown dataset source for item: {item}")
            continue

        # Check if this specific protein in this specific dataset is completed
        if identifier not in completed_set:
            df_remaining.append(item)

    print(f"Remaining proteins to process: {len(df_remaining)}")

    thresholds = list(range(0, 11)) # 0 to 10 inclusive
    resolutions = [round(x * 0.1, 1) for x in range(1, 16)] # 0.1 to 1.5 inclusive with step 0.1

    if len(df_remaining) > 0:
        print(f"Processing {len(df_remaining)} remaining proteins...")

        # Prepare arguments for multiprocessing
        protein_args = []
        for item in df_remaining:
            protein_args.append((item, thresholds, resolutions, cath_data, pdb_data))

        # Use multiprocessing within this chunk
        num_processes = min(multiprocessing.cpu_count(), len(df_remaining), 4)
        print(f"Using {num_processes} processes for this chunk")

        with multiprocessing.Pool(processes=num_processes) as pool:
            results = []
            for result in pool.imap_unordered(evaluate_protein_all_params, protein_args):
                if result is not None:
                    results.append(result)
                    print(f"Completed protein, total done: {len(results)}/{len(df_remaining)}")

    print(f"Chunk {start_index + 1}-{end_index + 1} completed!")

if __name__ == '__main__':
    main()
