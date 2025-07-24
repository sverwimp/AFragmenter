import pickle
from glob import glob
import os
import pandas as pd

def load_pickle_files(checkpoint_dir):
    files = glob(os.path.join(checkpoint_dir, '*.pkl'))
    data = []
    for file in files:
        with open(file, 'rb') as f:
            data.append(pickle.load(f))
    return data


def prepare_data(data):
    """
    Extract and prepare data for visualization
    Assumes data is a list where each element is a dict with (threshold, resolution) keys
    """
    data_rows = []
    
    for protein_data in data:
        for (threshold, resolution), metrics in protein_data.items():
            data_rows.append({
                'threshold': threshold,
                'resolution': resolution,
                'iou': round(metrics['iou'], 3),
                'identifier': metrics['identifier'],
                'uniprot_id': metrics['uniprot_id']
            })
    
    return pd.DataFrame(data_rows)



if __name__ == "__main__":
    
    # Load data from CATH, ECOD PDB, and ECOD AFDB checkpoints
    cath = load_pickle_files('checkpoints_cath')
    ecod_pdb = load_pickle_files('checkpoints_ecod_pdb')
    ecod_afdb = load_pickle_files('checkpoints_ecod_afdb')

    print(f"CATH data loaded: {len(cath)} files")
    print(f"ECOD PDB data loaded: {len(ecod_pdb)} files")
    print(f"ECOD AFDB data loaded: {len(ecod_afdb)} files")

    df_cath = prepare_data(cath)
    df_ecod_pdb = prepare_data(ecod_pdb)
    df_ecod_afdb = prepare_data(ecod_afdb)

    os.makedirs('results_csv', exist_ok=True)

    df_cath.to_csv('results_csv/cath.csv', index=False)
    df_ecod_pdb.to_csv('results_csv/ecod_pdb.csv', index=False)
    df_ecod_afdb.to_csv('results_csv/ecod_afdb.csv', index=False)

