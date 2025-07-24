#!/bin/bash
#SBATCH --job-name=protein_chunks
#SBATCH --array=1-16
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64GB
#SBATCH --time=01:00:00
#SBATCH --output=logs/chunk_%A_%a.out
#SBATCH --error=logs/chunk_%A_%a.error
#SBATCH --partition=doduo


export PATH="${VSC_SCRATCH}/miniforge3/bin:${PATH}"
source ${VSC_SCRATCH}/miniforge3/etc/profile.d/conda.sh
conda activate benchmarking

# Calculate protein range for this chunk
# 1000 proteins divided into 16 chunks = ~63 proteins per chunk
PROTEINS_PER_CHUNK=64
START_PROTEIN=$(( (SLURM_ARRAY_TASK_ID - 1) * PROTEINS_PER_CHUNK + 1 ))
END_PROTEIN=$(( SLURM_ARRAY_TASK_ID * PROTEINS_PER_CHUNK ))

# Handle the last chunk which might have fewer proteins
if [ $END_PROTEIN -gt 1000 ]; then
    END_PROTEIN=1000
fi

# Process this chunk of proteins
echo "Processing proteins $START_PROTEIN to $END_PROTEIN"
${VSC_SCRATCH}/miniforge3/envs/benchmarking/bin/python 03_parameter_testing_chunks.py $START_PROTEIN $END_PROTEIN

# Check if the chunk succeeded
if [ $? -ne 0 ]; then
    echo "Failed to process chunk $START_PROTEIN to $END_PROTEIN"
    exit 1
fi

echo "Completed chunk $SLURM_ARRAY_TASK_ID (proteins $START_PROTEIN-$END_PROTEIN)"