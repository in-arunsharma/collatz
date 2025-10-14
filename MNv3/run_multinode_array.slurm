#!/bin/bash
#SBATCH --job-name=v15_array
#SBATCH --output=v15_array_%A_%a.out
#SBATCH --error=v15_array_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=112
#SBATCH --time=08:00:00
#SBATCH --partition=gpp
#SBATCH --qos=gp_bsccs
#SBATCH --array=0-9          # 10 nodes (array indices 0-9)

# Load modules
module purge
module load gcc/13.2.0

# OpenMP configuration
export OMP_NUM_THREADS=112
export OMP_PLACES=cores
export OMP_PROC_BIND=spread

# Work distribution
# Total range: 2^71 to 2^71 + 10 billion
# Each array task gets 1 billion numbers (at 137M/sec = ~7 seconds per billion)
# For 8-hour run, each node can process ~3,936 billion numbers
# So we'll give each node a much larger chunk

TOTAL_RANGE=100000000000    # 100 billion total
NUM_NODES=10
CHUNK_SIZE=$((TOTAL_RANGE / NUM_NODES))  # 10 billion per node

START_OFFSET=$((SLURM_ARRAY_TASK_ID * CHUNK_SIZE))
COUNT=$CHUNK_SIZE

# For last node, take any remainder
if [ $SLURM_ARRAY_TASK_ID -eq $((NUM_NODES - 1)) ]; then
    COUNT=$((TOTAL_RANGE - START_OFFSET))
fi

echo "=== V1.5-openmp: Multi-Node Array Job ==="
echo "Array Task ID: $SLURM_ARRAY_TASK_ID"
echo "Node: $SLURM_NODELIST"
echo "Date: $(date)"
echo "Range: [$START_OFFSET, +$COUNT)"
echo ""

# Run with unique tag for this array task
./V1.5-openmp $START_OFFSET $COUNT \
    --threads 112 \
    --tag v15_node${SLURM_ARRAY_TASK_ID}_job${SLURM_ARRAY_JOB_ID}

echo ""
echo "=== COMPLETED ==="
echo "Array Task: $SLURM_ARRAY_TASK_ID"
echo "Exit code: $?"
