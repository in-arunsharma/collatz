#!/bin/bash
#SBATCH --job-name=collatz-omp-array
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=80
#SBATCH --time=01:00:00
#SBATCH --array=0-99          # 100 nodes, adjust as needed
#SBATCH --output=collatz_omp_%A_%a.out
#SBATCH --error=collatz_omp_%A_%a.err

# MareNostrum 5 Multi-Node OpenMP Array Job
# Each array task processes a disjoint range
# Target: 100 nodes × 80 cores × 2.4M ≈ 19.2B nums/sec

echo "=== MareNostrum 5 OpenMP Array Job ==="
echo "Job ID: $SLURM_ARRAY_JOB_ID"
echo "Array Task ID: $SLURM_ARRAY_TASK_ID"
echo "Node: $SLURM_NODELIST"
echo "Cores: $SLURM_CPUS_PER_TASK"

# OpenMP configuration
export OMP_NUM_THREADS=80
export OMP_PROC_BIND=close
export OMP_PLACES=cores
export OMP_SCHEDULE=static

# Disjoint range per task (1B seeds per task)
RANGE_SIZE=1000000000
START_OFFSET=$((SLURM_ARRAY_TASK_ID * RANGE_SIZE))
COUNT=$RANGE_SIZE
TAG="mn5_job${SLURM_ARRAY_JOB_ID}_r${SLURM_ARRAY_TASK_ID}"

echo "Task: $SLURM_ARRAY_TASK_ID / $SLURM_ARRAY_TASK_MAX"
echo "Range: [2^71 + $START_OFFSET, 2^71 + $START_OFFSET + $COUNT)"
echo "Tag: $TAG"
echo ""

# Execute
./V1.4b-openmp $START_OFFSET $COUNT --tag $TAG

echo ""
echo "=== Task $SLURM_ARRAY_TASK_ID Complete ==="
echo "Check logs: *_${TAG}_t*.txt"
echo "Metadata: run_metadata_${TAG}.json"
