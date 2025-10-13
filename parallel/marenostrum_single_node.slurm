#!/bin/bash
#SBATCH --job-name=collatz-omp
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=80
#SBATCH --time=01:00:00
#SBATCH --output=collatz_omp_%j.out
#SBATCH --error=collatz_omp_%j.err

# MareNostrum 5 Single-Node OpenMP Job
# Target: 80 cores × 2.4M ≈ 192M nums/sec

echo "=== MareNostrum 5 OpenMP Collatz Job ==="
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo "Cores: $SLURM_CPUS_PER_TASK"

# OpenMP configuration for NUMA-aware execution
export OMP_NUM_THREADS=80
export OMP_PROC_BIND=close     # Bind threads to cores
export OMP_PLACES=cores        # One thread per core
export OMP_SCHEDULE=static     # Static scheduling (deterministic)

echo "OMP_NUM_THREADS=$OMP_NUM_THREADS"
echo "OMP_PROC_BIND=$OMP_PROC_BIND"
echo "OMP_PLACES=$OMP_PLACES"

# Run parameters
START_OFFSET=0
COUNT=1000000000  # 1 billion seeds
TAG="mn5_job${SLURM_JOB_ID}"

echo "Range: [2^71 + $START_OFFSET, 2^71 + $START_OFFSET + $COUNT)"
echo "Tag: $TAG"
echo ""

# Execute
./V1.4b-openmp $START_OFFSET $COUNT --tag $TAG

echo ""
echo "=== Job Complete ==="
echo "Check logs: *_${TAG}_t*.txt"
echo "Metadata: run_metadata_${TAG}.json"
