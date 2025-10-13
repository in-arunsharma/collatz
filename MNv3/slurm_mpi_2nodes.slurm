#!/bin/bash
#SBATCH --job-name=collatz_mpi2
#SBATCH --output=collatz_mpi2_%j.out
#SBATCH --error=collatz_mpi2_%j.err
#SBATCH --nodes=2
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=112
#SBATCH --time=00:30:00
#SBATCH --qos=gp_debug
#SBATCH --account=nct_352

# MareNostrum 5 GPP - MPI + OpenMP 2-Node Test
# 2 nodes × 112 cores = 224 cores total
# Expected: 274M nums/sec (based on 137M/node)

echo "========================================="
echo "Collatz V1.6-mpi-openmp - 2-Node Test"
echo "========================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Nodes: $SLURM_NNODES"
echo "Tasks: $SLURM_NTASKS"
echo "CPUs per task: $SLURM_CPUS_PER_TASK"
echo "Total cores: $((SLURM_NNODES * SLURM_CPUS_PER_TASK))"
echo "Date: $(date)"
echo ""

# Set OpenMP environment
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PROC_BIND=close
export OMP_PLACES=cores

# Test parameters - smaller for quick validation
START_OFFSET=0
COUNT=100000000  # 100M numbers (~0.4 seconds with 2 nodes)
MEMO_BITS=20
TAG="mn5_2nodes_${SLURM_JOB_ID}"

echo "Test Parameters:"
echo "  Start offset: $START_OFFSET"
echo "  Total count: $COUNT"
echo "  Per node: $((COUNT / SLURM_NNODES))"
echo "  Memo bits: $MEMO_BITS"
echo ""

echo "Starting computation..."
echo ""

time srun ./V1.6-mpi-openmp $START_OFFSET $COUNT \
    --threads $OMP_NUM_THREADS \
    --small-limit $MEMO_BITS \
    --tag $TAG

echo ""
echo "========================================="
echo "Job completed at $(date)"
echo "========================================="
