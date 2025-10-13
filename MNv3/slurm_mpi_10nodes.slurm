#!/bin/bash
#SBATCH --job-name=collatz_mpi
#SBATCH --output=collatz_mpi_%j.out
#SBATCH --error=collatz_mpi_%j.err
#SBATCH --nodes=10
#SBATCH --ntasks=10
#SBATCH --cpus-per-task=112
#SBATCH --time=02:00:00
#SBATCH --qos=gp_bsccs
#SBATCH --account=nct_352

# MareNostrum 5 GPP - MPI + OpenMP Multi-Node
# 10 nodes × 112 cores = 1,120 cores total
# Expected: 1.37 BILLION nums/sec (based on 137M/node)

echo "========================================="
echo "Collatz V1.6-mpi-openmp - Multi-Node Test"
echo "========================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Nodes: $SLURM_NNODES"
echo "Tasks: $SLURM_NTASKS"
echo "CPUs per task: $SLURM_CPUS_PER_TASK"
echo "Total cores: $((SLURM_NNODES * SLURM_CPUS_PER_TASK))"
echo "Nodelist: $SLURM_NODELIST"
echo "Date: $(date)"
echo ""

# Load modules
# module load intel/2023.2.0
# module load impi/2021.9.0

# Set OpenMP environment
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PROC_BIND=close
export OMP_PLACES=cores

echo "MPI Configuration:"
echo "  MPI tasks: $SLURM_NTASKS"
echo "  OMP_NUM_THREADS=$OMP_NUM_THREADS"
echo "  OMP_PROC_BIND=$OMP_PROC_BIND"
echo ""

# Test parameters
START_OFFSET=0
COUNT=1000000000  # 1 BILLION numbers (~7 seconds with 10 nodes)
MEMO_BITS=20
TAG="mn5_10nodes_${SLURM_JOB_ID}"

echo "Test Parameters:"
echo "  Start offset: $START_OFFSET"
echo "  Total count: $COUNT"
echo "  Per node: $((COUNT / SLURM_NNODES))"
echo "  Memo bits: $MEMO_BITS"
echo "  Tag: $TAG"
echo ""

# Run with srun (SLURM's MPI launcher)
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
