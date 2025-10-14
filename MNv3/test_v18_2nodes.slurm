#!/bin/bash
#SBATCH --job-name=v18_2nodes
#SBATCH --output=v18_2nodes_%j.out
#SBATCH --error=v18_2nodes_%j.err
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=1       # 1 MPI rank per node
#SBATCH --cpus-per-task=112       # Full GPP node
#SBATCH --time=00:10:00
#SBATCH --qos=gp_debug
#SBATCH --account=nct_352

echo "=== V1.8-mpi-simple: 2 Nodes with Proper Binding ==="
echo "Date: $(date)"
echo "Nodes: $SLURM_NODELIST"
echo "Job: $SLURM_JOB_ID"
echo ""

module load intel impi

# Critical: OpenMP configuration
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PROC_BIND=spread
export OMP_PLACES=cores

echo "Configuration:"
echo "  MPI ranks: $SLURM_NTASKS"
echo "  CPUs/task: $SLURM_CPUS_PER_TASK"
echo "  OMP_NUM_THREADS: $OMP_NUM_THREADS"
echo "  OMP_PROC_BIND: $OMP_PROC_BIND"
echo "  OMP_PLACES: $OMP_PLACES"
echo ""

# Use srun with proper CPU binding (as per friend's advice)
srun --cpu-bind=cores ./V1.8-mpi-simple 0 100000000 --threads 112 --small-limit 20

echo ""
echo "Expected: ~274M nums/sec (2 × 137M)"
