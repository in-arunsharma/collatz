#!/bin/bash
#SBATCH --job-name=collatz_full_cores
#SBATCH --account=nct_352
#SBATCH --qos=gp_debug
#SBATCH --partition=gpp
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=56
#SBATCH --time=00:10:00
#SBATCH --output=/gpfs/scratch/nct_352/nct01225/collatz_output/full_cores_%j.out
#SBATCH --error=/gpfs/scratch/nct_352/nct01225/collatz_output/full_cores_%j.err

# Load Intel compiler environment
source /etc/profile.d/modules.sh
module purge
module load oneapi/2023.2

echo "==========================================="
echo "FULL CORES UTILIZATION EXPERIMENT"
echo "==========================================="
echo "Testing different strategies to use all 112 cores per node"
echo ""

# Show available resources
echo "=== Node Resources ==="
echo "SLURM_CPUS_ON_NODE: $SLURM_CPUS_ON_NODE"
echo "SLURM_CPUS_PER_TASK: $SLURM_CPUS_PER_TASK" 
echo "SLURM_NTASKS: $SLURM_NTASKS"
echo "SLURM_NNODES: $SLURM_NNODES"
echo ""

# Show CPU info
echo "=== CPU Topology ==="
lscpu | grep -E "(CPU\(s\)|Thread|Core|Socket)"
echo ""

# Strategy 1: Try OMP_NUM_THREADS=112 (override SLURM setting)
echo "=== Strategy 1: Force OMP_NUM_THREADS=112 ==="
export OMP_NUM_THREADS=112
export OMP_PROC_BIND=close
export OMP_PLACES=cores

echo "OMP_NUM_THREADS: $OMP_NUM_THREADS"
echo "Testing if OpenMP can access all 112 cores..."

# Quick test with 10M seeds - Intel MPI syntax
mpirun -np 2 ./collatz_mpi_gpp 0 10000000 strategy1

echo ""
echo "=== Strategy 1 Results Above ==="
echo ""

# Strategy 2: Reset to safe 56 threads for comparison
echo "=== Strategy 2: Safe 56 threads (baseline) ==="
export OMP_NUM_THREADS=56
echo "OMP_NUM_THREADS: $OMP_NUM_THREADS"

# Same 10M test for comparison
mpirun -np 2 ./collatz_mpi_gpp 0 10000000 strategy2

echo ""
echo "=== Strategy 2 Results Above ==="
echo ""

echo "==========================================="
echo "Full cores experiment complete: $(date)"
echo "==========================================="