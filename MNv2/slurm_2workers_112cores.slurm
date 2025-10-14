#!/bin/bash
#SBATCH --job-name=collatz_2workers
#SBATCH --account=nct_352
#SBATCH --qos=gp_debug
#SBATCH --partition=gpp
#SBATCH --nodes=1
#SBATCH --ntasks=3
#SBATCH --cpus-per-task=56
#SBATCH --time=00:10:00
#SBATCH --output=/gpfs/scratch/nct_352/nct01225/collatz_output/2workers_%j.out
#SBATCH --error=/gpfs/scratch/nct_352/nct01225/collatz_output/2workers_%j.err

# Load Intel compiler environment
module purge 2>/dev/null
module load oneapi/2023.2

echo "==========================================="
echo "2 WORKERS × 56 THREADS = 112 CORES TOTAL"
echo "==========================================="
echo "Optimal configuration: Each worker on separate NUMA node"
echo ""

# Show resources
echo "=== SLURM Resources ==="
echo "SLURM_CPUS_ON_NODE: $SLURM_CPUS_ON_NODE"
echo "SLURM_CPUS_PER_TASK: $SLURM_CPUS_PER_TASK"
echo "SLURM_NTASKS: $SLURM_NTASKS (1 master + 2 workers)"
echo ""

# Configure OpenMP for 56 threads per worker
export OMP_NUM_THREADS=56
export OMP_PROC_BIND=close
export OMP_PLACES=cores

echo "=== Configuration ==="
echo "MPI ranks: 3 (1 master + 2 workers)"
echo "Threads per worker: 56"
echo "Total cores utilized: 112"
echo "Architecture: Master coordinates, workers compute on separate NUMA nodes"
echo ""

echo "=== Running 100M seeds test ==="
# Run with proper CPU binding for NUMA awareness
mpirun -np 3 ./collatz_mpi_gpp 0 100000000 2workers_test

echo ""
echo "==========================================="
echo "2 workers test complete: $(date)"
echo "==========================================="