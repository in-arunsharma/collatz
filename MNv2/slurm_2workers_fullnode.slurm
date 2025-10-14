#!/bin/bash
#SBATCH --job-name=collatz_2workers_fullnode
#SBATCH --account=nct_352
#SBATCH --qos=gp_debug
#SBATCH --partition=gpp
#SBATCH --nodes=1
#SBATCH --ntasks=3
#SBATCH --exclusive
#SBATCH --time=00:10:00
#SBATCH --output=/gpfs/scratch/nct_352/nct01225/collatz_output/2workers_fullnode_%j.out
#SBATCH --error=/gpfs/scratch/nct_352/nct01225/collatz_output/2workers_fullnode_%j.err

# Load Intel compiler environment
module purge 2>/dev/null
module load oneapi/2023.2

echo "==========================================="
echo "2 WORKERS WITH EXCLUSIVE NODE ACCESS"
echo "==========================================="
echo "Strategy: Request exclusive node, let workers share all 112 cores"
echo ""

# Show resources
echo "=== SLURM Resources ==="
echo "SLURM_CPUS_ON_NODE: $SLURM_CPUS_ON_NODE"
echo "SLURM_JOB_NUM_NODES: $SLURM_JOB_NUM_NODES"
echo "SLURM_NTASKS: $SLURM_NTASKS"
echo ""

# Each worker uses 56 threads
export OMP_NUM_THREADS=56
export OMP_PROC_BIND=close
export OMP_PLACES=cores

echo "=== Configuration ==="
echo "Node allocation: EXCLUSIVE (all 112 cores)"
echo "MPI ranks: 3 (1 master + 2 workers)"
echo "OMP_NUM_THREADS: $OMP_NUM_THREADS per worker"
echo "Total thread capacity: 2 × 56 = 112 threads"
echo ""

echo "=== Running 100M seeds test ==="
mpirun -np 3 ./collatz_mpi_gpp 0 100000000 2workers_fullnode

echo ""
echo "==========================================="
echo "2 workers fullnode test complete: $(date)"
echo "==========================================="