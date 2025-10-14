#!/bin/bash
#SBATCH --job-name=collatz_2workers_asymmetric
#SBATCH --account=nct_352
#SBATCH --qos=gp_debug
#SBATCH --partition=gpp
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=3
#SBATCH --cpus-per-task=1
#SBATCH --time=00:10:00
#SBATCH --output=/gpfs/scratch/nct_352/nct01225/collatz_output/2workers_asymmetric_%j.out
#SBATCH --error=/gpfs/scratch/nct_352/nct01225/collatz_output/2workers_asymmetric_%j.err

# Load Intel compiler environment
module purge 2>/dev/null
module load oneapi/2023.2

echo "==========================================="
echo "2 WORKERS ASYMMETRIC CONFIGURATION"
echo "==========================================="
echo "Strategy: Allocate minimally, let OpenMP expand to all cores"
echo ""

# Show resources
echo "=== SLURM Resources ==="
echo "SLURM_CPUS_ON_NODE: $SLURM_CPUS_ON_NODE"
echo "SLURM_CPUS_PER_TASK: $SLURM_CPUS_PER_TASK"
echo "SLURM_NTASKS: $SLURM_NTASKS"
echo "SLURM_JOB_CPUS_PER_NODE: $SLURM_JOB_CPUS_PER_NODE"
echo ""

# Let OpenMP use all available cores
export OMP_NUM_THREADS=56
export OMP_PROC_BIND=close
export OMP_PLACES=cores

echo "=== Configuration ==="
echo "MPI ranks: 3 (1 master + 2 workers)"
echo "OMP_NUM_THREADS: $OMP_NUM_THREADS per worker"
echo "Strategy: Let workers expand to available cores"
echo ""

echo "=== Running 100M seeds test ==="
# Use Intel MPI process placement to distribute workers across sockets
export I_MPI_PIN_DOMAIN=auto
export I_MPI_PIN_PROCESSOR_EXCLUDE_LIST=0-55:56-111

mpirun -np 3 ./collatz_mpi_gpp 0 100000000 2workers_asymmetric

echo ""
echo "==========================================="
echo "2 workers asymmetric test complete: $(date)"
echo "==========================================="