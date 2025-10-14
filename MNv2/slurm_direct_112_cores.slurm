#!/bin/bash
#SBATCH --job-name=collatz_112_cores  
#SBATCH --account=nct_352
#SBATCH --qos=gp_debug
#SBATCH --partition=gpp
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=112
#SBATCH --time=00:10:00
#SBATCH --output=/gpfs/scratch/nct_352/nct01225/collatz_output/cores_112_%j.out
#SBATCH --error=/gpfs/scratch/nct_352/nct01225/collatz_output/cores_112_%j.err

# Load Intel compiler environment
source /etc/profile.d/modules.sh 2>/dev/null || module purge
module purge
module load oneapi/2023.2

echo "==========================================="
echo "DIRECT 112 CORES REQUEST TEST"
echo "==========================================="
echo "SLURM allocation: 1 task with 112 CPUs per task"
echo ""

# Show resources
echo "=== SLURM Resources ==="
echo "SLURM_CPUS_ON_NODE: $SLURM_CPUS_ON_NODE"
echo "SLURM_CPUS_PER_TASK: $SLURM_CPUS_PER_TASK"
echo "SLURM_NTASKS: $SLURM_NTASKS"
echo ""

# Test with 112 threads in single MPI process
export OMP_NUM_THREADS=112
export OMP_PROC_BIND=close
export OMP_PLACES=cores

echo "=== Single Process, 112 Threads ==="
echo "OMP_NUM_THREADS: $OMP_NUM_THREADS"

# Run with single MPI process using all 112 cores
mpirun -np 1 ./collatz_mpi_gpp 0 50000000 single_112cores

echo ""
echo "==========================================="
echo "Direct 112 cores test complete: $(date)"
echo "==========================================="