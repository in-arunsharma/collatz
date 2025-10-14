#!/bin/bash
#SBATCH --job-name=mnv2_numa_test
#SBATCH --account=nct_352
#SBATCH --qos=gp_debug
#SBATCH --partition=gpp
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=28
#SBATCH --time=00:05:00
#SBATCH --output=/gpfs/scratch/nct_352/nct01225/collatz_output/numa_test_%j.out
#SBATCH --error=/gpfs/scratch/nct_352/nct01225/collatz_output/numa_test_%j.err

# Test with single socket (28 cores) to avoid NUMA effects

echo "========================================="
echo "NUMA Test - Single Socket (28 cores)"
echo "========================================="

module purge
module load intel/2023.2
module load impi/2021.10

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PROC_BIND=close
export OMP_PLACES=cores

# Bind to socket 0 only
export OMP_PLACES="{0:28:1}"

echo "Configuration:"
echo "  Socket:    0 (first 28 cores)"
echo "  Threads:   $OMP_NUM_THREADS per task"
echo ""

# Same workload as before
START_OFFSET=0
COUNT=100000000
RUN_TAG="numa_test_${SLURM_JOB_ID}"

srun -n 2 --cpus-per-task=28 --cpu-bind=cores numactl --cpunodebind=0 --membind=0 ./collatz_mpi_gpp $START_OFFSET $COUNT $RUN_TAG

echo ""
echo "========================================="
echo "Test complete: $(date)"
echo "========================================="