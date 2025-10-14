#!/bin/bash
#SBATCH --job-name=mnv2_cache_test
#SBATCH --account=nct_352
#SBATCH --qos=gp_debug  
#SBATCH --partition=gpp
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=56
#SBATCH --time=00:10:00
#SBATCH --output=/gpfs/scratch/nct_352/nct01225/collatz_output/cache_test_%j.out
#SBATCH --error=/gpfs/scratch/nct_352/nct01225/collatz_output/cache_test_%j.err

# Cache optimization test: Compare 4MB vs 1MB memo table

echo "========================================="
echo "Cache Optimization Test"
echo "========================================="

module purge
module load intel/2023.2
module load impi/2021.10

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PROC_BIND=close
export OMP_PLACES=cores

echo "Configuration:"
echo "  New memo table: 2^18 = 256K entries = 1MB (L2 cache optimized)"
echo "  Old memo table: 2^20 = 1M entries = 4MB (L3 cache, shared)"
echo "  Threads: $OMP_NUM_THREADS per task"
echo ""

# Same workload as successful test
START_OFFSET=0
COUNT=100000000
RUN_TAG="cache_optimized_${SLURM_JOB_ID}"

echo "Workload: $COUNT seeds"
echo ""

# Test with cache-optimized memo table
echo "=== Cache Optimized Version (1MB memo table) ==="
srun -n 2 --cpus-per-task=56 --cpu-bind=cores ./collatz_mpi_gpp $START_OFFSET $COUNT $RUN_TAG

echo ""
echo "Expected improvement: 20-40% throughput increase due to L2 cache efficiency"
echo ""
echo "========================================="
echo "Cache test complete: $(date)"
echo "========================================="