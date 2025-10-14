#!/bin/bash
#SBATCH --job-name=mnv2_optimized_56t
#SBATCH --account=nct_352
#SBATCH --qos=gp_debug
#SBATCH --partition=gpp
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=56
#SBATCH --time=00:10:00
#SBATCH --output=/gpfs/scratch/nct_352/nct01225/collatz_output/optimized_56t_%j.out
#SBATCH --error=/gpfs/scratch/nct_352/nct01225/collatz_output/optimized_56t_%j.err

# OPTIMIZED CODE TEST: Same 56 threads, but with all performance fixes applied

echo "========================================="
echo "OPTIMIZED CODE TEST (56 threads)"
echo "========================================="

module purge
module load intel/2023.2
module load impi/2021.10

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PROC_BIND=spread  
export OMP_PLACES=cores

echo "Configuration:"
echo "  OPTIMIZATIONS APPLIED:"
echo "  ✅ Peak tracking: Only on 3n+1 spikes (not every iteration)" 
echo "  ✅ Schedule: static (no chunk overhead)"
echo "  ✅ Critical section: Minimal (no vector operations)"
echo "  ✅ NUMA: Memory interleaving"
echo "  Threads: $OMP_NUM_THREADS per worker"
echo ""

# Same workload for direct comparison
START_OFFSET=0
COUNT=100000000
RUN_TAG="optimized_56t_${SLURM_JOB_ID}"

echo "Workload: $COUNT seeds"
echo ""

echo "=== OPTIMIZED vs BASELINE COMPARISON ==="

srun -n 2 --cpus-per-task=56 --cpu-bind=cores \
     numactl --interleave=all \
     ./collatz_mpi_gpp $START_OFFSET $COUNT $RUN_TAG

echo ""
echo "Performance Comparison:"
echo "  Baseline (56 threads, unoptimized): 68.33 M/sec"  
echo "  Target (56 threads, optimized):     85-95 M/sec"
echo "  Expected improvement:               +25-40%"
echo ""
echo "Optimizations that should improve performance:"
echo "  1. Peak tracking fix: +10-20% (major hot loop optimization)"
echo "  2. Schedule optimization: +5-10% (reduced overhead)"  
echo "  3. Critical section removal: +5-15% (less contention)"
echo "  4. NUMA optimization: +5-10% (better memory access)"
echo ""
echo "========================================="
echo "Optimized test complete: $(date)" 
echo "========================================="