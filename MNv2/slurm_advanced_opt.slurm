#!/bin/bash
#SBATCH --job-name=mnv2_advanced_opt
#SBATCH --account=nct_352
#SBATCH --qos=gp_debug
#SBATCH --partition=gpp
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=56
#SBATCH --time=00:10:00
#SBATCH --output=/gpfs/scratch/nct_352/nct01225/collatz_output/advanced_opt_%j.out
#SBATCH --error=/gpfs/scratch/nct_352/nct01225/collatz_output/advanced_opt_%j.err

# Advanced Compiler Optimization Test

echo "========================================="
echo "Advanced Compiler Optimization Test"
echo "========================================="

module purge
module load intel/2023.2
module load impi/2021.10

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PROC_BIND=close
export OMP_PLACES=cores

echo "Configuration:"
echo "  Advanced Intel compiler flags: -march=sapphirerapids -xCORE-AVX512"
echo "  Optimization level: -O3 -ipo -qopt-prefetch=4"
echo "  Memo table: 2^19 = 512K entries = 2MB (optimal from perf analysis)"
echo "  Threads: $OMP_NUM_THREADS per task"
echo ""

echo "=== Building with Advanced Optimizations ==="
bash build_gpp_advanced.sh

if [ $? -ne 0 ]; then
    echo "❌ Advanced build failed, falling back to standard build"
    bash build_gpp.sh
    EXECUTABLE="./collatz_mpi_gpp"
else
    echo "✅ Advanced build successful!"
    EXECUTABLE="./collatz_mpi_gpp_advanced"
fi

echo ""
echo "Using executable: $EXECUTABLE"

# Same workload as previous tests for comparison
START_OFFSET=0
COUNT=100000000
RUN_TAG="advanced_${SLURM_JOB_ID}"

echo "Workload: $COUNT seeds (same as baseline for comparison)"
echo ""

echo "=== Performance Test: Advanced vs Standard ==="

# Test 1: Advanced optimized version
echo "--- Test 1: Advanced Optimizations ---"
srun -n 2 --cpus-per-task=56 --cpu-bind=cores $EXECUTABLE $START_OFFSET $COUNT $RUN_TAG

echo ""
echo "Expected improvements:"
echo "  - Intel Sapphire Rapids optimizations: +5-10%"  
echo "  - AVX-512 vectorization: +5-15%"
echo "  - Prefetch optimizations: +3-8%"
echo "  - Total expected gain: +10-25% over baseline"
echo ""
echo "Baseline (2MB memo): 68.33 M/sec"
echo "Target (advanced):   75-85 M/sec"
echo ""
echo "========================================="
echo "Advanced test complete: $(date)"
echo "========================================="