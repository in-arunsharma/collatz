#!/bin/bash
#SBATCH --job-name=mnv2_profile
#SBATCH --account=nct_352
#SBATCH --qos=gp_debug
#SBATCH --partition=gpp
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=56
#SBATCH --time=00:10:00
#SBATCH --output=/gpfs/scratch/nct_352/nct01225/collatz_output/profile_%j.out
#SBATCH --error=/gpfs/scratch/nct_352/nct01225/collatz_output/profile_%j.err

# Performance profiling test with perf

echo "========================================="
echo "Performance Profiling Test"
echo "========================================="

module purge
module load intel/2023.2
module load impi/2021.10

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PROC_BIND=close
export OMP_PLACES=cores

echo "Configuration:"
echo "  Threads:   $OMP_NUM_THREADS per task"
echo "  Profiling: perf record + report"
echo ""

# Smaller workload for profiling (10M seeds instead of 100M)
START_OFFSET=0
COUNT=10000000
RUN_TAG="profile_${SLURM_JOB_ID}"

echo "Workload: $COUNT seeds (smaller for profiling)"
echo ""

# Record performance counters during execution
echo "=== Running with perf record ==="
perf record -g -e cycles,instructions,cache-misses,cache-references,L1-dcache-load-misses,LLC-load-misses \
    srun -n 2 --cpus-per-task=56 --cpu-bind=cores ./collatz_mpi_gpp $START_OFFSET $COUNT $RUN_TAG

echo ""
echo "=== Performance Report ==="
perf report --stdio > /gpfs/scratch/nct_352/nct01225/collatz_output/perf_report_${SLURM_JOB_ID}.txt

echo ""
echo "=== Performance Statistics ==="
perf stat -e cycles,instructions,cache-misses,cache-references,L1-dcache-load-misses,LLC-load-misses,stalled-cycles-frontend,stalled-cycles-backend \
    srun -n 2 --cpus-per-task=56 --cpu-bind=cores ./collatz_mpi_gpp $START_OFFSET $COUNT ${RUN_TAG}_stat

echo ""
echo "Performance report saved to: /gpfs/scratch/nct_352/nct01225/collatz_output/perf_report_${SLURM_JOB_ID}.txt"
echo ""
echo "========================================="
echo "Profile complete: $(date)"
echo "========================================="