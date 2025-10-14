#!/bin/bash
#SBATCH --job-name=mnv2_perf_stat
#SBATCH --account=nct_352
#SBATCH --qos=gp_debug
#SBATCH --partition=gpp
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=56
#SBATCH --time=00:10:00
#SBATCH --output=/gpfs/scratch/nct_352/nct01225/collatz_output/perf_stat_%j.out
#SBATCH --error=/gpfs/scratch/nct_352/nct01225/collatz_output/perf_stat_%j.err

# Detailed performance counter analysis

echo "========================================="
echo "Performance Counter Analysis"
echo "========================================="

module purge
module load intel/2023.2
module load impi/2021.10

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PROC_BIND=close
export OMP_PLACES=cores

echo "Configuration:"
echo "  Memo table: 2^19 = 512K entries = 2MB"
echo "  Threads: $OMP_NUM_THREADS per task"
echo ""

# Smaller workload for accurate measurement
START_OFFSET=0
COUNT=50000000  # 50M seeds for good measurement
RUN_TAG="perf_stat_${SLURM_JOB_ID}"

echo "Workload: $COUNT seeds (optimized for measurement)"
echo ""

echo "=== Detailed Performance Counters ==="

# Run with comprehensive performance counters
perf stat -d -d -d \
  -e cycles,instructions,cache-references,cache-misses \
  -e L1-dcache-loads,L1-dcache-load-misses \
  -e LLC-loads,LLC-load-misses \
  -e node-loads,node-load-misses \
  -e dTLB-loads,dTLB-load-misses \
  -e branch-loads,branch-load-misses \
  -e stalled-cycles-frontend,stalled-cycles-backend \
  -e context-switches,cpu-migrations \
  -e page-faults,minor-faults,major-faults \
  srun -n 2 --cpus-per-task=56 --cpu-bind=cores ./collatz_mpi_gpp $START_OFFSET $COUNT $RUN_TAG 2>&1

echo ""
echo "=== Memory Bandwidth Analysis ==="

# Additional memory-specific counters
perf stat \
  -e mem_load_retired.l1_hit \
  -e mem_load_retired.l1_miss \
  -e mem_load_retired.l2_hit \
  -e mem_load_retired.l2_miss \
  -e mem_load_retired.l3_hit \
  -e mem_load_retired.l3_miss \
  -e mem_inst_retired.all_loads \
  -e mem_inst_retired.all_stores \
  srun -n 2 --cpus-per-task=56 --cpu-bind=cores ./collatz_mpi_gpp $START_OFFSET $COUNT ${RUN_TAG}_mem 2>&1

echo ""
echo "=== CPU Utilization Breakdown ==="

# CPU efficiency counters
perf stat \
  -e cpu_clk_unhalted.thread \
  -e cpu_clk_unhalted.ref_tsc \
  -e uops_issued.any \
  -e uops_retired.retire_slots \
  -e resource_stalls.any \
  -e resource_stalls.memory \
  srun -n 2 --cpus-per-task=56 --cpu-bind=cores ./collatz_mpi_gpp $START_OFFSET $COUNT ${RUN_TAG}_cpu 2>&1

echo ""
echo "Performance analysis complete."
echo ""
echo "========================================="
echo "Stat complete: $(date)"
echo "========================================="