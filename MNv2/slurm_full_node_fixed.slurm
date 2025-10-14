#!/bin/bash
#SBATCH --job-name=mnv2_full_node_fixed
#SBATCH --account=nct_352
#SBATCH --qos=gp_debug
#SBATCH --partition=gpp
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=3
#SBATCH --cpus-per-task=37
#SBATCH --time=00:10:00
#SBATCH --output=/gpfs/scratch/nct_352/nct01225/collatz_output/full_node_fixed_%j.out
#SBATCH --error=/gpfs/scratch/nct_352/nct01225/collatz_output/full_node_fixed_%j.err

# FIXED: Use 56 threads per task, but run 2 workers to use full node

echo "========================================="
echo "FULL NODE OPTIMIZATION TEST (FIXED)"
echo "========================================="

module purge
module load intel/2023.2
module load impi/2021.10

# Use 56 threads per MPI task (what SLURM gives us)
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PROC_BIND=spread  
export OMP_PLACES=cores

echo "Configuration:"
echo "  STRATEGY: 1 master + 2 workers × 37 threads = 74 total compute threads"
echo "  Peak tracking: Only on 3n+1 spikes" 
echo "  Schedule: static (no chunk overhead)"
echo "  Critical section: Minimal (no vectors)"
echo "  Threads per task: $OMP_NUM_THREADS"
echo "  Total compute threads: $(($OMP_NUM_THREADS * 2)) (2 workers)"
echo ""

# Same workload for comparison
START_OFFSET=0
COUNT=100000000
RUN_TAG="full_node_fixed_${SLURM_JOB_ID}"

echo "Workload: $COUNT seeds"
echo ""

echo "=== FULL NODE PERFORMANCE TEST (2 WORKERS) ==="

# Use 3 MPI ranks: 1 master + 2 workers
srun -n 3 --cpus-per-task=56 --cpu-bind=cores \
     ./collatz_mpi_gpp $START_OFFSET $COUNT $RUN_TAG

echo ""
echo "Performance Targets:"
echo "  Previous (1 worker, 56 threads): 68.33 M/sec"  
echo "  Expected (2 workers, 112 threads): 130-140 M/sec"
echo "  Theoretical max: ~2.2M/sec × 112 cores = 246 M/sec"
echo ""
echo "========================================="
echo "Full node test complete: $(date)" 
echo "========================================="