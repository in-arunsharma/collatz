#!/bin/bash
#SBATCH --job-name=mnv2_full_node
#SBATCH --account=nct_352
#SBATCH --qos=gp_debug
#SBATCH --partition=gpp
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=112
#SBATCH --ntasks-per-node=2
#SBATCH --time=00:10:00
#SBATCH --output=/gpfs/scratch/nct_352/nct01225/collatz_output/full_node_%j.out
#SBATCH --error=/gpfs/scratch/nct_352/nct01225/collatz_output/full_node_%j.err

# FULL NODE TEST: 1 master + 1 worker with 112 threads

echo "========================================="
echo "FULL NODE OPTIMIZATION TEST"
echo "========================================="

module purge
module load intel/2023.2
module load impi/2021.10

# CRITICAL: Let worker use ALL 112 cores
export OMP_NUM_THREADS=112
export OMP_PROC_BIND=spread  
export OMP_PLACES=cores

echo "Configuration:"
echo "  OPTIMIZATION LEVEL: MAXIMUM"
echo "  Thread cap removed: Using full 112 threads"
echo "  Peak tracking: Only on 3n+1 spikes" 
echo "  Schedule: static (no chunk overhead)"
echo "  Critical section: Minimal (no vectors)"
echo "  Threads: $OMP_NUM_THREADS (FULL NODE!)"
echo ""

# Same workload for comparison
START_OFFSET=0
COUNT=100000000
RUN_TAG="full_node_${SLURM_JOB_ID}"

echo "Workload: $COUNT seeds"
echo ""

echo "=== FULL NODE PERFORMANCE TEST ==="

# Use optimal CPU binding and memory policy
srun -n 2 --cpus-per-task=112 --ntasks-per-node=2 --cpu-bind=cores \
     numactl --interleave=all \
     ./collatz_mpi_gpp $START_OFFSET $COUNT $RUN_TAG

echo ""
echo "Performance Targets:"
echo "  Previous (56 threads): 68.33 M/sec"  
echo "  Expected (112 threads): 200-250 M/sec"
echo "  Theoretical max: ~2.2M/sec × 112 cores = 246 M/sec"
echo ""
echo "========================================="
echo "Full node test complete: $(date)" 
echo "========================================="