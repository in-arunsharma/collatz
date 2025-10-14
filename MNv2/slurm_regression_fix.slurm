#!/bin/bash
#SBATCH --job-name=mnv2_regression_fix
#SBATCH --account=nct_352
#SBATCH --qos=gp_debug
#SBATCH --partition=gpp
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=56
#SBATCH --time=00:10:00
#SBATCH --output=/gpfs/scratch/nct_352/nct01225/collatz_output/regression_fix_%j.out
#SBATCH --error=/gpfs/scratch/nct_352/nct01225/collatz_output/regression_fix_%j.err

# REGRESSION FIX: Revert broken optimizations, keep only safe ones

echo "========================================="
echo "REGRESSION FIX TEST"
echo "========================================="

module purge
module load intel/2023.2
module load impi/2021.10

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PROC_BIND=spread  
export OMP_PLACES=cores

echo "Configuration:"
echo "  REVERTED CHANGES:"
echo "  ❌ Peak tracking: Reverted to original (was breaking algorithm)"
echo "  ❌ Schedule: Reverted to static,10000 (static alone caused issues)"
echo "  ✅ Critical section: Kept minimal (safe optimization)"
echo "  ✅ Config: Fixed GPP_NODES=1 (was incorrectly set to 2)"
echo "  Threads: $OMP_NUM_THREADS per worker"
echo ""

# Same workload for comparison
START_OFFSET=0
COUNT=100000000
RUN_TAG="regression_fix_${SLURM_JOB_ID}"

echo "Workload: $COUNT seeds"
echo ""

echo "=== REGRESSION FIX TEST ==="

srun -n 2 --cpus-per-task=56 --cpu-bind=cores \
     ./collatz_mpi_gpp $START_OFFSET $COUNT $RUN_TAG

echo ""
echo "Performance Targets:"
echo "  Original working: 68.33 M/sec"  
echo "  Broken optimized: 32.23 M/sec"
echo "  Target (fixed):   65-70 M/sec (back to baseline)"
echo ""
echo "========================================="
echo "Regression fix complete: $(date)" 
echo "========================================="