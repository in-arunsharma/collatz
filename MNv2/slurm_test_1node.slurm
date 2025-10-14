#!/bin/bash
#SBATCH --job-name=mnv2_test_1node
#SBATCH --account=nct_352
#SBATCH --qos=gp_debug
#SBATCH --partition=gpp
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=56
#SBATCH --time=00:10:00
#SBATCH --output=/gpfs/scratch/nct_352/nct01225/collatz_output/mnv2_test_%j.out
#SBATCH --error=/gpfs/scratch/nct_352/nct01225/collatz_output/mnv2_test_%j.err

# Simple 1-node test: 1 master + 1 worker
# Each task gets 56 threads (half a node)

echo "========================================="
echo "MNv2 Test - 1 Node"
echo "========================================="
echo "Job ID:       $SLURM_JOB_ID"
echo "Start time:   $(date)"
echo ""

module purge
module load intel/2023.2
module load impi/2021.10

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PROC_BIND=close
export OMP_PLACES=cores

echo "Configuration:"
echo "  Nodes:     1"
echo "  Tasks:     2 (1 master + 1 worker)"
echo "  Threads:   $OMP_NUM_THREADS per task"
echo ""

# Test with 100M seeds
START_OFFSET=0
COUNT=100000000
RUN_TAG="test_${SLURM_JOB_ID}"

echo "Workload: $COUNT seeds"
echo ""

srun -n 2 --cpus-per-task=56 --cpu-bind=cores ./collatz_mpi_gpp $START_OFFSET $COUNT $RUN_TAG

echo ""
echo "========================================="
echo "Test complete: $(date)"
echo "========================================="
