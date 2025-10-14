#!/bin/bash
#SBATCH --job-name=v15_overnight
#SBATCH --output=v15_overnight_%j.out
#SBATCH --error=v15_overnight_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=112
#SBATCH --time=02:00:00
#SBATCH --qos=acc_debug
#SBATCH --account=nct_352

# Load modules
module purge
module load gcc/13.2.0

# OpenMP configuration
export OMP_NUM_THREADS=112
export OMP_PLACES=cores
export OMP_PROC_BIND=spread

# Configuration
START_OFFSET=0
TOTAL_RANGE=3447600000000   # 3,447,600,000,000 numbers
                           # At 137M/sec = ~7 hours

echo "=== V1.5-openmp-checkpoint: Overnight Run with Progress Logging ==="
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo "Date: $(date)"
echo "Range: 2^71 + [$START_OFFSET, +$TOTAL_RANGE)"
echo "Expected time: ~7 hours (3,447,600,000,000 numbers at 137M/sec)"
echo "Checkpoint file: checkpoint_overnight_job${SLURM_JOB_ID}.log"
echo ""

# Run with progress logging (checkpoints every 1 billion numbers = ~7 seconds)
./V1.5-openmp-checkpoint $START_OFFSET $TOTAL_RANGE \
    --threads 112 \
    --tag overnight_job${SLURM_JOB_ID}

EXIT_CODE=$?

echo ""
echo "=== COMPLETED ==="
echo "Exit code: $EXIT_CODE"
echo "End time: $(date)"

# Create completion marker
if [ $EXIT_CODE -eq 0 ]; then
    echo "SUCCESS: Range [$START_OFFSET, +$TOTAL_RANGE) completed" > completion_${SLURM_JOB_ID}.txt
    echo "Completion time: $(date)" >> completion_${SLURM_JOB_ID}.txt
    echo "Checkpoint log: checkpoint_overnight_job${SLURM_JOB_ID}.log" >> completion_${SLURM_JOB_ID}.txt
fi
