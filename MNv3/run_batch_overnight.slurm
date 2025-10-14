#!/bin/bash
#SBATCH --job-name=v15_batch
#SBATCH --output=v15_batch_%A_%a.out
#SBATCH --error=v15_batch_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=112
#SBATCH --time=02:00:00
#SBATCH --qos=acc_debug
#SBATCH --account=nct_352
#SBATCH --array=0-35    # 36 jobs × 2 hours = 72 hours of work (~3 days if sequential)

# Load modules
module purge
module load gcc/13.2.0

# OpenMP configuration
export OMP_NUM_THREADS=112
export OMP_PLACES=cores
export OMP_PROC_BIND=spread

# Configuration - Each job processes 980 billion numbers (~2 hours each)
CHUNK_SIZE=980000000000      # 980 billion per job
START_OFFSET=$((SLURM_ARRAY_TASK_ID * CHUNK_SIZE))
TOTAL_RANGE=$CHUNK_SIZE

echo "=== V1.5-openmp-checkpoint: Batch Job ${SLURM_ARRAY_TASK_ID} of 36 ==="
echo "Job ID: ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
echo "Node: $SLURM_NODELIST"
echo "Date: $(date)"
echo "Range: 2^71 + [$START_OFFSET, +$TOTAL_RANGE)"
echo "Expected time: ~2 hours"
echo "Checkpoint file: checkpoint_batch${SLURM_ARRAY_TASK_ID}_job${SLURM_ARRAY_JOB_ID}.log"
echo ""
echo "Progress: Job $((SLURM_ARRAY_TASK_ID + 1))/36"
echo "Total coverage: 0 to $((36 * CHUNK_SIZE)) (35.28 trillion numbers)"
echo ""

# Run with progress logging
./V1.5-openmp-checkpoint $START_OFFSET $TOTAL_RANGE \
    --threads 112 \
    --tag batch${SLURM_ARRAY_TASK_ID}_job${SLURM_ARRAY_JOB_ID}

EXIT_CODE=$?

echo ""
echo "=== COMPLETED ==="
echo "Batch: ${SLURM_ARRAY_TASK_ID}/5"
echo "Range: [$START_OFFSET, +$TOTAL_RANGE)"
echo "Exit code: $EXIT_CODE"
echo "End time: $(date)"

# Create completion marker
if [ $EXIT_CODE -eq 0 ]; then
    echo "SUCCESS: Batch ${SLURM_ARRAY_TASK_ID} - Range [$START_OFFSET, +$TOTAL_RANGE) completed" > completion_batch${SLURM_ARRAY_TASK_ID}_${SLURM_ARRAY_JOB_ID}.txt
    echo "Completion time: $(date)" >> completion_batch${SLURM_ARRAY_TASK_ID}_${SLURM_ARRAY_JOB_ID}.txt
    echo "Checkpoint log: checkpoint_batch${SLURM_ARRAY_TASK_ID}_job${SLURM_ARRAY_JOB_ID}.log" >> completion_batch${SLURM_ARRAY_TASK_ID}_${SLURM_ARRAY_JOB_ID}.txt
fi
