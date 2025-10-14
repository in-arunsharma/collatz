#!/bin/bash
# Monitor progress of batch jobs

echo "=== Batch Job Monitor ==="
echo ""

# Find the job ID from output files
JOB_ID=$(ls v15_batch_*.out 2>/dev/null | head -1 | sed 's/v15_batch_\([0-9]*\)_.*/\1/')

if [ -z "$JOB_ID" ]; then
    echo "No batch jobs found. Submit with: sbatch run_batch_overnight.slurm"
    exit 1
fi

echo "Batch Job ID: $JOB_ID"
echo ""

# Check job status
echo "--- Job Status ---"
squeue -j $JOB_ID

echo ""
echo "--- Completed Jobs ---"
ls completion_batch*_${JOB_ID}.txt 2>/dev/null | wc -l | xargs echo "Completed:"
echo ""

# Show latest checkpoint from each job
echo "--- Latest Checkpoints ---"
for i in {0..4}; do
    CHECKPOINT="checkpoint_batch${i}_job${JOB_ID}.log"
    if [ -f "$CHECKPOINT" ]; then
        echo "Job $i:"
        tail -1 "$CHECKPOINT" 2>/dev/null || echo "  No checkpoints yet"
    fi
done

echo ""
echo "--- Quick Summary ---"
TOTAL_COMPLETED=$(cat completion_batch*_${JOB_ID}.txt 2>/dev/null | wc -l)
echo "Jobs completed: $TOTAL_COMPLETED / 5"
echo "Numbers tested: $((TOTAL_COMPLETED * 980000000000)) / 4900000000000"
echo ""
echo "To see live updates from a specific job:"
echo "  tail -f checkpoint_batch0_job${JOB_ID}.log"
