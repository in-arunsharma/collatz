#!/bin/bash
#SBATCH --job-name=v15_array2
#SBATCH --output=v15_array2_%A_%a.out
#SBATCH --error=v15_array2_%A_%a.err
#SBATCH --array=0-1               # 2 tasks = 2 nodes
#SBATCH --nodes=1                 # Each task gets 1 node
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=112
#SBATCH --time=00:10:00
#SBATCH --qos=gp_debug
#SBATCH --account=nct_352

echo "=== V1.5 Array Job - Task $SLURM_ARRAY_TASK_ID ==="
echo "Date: $(date)"
echo "Node: $(hostname)"
echo "Job: ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
echo ""

# Split 100M range across 2 tasks
TOTAL_RANGE=100000000
TASKS=2
RANGE_PER_TASK=$((TOTAL_RANGE / TASKS))
START_OFFSET=$((SLURM_ARRAY_TASK_ID * RANGE_PER_TASK))

echo "Work Distribution:"
echo "  Task ID: $SLURM_ARRAY_TASK_ID"
echo "  Start:   $START_OFFSET"
echo "  Count:   $RANGE_PER_TASK"
echo ""

# OpenMP configuration (same as proven baseline)
export OMP_NUM_THREADS=112
export OMP_PROC_BIND=spread
export OMP_PLACES=cores

echo "Configuration:"
echo "  Threads: $OMP_NUM_THREADS"
echo "  Binding: $OMP_PROC_BIND"
echo ""

time ./V1.5-openmp $START_OFFSET $RANGE_PER_TASK --threads 112 --small-limit 20 --tag "task${SLURM_ARRAY_TASK_ID}"

echo ""
echo "=== Task $SLURM_ARRAY_TASK_ID Complete ==="
echo "Expected per-node: 137M nums/sec"
echo "Expected total (2 nodes): 274M nums/sec"
