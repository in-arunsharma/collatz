#!/bin/bash
#SBATCH --job-name=collatz_omp
#SBATCH --output=collatz_omp_%j.out
#SBATCH --error=collatz_omp_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=112
#SBATCH --time=02:00:00
#SBATCH --qos=gp_debug
#SBATCH --account=nct_352

# MareNostrum 5 GPP Node - OpenMP Single Node Test
# 1 node × 112 cores = 112 threads
# Expected: 100-200M nums/sec (50-100× speedup)

echo "========================================="
echo "Collatz V1.5-openmp - Single Node Test"
echo "========================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo "Cores: $SLURM_CPUS_PER_TASK"
echo "Date: $(date)"
echo ""

# Load modules (adjust as needed for MareNostrum 5)
# module load intel/2023.2.0

# Set OpenMP environment
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PROC_BIND=close
export OMP_PLACES=cores

echo "OpenMP Configuration:"
echo "  OMP_NUM_THREADS=$OMP_NUM_THREADS"
echo "  OMP_PROC_BIND=$OMP_PROC_BIND"
echo "  OMP_PLACES=$OMP_PLACES"
echo ""

# Test parameters
START_OFFSET=0
COUNT=100000000  # 100M numbers (about 1 minute with 112 cores)
MEMO_BITS=20     # 2^20 = 1M entries (4 MB - fits in L3)
TAG="mn5_1node_${SLURM_JOB_ID}"

echo "Test Parameters:"
echo "  Start offset: $START_OFFSET"
echo "  Count: $COUNT"
echo "  Memo bits: $MEMO_BITS"
echo "  Tag: $TAG"
echo ""

# Run the program
echo "Starting computation..."
echo ""

time ./V1.5-openmp $START_OFFSET $COUNT \
    --threads $OMP_NUM_THREADS \
    --small-limit $MEMO_BITS \
    --tag $TAG

echo ""
echo "========================================="
echo "Job completed at $(date)"
echo "========================================="

# Show output files
echo ""
echo "Generated files:"
ls -lh *${TAG}* 2>/dev/null || echo "No output files found"
