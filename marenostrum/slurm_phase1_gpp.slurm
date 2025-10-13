#!/bin/bash
#SBATCH --job-name=collatz_phase1_gpp
#SBATCH --partition=gpp
#SBATCH --nodes=2
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=112
#SBATCH --time=00:15:00
#SBATCH --output=collatz_phase1_%j.out
#SBATCH --error=collatz_phase1_%j.err
#SBATCH --exclusive

# MareNostrum 5 - Phase 1: GPP Nodes Only (MPI + OpenMP)
# Shared resource pool: 150 GPP nodes for all participants
# Using 2 nodes for initial validation (very conservative, 1.3% of pool)
# Expected throughput: ~500M numbers/sec (2 nodes × 112 cores × 2.3M/sec)
# Scale up: 2 → 3 → 5 → 10 nodes as validation succeeds

echo "========================================="
echo "Collatz Phase 1 - GPP Nodes Only"
echo "========================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Nodes: $SLURM_JOB_NUM_NODES"
echo "MPI Tasks: $SLURM_NTASKS"
echo "CPUs per task: $SLURM_CPUS_PER_TASK"
echo "Start time: $(date)"
echo "========================================="

# Load required modules
module purge
module load intel/2023.2
module load openmpi/4.1.5

# Environment setup
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PROC_BIND=close
export OMP_PLACES=cores

# Show configuration
echo ""
echo "OpenMP Configuration:"
echo "  OMP_NUM_THREADS: $OMP_NUM_THREADS"
echo "  OMP_PROC_BIND: $OMP_PROC_BIND"
echo "  OMP_PLACES: $OMP_PLACES"
echo ""

# Node information
echo "Allocated nodes:"
scontrol show hostname $SLURM_JOB_NODELIST
echo ""

# Change to working directory
cd $SLURM_SUBMIT_DIR

# Verify executable exists
if [ ! -f ./collatz_mpi_gpp ]; then
    echo "ERROR: Executable 'collatz_mpi_gpp' not found!"
    echo "Please run build_phase1_gpp.sh first."
    exit 1
fi

# Test parameters (adjust as needed)
START_SEED=0
END_SEED=100000000  # 100 million for initial validation (ultra-fast test)
OUTPUT_PREFIX="phase1_gpp_validation"

echo "========================================="
echo "Running Collatz MPI computation"
echo "========================================="
echo "Seed range: $START_SEED to $END_SEED"
echo "Output prefix: $OUTPUT_PREFIX"
echo ""

# Run the MPI program
# Each MPI rank will spawn $OMP_NUM_THREADS threads
# Total parallelism: 2 nodes × 112 threads = 224 parallel threads
mpirun -np $SLURM_NTASKS \
       --map-by ppr:1:node \
       --bind-to none \
       ./collatz_mpi_gpp $START_SEED $END_SEED $OUTPUT_PREFIX

EXIT_CODE=$?

echo ""
echo "========================================="
echo "Job completed"
echo "========================================="
echo "Exit code: $EXIT_CODE"
echo "End time: $(date)"
echo "========================================="

# Show output files
echo ""
echo "Generated files:"
ls -lh ${OUTPUT_PREFIX}*.json 2>/dev/null || echo "No output files found"

exit $EXIT_CODE
