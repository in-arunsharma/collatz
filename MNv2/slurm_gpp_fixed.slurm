#!/bin/bash
#SBATCH --job-name=mnv2_gpp_phase1
#SBATCH --account=nct_352
#SBATCH --qos=gp_debug
#SBATCH --partition=gpp
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=112
#SBATCH --time=00:15:00
#SBATCH --output=/gpfs/scratch/nct_352/nct01225/collatz_output/mnv2_gpp_%j.out
#SBATCH --error=/gpfs/scratch/nct_352/nct01225/collatz_output/mnv2_gpp_%j.err
#SBATCH --exclusive

# MNv2 Phase 1: GPP Nodes (MPI + OpenMP)
# 
# Configuration: 2 GPP nodes, 1 task per node (worker), 112 threads per task
# Total: 2 workers + 1 master (master runs on login/first node)
# Expected throughput: ~500M nums/sec

echo "========================================="
echo "MNv2 Phase 1 - GPP Nodes"
echo "========================================="
echo "Job ID:       $SLURM_JOB_ID"
echo "Nodes:        $SLURM_JOB_NUM_NODES"
echo "MPI tasks:    $SLURM_NTASKS"
echo "CPUs/task:    $SLURM_CPUS_PER_TASK"
echo "Partition:    $SLURM_JOB_PARTITION"
echo "Start time:   $(date)"
echo "========================================="
echo ""

# Load modules (MareNostrum 5 specific)
module purge
module load intel/2023.2
module load impi/2021.10

# OpenMP configuration (critical for performance!)
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PROC_BIND=close
export OMP_PLACES=cores
export OMP_DISPLAY_ENV=verbose

echo "OpenMP Configuration:"
echo "  OMP_NUM_THREADS:  $OMP_NUM_THREADS"
echo "  OMP_PROC_BIND:    $OMP_PROC_BIND"
echo "  OMP_PLACES:       $OMP_PLACES"
echo ""

# Node information
echo "Allocated nodes:"
scontrol show hostname $SLURM_JOB_NODELIST
echo ""

# MPI configuration: 3 total ranks (1 master + 2 workers)
TOTAL_RANKS=3

echo "MPI Configuration:"
echo "  Total ranks:  $TOTAL_RANKS (1 master + 2 workers)"
echo ""

# Workload configuration
START_OFFSET=0
COUNT=1000000000         # 1 billion seeds (~2 seconds on 2 nodes)
RUN_TAG="mnv2_gpp_${SLURM_JOB_ID}"

echo "Workload:"
echo "  Base:         2^71"
echo "  Start offset: $START_OFFSET"
echo "  Count:        $COUNT"
echo "  Run tag:      $RUN_TAG"
echo ""

# Run the MPI+OpenMP executable
echo "========================================="
echo "Starting MPI+OpenMP execution..."
echo "========================================="
echo ""

# Use srun to launch MPI job across allocated nodes
srun -n $TOTAL_RANKS ./collatz_mpi_gpp $START_OFFSET $COUNT $RUN_TAG

EXIT_CODE=$?

echo ""
echo "========================================="
if [ $EXIT_CODE -eq 0 ]; then
    echo "✅ JOB COMPLETED SUCCESSFULLY"
else
    echo "❌ JOB FAILED (exit code: $EXIT_CODE)"
fi
echo "========================================="
echo "End time: $(date)"
echo ""

# Performance summary
echo "Performance Summary:"
sacct -j $SLURM_JOB_ID --format=JobID,JobName,Partition,AllocNodes,State,ExitCode,Elapsed,CPUTime,MaxRSS,MaxVMSize

exit $EXIT_CODE
