#!/bin/bash
#SBATCH --job-name=mnv2_gpp_phase1
#SBATCH --account=nct_352
#SBATCH --qos=gp_debug
#SBATCH --partition=gpp
#SBATCH --nodes=2
#SBATCH --ntasks=3
#SBATCH --cpus-per-task=112
#SBATCH --time=00:15:00
#SBATCH --output=/gpfs/scratch/nct_352/nct01225/collatz_output/mnv2_gpp_%j.out
#SBATCH --error=/gpfs/scratch/nct_352/nct01225/collatz_output/mnv2_gpp_%j.err
#SBATCH --exclusive

# MNv2 Phase 1: GPP Nodes (MPI + OpenMP)
# 
# HACKATHON TEAM: Modify these SBATCH directives as needed:
# - --nodes:          Number of GPP nodes (matches config.hpp GPP_NODES)
# - --ntasks:         nodes + 1 (1 master + N workers)
# - --cpus-per-task:  112 (GPP cores) or 80 (ACC cores)
# - --time:           Walltime limit (00:15:00 = 15 minutes for testing)
# - --output/--error: Put logs in /gpfs/scratch (temporary, no backup)
#
# FILESYSTEM USAGE:
# - Source code:  /gpfs/projects/nct_352/$USER/collatz/MNv2/  (backed up)
# - Job outputs:  /gpfs/scratch/nct_352/$USER/collatz_output/ (temporary)
# - Temp files:   $TMPDIR (local SSD, auto-cleaned)
#
# Current config: 2 nodes = 1 master + 2 workers
# Expected throughput: ~500M nums/sec (2 × 112 cores × 2.2M/sec)

# Create output directory if needed
mkdir -p /gpfs/scratch/nct_352/nct01225/collatz_output

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

# MPI configuration
echo "MPI Configuration:"
echo "  Master rank:  0 (coordinator)"
echo "  Worker ranks: 1..$((SLURM_NTASKS - 1))"
echo ""

# Workload configuration (adjust for your run)
START_OFFSET=0           # Offset from 2^71
COUNT=1000000000         # 1 billion seeds (adjust as needed)
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

mpirun -np $SLURM_NTASKS ./collatz_mpi_gpp $START_OFFSET $COUNT $RUN_TAG

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
