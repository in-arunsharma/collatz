#!/bin/bash
# Local MPI testing script for Phase 1 GPP code
# Run this before submitting to MareNostrum to verify everything works

set -e  # Exit on error

echo "========================================="
echo "Collatz Phase 1 - Local MPI Test"
echo "========================================="

# Check for build
if [ ! -f ./collatz_mpi_gpp ]; then
    echo "Executable not found. Building..."
    ./build_phase1_gpp.sh
fi

# Test parameters (small range for local testing)
START_SEED=0
END_SEED=1000000  # 1 million numbers (quick test)
OUTPUT_PREFIX="local_test"

# Determine number of MPI ranks to use
if [ $# -gt 0 ]; then
    NUM_RANKS=$1
else
    # Auto-detect: use 4 ranks or number of physical cores, whichever is smaller
    NUM_CORES=$(nproc)
    NUM_RANKS=$((NUM_CORES < 4 ? NUM_CORES : 4))
fi

# Set OpenMP threads per rank
# For local testing: divide available cores among ranks
THREADS_PER_RANK=$((NUM_CORES / NUM_RANKS))
export OMP_NUM_THREADS=$THREADS_PER_RANK
export OMP_PROC_BIND=close
export OMP_PLACES=cores

echo ""
echo "Test Configuration:"
echo "  MPI ranks: $NUM_RANKS"
echo "  OpenMP threads per rank: $OMP_NUM_THREADS"
echo "  Total parallelism: $((NUM_RANKS * OMP_NUM_THREADS)) threads"
echo "  Seed range: $START_SEED to $END_SEED"
echo "  Available cores: $NUM_CORES"
echo ""

# Clean old test files
rm -f ${OUTPUT_PREFIX}*.json

echo "========================================="
echo "Running MPI test..."
echo "========================================="

# Run with mpirun (adjust for your MPI installation)
mpirun -np $NUM_RANKS ./collatz_mpi_gpp $START_SEED $END_SEED $OUTPUT_PREFIX

EXIT_CODE=$?

echo ""
echo "========================================="
echo "Test Results"
echo "========================================="

if [ $EXIT_CODE -eq 0 ]; then
    echo "✓ Test PASSED"
    echo ""
    echo "Generated files:"
    ls -lh ${OUTPUT_PREFIX}*.json 2>/dev/null || echo "No output files found"
    echo ""
    echo "Ready for MareNostrum deployment!"
    echo "Next steps:"
    echo "  1. Copy code to MareNostrum: scp -r marenostrum/ <username>@mn5.bsc.es:~/collatz/"
    echo "  2. SSH to MareNostrum: ssh <username>@mn5.bsc.es"
    echo "  3. Build there: cd collatz/marenostrum && ./build_phase1_gpp.sh"
    echo "  4. Submit job: sbatch slurm_phase1_gpp.slurm"
else
    echo "✗ Test FAILED (exit code: $EXIT_CODE)"
    echo "Check error messages above"
    exit 1
fi

echo "========================================="
