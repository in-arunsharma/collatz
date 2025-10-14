#!/bin/bash
# Local MPI testing script for MNv2

set -e

echo "========================================="
echo "MNv2: Local MPI Testing"
echo "========================================="
echo ""

# Check if executable exists
if [ ! -f "./collatz_mpi_gpp" ]; then
    echo "[ERROR] Executable not found: ./collatz_mpi_gpp"
    echo "[INFO] Run ./build_gpp.sh first"
    exit 1
fi

# Check if mpirun is available
if ! command -v mpirun &> /dev/null; then
    echo "[ERROR] mpirun not found!"
    echo "[INFO] Install OpenMPI: sudo apt install openmpi-bin libopenmpi-dev"
    exit 1
fi

# Test configuration
NUM_RANKS=${1:-4}        # Default: 4 ranks (1 master + 3 workers)
START_OFFSET=0
COUNT=${2:-10000000}     # Default: 10M seeds (realistic test, not too small!)
RUN_TAG="local_test"

echo "Test Configuration:"
echo "  MPI ranks:    $NUM_RANKS (1 master + $((NUM_RANKS - 1)) workers)"
echo "  Start offset: $START_OFFSET"
echo "  Count:        $COUNT seeds"
echo "  Run tag:      $RUN_TAG"
echo ""
echo "⚠️  Note: MPI has ~50ms overhead. Use ≥10M seeds for realistic throughput!"
echo ""

# Set OpenMP threads per MPI rank (conservative for local testing)
THREADS_PER_RANK=${OMP_NUM_THREADS:-4}
export OMP_NUM_THREADS=$THREADS_PER_RANK
export OMP_PROC_BIND=close
export OMP_PLACES=cores

echo "OpenMP Configuration:"
echo "  Threads/rank: $THREADS_PER_RANK"
echo "  Proc bind:    $OMP_PROC_BIND"
echo "  Places:       $OMP_PLACES"
echo ""

echo "========================================="
echo "Running MPI test..."
echo "========================================="
echo ""

# Run MPI job
mpirun -np $NUM_RANKS ./collatz_mpi_gpp $START_OFFSET $COUNT $RUN_TAG

EXIT_CODE=$?

echo ""
echo "========================================="
if [ $EXIT_CODE -eq 0 ]; then
    echo "✅ TEST PASSED"
else
    echo "❌ TEST FAILED (exit code: $EXIT_CODE)"
fi
echo "========================================="
echo ""

exit $EXIT_CODE
