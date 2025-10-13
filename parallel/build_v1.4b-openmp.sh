#!/bin/bash
# Build script for V1.4b-openmp (OpenMP parallelization)

set -e  # Exit on error

echo "=== Building V1.4b-openmp (OpenMP) ==="

# Get git hash for reproducibility
GIT_HASH=$(git rev-parse --short HEAD 2>/dev/null || echo "unknown")
echo "Git commit: $GIT_HASH"

# Compiler flags
CXX=${CXX:-g++}
CXXFLAGS="-O3 -march=native -std=c++17 -fopenmp -flto -fno-exceptions -fno-rtti -funroll-loops -fno-asynchronous-unwind-tables -DNDEBUG -DGIT_HASH='\"$GIT_HASH\"'"

echo "Compiler: $CXX"
echo "Flags: $CXXFLAGS"

# Build
$CXX $CXXFLAGS V1.4b-openmp.cpp -o V1.4b-openmp

echo "✅ Build successful: V1.4b-openmp"
echo ""
echo "Usage examples:"
echo "  OMP_NUM_THREADS=1 ./V1.4b-openmp 0 1000000  # Single-threaded (verify determinism)"
echo "  OMP_NUM_THREADS=2 ./V1.4b-openmp 0 1000000  # 2 threads"
echo "  OMP_NUM_THREADS=4 ./V1.4b-openmp 0 1000000  # 4 threads"
echo "  OMP_NUM_THREADS=8 ./V1.4b-openmp 0 10000000 # 8 threads, larger workload"
echo ""
echo "MareNostrum 5 deployment:"
echo "  export OMP_NUM_THREADS=80"
echo "  export OMP_PROC_BIND=close"
echo "  export OMP_PLACES=cores"
echo "  ./V1.4b-openmp 0 100000000 --tag mn5_run1"
