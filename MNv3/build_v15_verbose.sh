#!/bin/bash
# Diagnostic build script - verbose output

echo "=== Checking GCC and OpenMP ==="
module purge
module load gcc/13.2.0 openmpi/4.1.5-gcc

echo ""
echo "Compiler versions:"
g++ --version | head -1
mpicxx --show

echo ""
echo "=== Building with verbose flags ==="
set -x
mpicxx -v -O3 -march=native -flto -fno-exceptions -fno-rtti -funroll-loops -DNDEBUG -fopenmp \
       -o V1.5-mpi-lean-verbose V1.5-mpi-lean.cpp 2>&1 | tee build_verbose.log
set +x

echo ""
echo "=== Checking OpenMP linkage ==="
ldd V1.5-mpi-lean-verbose | grep -i omp

echo ""
echo "=== Checking symbols ==="
nm V1.5-mpi-lean-verbose | grep -i "omp_get_num_threads"

echo ""
if [ -f V1.5-mpi-lean-verbose ]; then
    echo "✅ Build SUCCESS: V1.5-mpi-lean-verbose"
    ls -lh V1.5-mpi-lean-verbose
else
    echo "❌ Build FAILED"
    exit 1
fi
