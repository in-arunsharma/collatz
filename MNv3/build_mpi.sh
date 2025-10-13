#!/bin/bash
# Build script for V1.6-mpi-openmp on MareNostrum 5
# Requires MPI + OpenMP support

set -e

echo "=== Building V1.6-mpi-openmp for MareNostrum 5 GPP Multi-Node ==="

# Load modules (adjust for MareNostrum 5)
# module load intel/2023.2.0
# module load impi/2021.9.0

# Compilation
MPICXX=${MPICXX:-mpicxx}
CXXFLAGS="-O3 -march=native -mtune=native -fopenmp"
CXXFLAGS="$CXXFLAGS -Wall -Wextra"
CXXFLAGS="$CXXFLAGS -ffast-math -funroll-loops"

# For Intel MPI compiler:
# MPICXX=mpiicpx
# CXXFLAGS="-O3 -xHost -qopenmp -ipo -no-prec-div -fp-model fast=2"

OUTPUT="V1.6-mpi-openmp"

echo "MPI Compiler: $MPICXX"
echo "Flags: $CXXFLAGS"
echo "Output: $OUTPUT"
echo ""

# Compile
$MPICXX $CXXFLAGS V1.6-mpi-openmp.cpp -o $OUTPUT

# Check result
if [ -f "$OUTPUT" ]; then
    echo "✓ Build successful: $OUTPUT"
    ls -lh $OUTPUT
    echo ""
    echo "To test locally (2 ranks × 4 threads):"
    echo "  mpirun -np 2 ./$OUTPUT 0 1000000 --threads 4 --tag local_mpi"
    echo ""
    echo "To submit to MareNostrum (10 nodes):"
    echo "  sbatch slurm_mpi_10nodes.slurm"
else
    echo "✗ Build failed"
    exit 1
fi
