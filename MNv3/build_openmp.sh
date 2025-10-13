#!/bin/bash
# Build script for V1.5-openmp on MareNostrum 5
# Target: GPP nodes with Intel Sapphire Rapids

set -e

echo "=== Building V1.5-openmp for MareNostrum 5 GPP ==="

# Load Intel compiler (adjust module names if needed)
# module load intel/2023.2.0  # Uncomment on MareNostrum

# Compilation flags optimized for Sapphire Rapids
CXX=${CXX:-g++}
CXXFLAGS="-O3 -march=native -mtune=native -fopenmp"
CXXFLAGS="$CXXFLAGS -Wall -Wextra"
CXXFLAGS="$CXXFLAGS -ffast-math -funroll-loops"

# For Intel compiler, use:
# CXX=icpx
# CXXFLAGS="-O3 -xHost -qopenmp -ipo -no-prec-div -fp-model fast=2"

OUTPUT="V1.5-openmp"

echo "Compiler: $CXX"
echo "Flags: $CXXFLAGS"
echo "Output: $OUTPUT"
echo ""

# Compile
$CXX $CXXFLAGS V1.5-openmp.cpp -o $OUTPUT

# Check result
if [ -f "$OUTPUT" ]; then
    echo "✓ Build successful: $OUTPUT"
    ls -lh $OUTPUT
    echo ""
    echo "To run locally (4 threads, small test):"
    echo "  ./$OUTPUT 0 1000000 --threads 4 --small-limit 20 --tag test"
    echo ""
    echo "To submit to MareNostrum:"
    echo "  sbatch slurm_openmp_1node.slurm"
else
    echo "✗ Build failed"
    exit 1
fi
