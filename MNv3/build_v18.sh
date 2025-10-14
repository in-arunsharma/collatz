#!/bin/bash
# Build V1.8-mpi-simple with proper flags

module load intel impi

echo "Building V1.8-mpi-simple..."
echo "Compiler: mpicxx (Intel MPI)"
echo ""

mpicxx -O3 -xHost -qopenmp -march=native -funroll-loops \
       -o V1.8-mpi-simple V1.8-mpi-simple.cpp

if [ $? -eq 0 ]; then
    echo ""
    echo "✓ Build successful!"
    ls -lh V1.8-mpi-simple
else
    echo ""
    echo "✗ Build failed!"
    exit 1
fi
