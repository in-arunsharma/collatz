#!/bin/bash
# rebuild_v16b.sh - Rebuild with correct modules on MareNostrum

echo "=== Rebuilding V1.6b with Intel MPI (default on MN5) ==="

# Use default Intel MPI that's already loaded
echo "Using Intel oneAPI compilers (icpx/mpicxx)"
echo ""

# Build with Intel MPI
echo "Building V1.6b-mpi-openmp..."
module load intel impi
mpicxx -O3 -xHost -qopenmp -o V1.6b-mpi-openmp V1.6b-mpi-openmp.cpp

if [ $? -eq 0 ]; then
    echo "✅ Build successful!"
    ls -lh V1.6b-mpi-openmp
    echo ""
    echo "Note: Built with Intel MPI (default on MareNostrum 5)"
else
    echo "❌ Build failed!"
    exit 1
fi
