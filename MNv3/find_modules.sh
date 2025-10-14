#!/bin/bash
# Find available GCC and MPI modules on MareNostrum 5

echo "=== Available GCC modules ==="
module avail gcc 2>&1 | grep -i gcc

echo ""
echo "=== Available MPI modules ==="
module avail 2>&1 | grep -i mpi

echo ""
echo "=== Recommended combinations ==="
echo "Option 1 (OpenMPI + GCC):"
echo "  module load gcc/13.1.0 openmpi/4.1.5"
echo ""
echo "Option 2 (Intel MPI + GCC backend):"
echo "  module load gcc impi/2021.10.0"
echo "  export I_MPI_CXX=g++"
echo ""
echo "Run this script on MareNostrum to see actual available versions"
