#!/bin/bash
# Build V1.5-mpi-lean with GCC (NOT icpc - icpc kills __int128 performance!)
# Friend's recommended flags: -O3 -march=native -flto -fno-exceptions -fno-rtti -funroll-loops -DNDEBUG

# MareNostrum 5 modules - OPTION 1: OpenMPI + GCC (RECOMMENDED)
module purge
module load gcc/13.2.0 openmpi/4.1.5-gcc

# OPTION 2: Intel MPI + GCC backend (uncomment if preferred)
# module purge
# module load gcc/13.2.0 impi/2021.10.0
# export I_MPI_CXX=g++

# Build with GCC
mpicxx -O3 -march=native -flto -fno-exceptions -fno-rtti -funroll-loops -DNDEBUG -fopenmp \
       -o V1.5-mpi-lean V1.5-mpi-lean.cpp

if [ $? -eq 0 ]; then
    echo "✅ Build SUCCESS: V1.5-mpi-lean"
    ls -lh V1.5-mpi-lean
else
    echo "❌ Build FAILED"
    exit 1
fi
