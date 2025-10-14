#!/bin/bash
module purge
module load gcc/13.2.0

echo "Building V1.5-openmp-checkpoint with GCC..."

g++ -O3 -march=native -flto -fno-exceptions -fno-rtti -funroll-loops -DNDEBUG -fopenmp \
    -o V1.5-openmp-checkpoint V1.5-openmp-checkpoint.cpp

if [ $? -eq 0 ]; then
    echo "✅ Build SUCCESS: V1.5-openmp-checkpoint"
    ls -lh V1.5-openmp-checkpoint
else
    echo "❌ Build FAILED"
    exit 1
fi
