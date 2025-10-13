#!/bin/bash
# Build script for V1.5-cuda-hybrid (CUDA + OpenMP)

set -e

echo "=== Building V1.5-cuda-hybrid (CUDA + OpenMP) ==="

# Check for nvcc
if ! command -v nvcc &> /dev/null; then
    echo "ERROR: nvcc not found. Please install CUDA toolkit."
    exit 1
fi

# Get git hash
GIT_HASH=$(git rev-parse --short HEAD 2>/dev/null || echo "unknown")
echo "Git commit: $GIT_HASH"

# Detect GPU architecture (SM version)
# RTX 3060 Laptop = sm_86 (Ampere)
# For auto-detection, use: -gencode arch=compute_86,code=sm_86
SM_ARCH=${SM_ARCH:-sm_86}

echo "Target architecture: $SM_ARCH"
echo "Compiler: nvcc"

# CUDA compile flags
NVCC_FLAGS="-O3 -std=c++17 -Xcompiler -fopenmp -Xptxas -O3 --generate-line-info"
NVCC_FLAGS="$NVCC_FLAGS -arch=$SM_ARCH -DGIT_HASH='\"$GIT_HASH\"'"

echo "Flags: $NVCC_FLAGS"

# Build
nvcc $NVCC_FLAGS V1.5-cuda-hybrid.cu -o V1.5-cuda-hybrid

echo "✅ Build successful: V1.5-cuda-hybrid"
echo ""
echo "Usage:"
echo "  ./V1.5-cuda-hybrid 0 1000000 --tag cuda_test"
echo ""
echo "Note: This is a starter template. Full implementation pending."
