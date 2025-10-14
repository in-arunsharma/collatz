#!/bin/bash
# Build script for MNv2 Phase 1 (GPP nodes - MPI + OpenMP)

set -e  # Exit on error

echo "========================================="
echo "MNv2: Building Phase 1 (GPP Nodes)"
echo "========================================="
echo ""

# Detect environment
if command -v module &> /dev/null; then
    echo "[INFO] Module system detected (MareNostrum 5)"
    echo "[INFO] Loading Intel OneAPI and OpenMPI modules..."
    module purge
    module load intel/2023.2 2>/dev/null || module load intel/oneapi 2>/dev/null || echo "[WARN] Intel module not found"
    module load impi/2021.10 2>/dev/null || module load openmpi 2>/dev/null || echo "[WARN] MPI module not found"
    echo ""
else
    echo "[INFO] No module system (local build)"
    echo ""
fi

# Check compiler
if ! command -v mpicxx &> /dev/null; then
    echo "[ERROR] mpicxx not found!"
    echo "[INFO] Install OpenMPI: sudo apt install openmpi-bin libopenmpi-dev"
    exit 1
fi

echo "Compiler: $(which mpicxx)"
mpicxx --version | head -n1
echo ""

# Extract git commit (if available)
GIT_HASH="unknown"
if command -v git &> /dev/null && git rev-parse --git-dir > /dev/null 2>&1; then
    GIT_HASH=$(git rev-parse --short HEAD)
    echo "Git commit: $GIT_HASH"
    echo ""
fi

# Compiler flags
CXX=mpicxx
TARGET=collatz_mpi_gpp
SOURCE=collatz_mpi_gpp.cpp

# Optimization flags (aggressive but safe)
CXXFLAGS=(
    -std=c++17
    -O3
    -march=native
    -fopenmp                    # OpenMP support
    -flto                       # Link-time optimization
    -ffast-math                 # Fast math (safe for Collatz)
    -funroll-loops
    -fno-asynchronous-unwind-tables  # Smaller code
    -DNDEBUG                    # Disable assertions
    -DGIT_COMMIT=\"$GIT_HASH\"  # Embed git hash
    -Wall -Wextra               # Warnings
)

echo "Building: $TARGET"
echo "Source:   $SOURCE"
echo "Flags:    ${CXXFLAGS[*]}"
echo ""

# Compile
echo "[COMPILE] $CXX ${CXXFLAGS[*]} $SOURCE -o $TARGET"
$CXX "${CXXFLAGS[@]}" $SOURCE -o $TARGET

if [ $? -eq 0 ]; then
    echo ""
    echo "========================================="
    echo "✅ BUILD SUCCESSFUL"
    echo "========================================="
    echo "Executable: ./$TARGET"
    ls -lh $TARGET
    echo ""
    echo "Test locally:"
    echo "  mpirun -np 4 ./$TARGET 0 1000000 test"
    echo ""
    echo "Submit to MareNostrum:"
    echo "  sbatch slurm_gpp.slurm"
    echo ""
else
    echo ""
    echo "========================================="
    echo "❌ BUILD FAILED"
    echo "========================================="
    exit 1
fi
