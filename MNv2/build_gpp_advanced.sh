#!/bin/bash

# Enhanced build script with advanced optimizations for Sapphire Rapids

echo "========================================="
echo "MNv2: Building Phase 1 (Advanced Optimizations)"
echo "========================================="

# Load Intel compiler if available
if command -v module >/dev/null 2>&1; then
    echo "[INFO] Loading Intel OneAPI modules"
    module purge
    module load intel/2023.2
    module load impi/2021.10
else
    echo "[INFO] No module system (local build)"
fi

# Determine compiler
if command -v icpx >/dev/null 2>&1; then
    COMPILER="icpx"
    echo "Compiler: Intel icpx (latest)"
elif command -v icpc >/dev/null 2>&1; then
    COMPILER="icpc"  
    echo "Compiler: Intel icpc (classic)"
elif command -v mpicxx >/dev/null 2>&1; then
    COMPILER="mpicxx"
    echo "Compiler: $(mpicxx --version | head -1)"
else
    echo "ERROR: No suitable compiler found"
    exit 1
fi

# Get git commit
GIT_COMMIT=$(git rev-parse --short HEAD 2>/dev/null || echo "unknown")
echo "Git commit: $GIT_COMMIT"

echo ""
echo "Building: collatz_mpi_gpp_advanced"
echo "Source:   collatz_mpi_gpp.cpp"

# ADVANCED OPTIMIZATION FLAGS for Sapphire Rapids
ADVANCED_FLAGS="-std=c++17 -O3 -DNDEBUG"

if [[ "$COMPILER" == "icpx" || "$COMPILER" == "icpc" ]]; then
    # Intel-specific optimizations (conservative)
    ADVANCED_FLAGS="$ADVANCED_FLAGS -march=native -mtune=native"
    ADVANCED_FLAGS="$ADVANCED_FLAGS -ipo -qopt-prefetch=2"
else
    # GCC-specific optimizations (conservative)
    ADVANCED_FLAGS="$ADVANCED_FLAGS -march=native -mtune=native"  
    ADVANCED_FLAGS="$ADVANCED_FLAGS -funroll-loops"
fi

# Common optimizations
ADVANCED_FLAGS="$ADVANCED_FLAGS -fopenmp -flto -ffast-math -funroll-loops"
ADVANCED_FLAGS="$ADVANCED_FLAGS -fno-asynchronous-unwind-tables -DGIT_COMMIT=\"$GIT_COMMIT\""
ADVANCED_FLAGS="$ADVANCED_FLAGS -Wall -Wextra"

echo "Flags:    $ADVANCED_FLAGS"

# Compile with advanced optimizations
echo "[COMPILE] $COMPILER $ADVANCED_FLAGS collatz_mpi_gpp.cpp -o collatz_mpi_gpp_advanced"

$COMPILER $ADVANCED_FLAGS collatz_mpi_gpp.cpp -o collatz_mpi_gpp_advanced

if [ $? -eq 0 ]; then
    echo ""
    echo "========================================="
    echo "✅ BUILD SUCCESSFUL (ADVANCED)"
    echo "========================================="
    echo "Executable: ./collatz_mpi_gpp_advanced"
    ls -lh collatz_mpi_gpp_advanced
    echo ""
    echo "Expected improvement: 5-15% over standard build"
    echo ""
    echo "Test locally:"
    echo "  mpirun -np 4 ./collatz_mpi_gpp_advanced 0 1000000 test"
    echo ""
    echo "Submit to MareNostrum:"
    echo "  sbatch slurm_gpp_advanced.slurm"
else
    echo ""
    echo "❌ BUILD FAILED"
    echo "Check compiler flags and dependencies"
    exit 1
fi