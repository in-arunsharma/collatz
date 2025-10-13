#!/bin/bash
# build_phase1_gpp.sh - Build Phase 1: GPP Nodes Only (MPI + OpenMP)

set -e  # Exit on error

echo "================================================================"
echo "  Building Phase 1: GPP Nodes (MPI + OpenMP)"
echo "================================================================"
echo ""

# Detect environment
if [[ -n "$SLURM_JOB_ID" ]]; then
    echo "Running on MareNostrum 5"
    ON_MARENOSTRUM=true
else
    echo "Running locally (testing)"
    ON_MARENOSTRUM=false
fi

# Load modules (MareNostrum 5 specific)
if [[ "$ON_MARENOSTRUM" == true ]]; then
    echo "Loading MareNostrum 5 modules..."
    module purge
    module load intel/oneapi
    module load openmpi
    module list
    echo ""
fi

# Compiler settings
if [[ "$ON_MARENOSTRUM" == true ]]; then
    CXX=mpic++
    CXXFLAGS="-O3 -march=native -std=c++17 -fopenmp"
else
    CXX=${MPICXX:-mpic++}
    CXXFLAGS="-O3 -march=native -std=c++17 -fopenmp -Wno-psabi"
fi

# Git hash for reproducibility
GIT_HASH=$(git rev-parse --short HEAD 2>/dev/null || echo "unknown")
CXXFLAGS="$CXXFLAGS -DGIT_HASH='\"$GIT_HASH\"'"

# Additional optimization flags
CXXFLAGS="$CXXFLAGS -flto -fno-exceptions -fno-rtti -funroll-loops"
CXXFLAGS="$CXXFLAGS -fno-asynchronous-unwind-tables -DNDEBUG"

echo "Compiler: $CXX"
echo "Flags: $CXXFLAGS"
echo "Git commit: $GIT_HASH"
echo ""

# Build
echo "Compiling collatz_mpi_gpp.cpp..."
$CXX $CXXFLAGS -I. collatz_mpi_gpp.cpp -o collatz_mpi_gpp

if [[ $? -eq 0 ]]; then
    echo ""
    echo "✅ Build successful: collatz_mpi_gpp"
    echo ""
    echo "Executable info:"
    ls -lh collatz_mpi_gpp
    echo ""
    echo "================================================================"
    echo "  Build complete!"
    echo "================================================================"
    echo ""
    echo "Next steps:"
    echo "  1. Test locally:  mpirun -np 4 ./collatz_mpi_gpp 0 1000000 test"
    echo "  2. Deploy to MN5: sbatch slurm_phase1_gpp.slurm"
    echo ""
else
    echo ""
    echo "❌ Build failed!"
    exit 1
fi
