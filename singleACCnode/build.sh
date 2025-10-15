#!/usr/bin/env bash
set -euo pipefail

# ---- Modules: fixed order for ACC ----
module purge
module load EB/install
module load intel
module load CUDA/12.6.0
module load impi/2021.10.0

# ---- Build ----
SRC="collatz_gpu_mpi_cycles.cu"
OBJ="${SRC%.cu}.o"
OUT="collatz_gpu_mpi_cycles"

# H100 => sm_90  (use: NV_ARCH=sm_80 bash build.sh  if you ever target A100)
NV_ARCH="${NV_ARCH:-sm_90}"

echo "[INFO] nvcc: $(nvcc --version | tail -1)"
echo "[INFO] mpicxx: $(mpicxx --version 2>/dev/null | head -1 || echo found)"
echo "[INFO] arch: ${NV_ARCH}"

# Compile .cu to object
nvcc -O3 -std=c++17 -arch="${NV_ARCH}" -Xcompiler -fPIC -c "${SRC}" -o "${OBJ}"

# Link with MPI wrapper + cudart (use mpicxx, NOT mpiicpc)
CUDA_HOME="${CUDA_HOME:-$(dirname "$(dirname "$(which nvcc)")")}"
mpicxx -Ofast -flto -std=c++17 \
       "${OBJ}" -L"${CUDA_HOME}/lib64" -lcudart \
       -o "${OUT}"

echo "[OK] Build complete: ./${OUT}"
