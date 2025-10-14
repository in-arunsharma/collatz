#!/bin/bash
# Build script for V1.5-mpi-lean.cpp (MPI + OpenMP)
set -e

SRC=V1.5-mpi-lean.cpp
BIN=V1.5-mpi-lean

mpicxx -O3 -march=native -fopenmp -std=c++17 -o $BIN $SRC

echo "Build complete: $BIN"