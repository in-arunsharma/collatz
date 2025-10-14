#!/bin/bash
# Quick diagnostic for the performance issue

echo "=== Checking latest job output ==="
ls -lt v15mpi_1node_sanity_*.err | head -1
echo ""

echo "=== Checking for thread diagnostics in stderr ==="
cat v15mpi_1node_sanity_*.err 2>/dev/null | grep -i "thread\|rank\|omp" || echo "No stderr output found or no thread messages"

echo ""
echo "=== Checking binary size and OpenMP linkage ==="
ls -lh V1.5-mpi-lean
ldd V1.5-mpi-lean | grep -i "omp\|gomp" || echo "⚠️  NO OpenMP library linked!"

echo ""
echo "=== Checking if OpenMP symbols exist ==="
nm V1.5-mpi-lean 2>/dev/null | grep -i "omp_get_num_threads" || echo "⚠️  NO OpenMP symbols found!"

echo ""
echo "=== Diagnosis ==="
echo "If 'NO OpenMP library linked' appears above:"
echo "  → OpenMP wasn't linked! The program ran single-threaded."
echo "  → Solution: Check compiler flags, try -fopenmp vs -qopenmp"
echo ""
echo "If stderr shows 'OpenMP started with 1 thread':"
echo "  → Environment variable issue (OMP_NUM_THREADS not set)"
echo ""
echo "If stderr shows 'OpenMP started with 112 threads':"
echo "  → Threads are running but something else is slow (investigate further)"
