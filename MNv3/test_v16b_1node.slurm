#!/bin/bash
#SBATCH --job-name=v16b_1node
#SBATCH --output=collatz_v16b_1node_%j.out
#SBATCH --error=collatz_v16b_1node_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=112
#SBATCH --time=00:10:00
#SBATCH --qos=gp_debug
#SBATCH --account=nct_352

module load intel impi

export OMP_NUM_THREADS=112
export OMP_PROC_BIND=spread
export OMP_PLACES=cores

echo "=== V1.6b Single-Node Validation Test ==="
echo "Date: $(date)"
echo "Node: $(hostname)"
echo "Cores: 112"
echo "Purpose: Verify bug fixes match Phase 1 baseline"
echo ""

mpirun -np 1 ./V1.6b-mpi-openmp 0 100000000 --tag v16b_1node

echo ""
echo "=== Expected Results ==="
echo "Throughput: 137-165M nums/sec"
echo "Baseline (V1.5): 137M nums/sec"
echo "Boost from L2 cache (2^19): +10-20%"
echo ""
echo "✅ If throughput > 130M: PASS - Proceed to 2-node test"
echo "❌ If throughput < 130M: FAIL - Debug V1.6b before multi-node"
