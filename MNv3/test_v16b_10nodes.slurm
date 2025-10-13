#!/bin/bash
#SBATCH --job-name=v16b_10nodes
#SBATCH --output=collatz_v16b_10nodes_%j.out
#SBATCH --error=collatz_v16b_10nodes_%j.err
#SBATCH --nodes=10
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=112
#SBATCH --time=00:15:00
#SBATCH --qos=gp_bsccs
#SBATCH --account=nct_352

module load intel impi

export OMP_NUM_THREADS=112
export OMP_PROC_BIND=spread
export OMP_PLACES=cores

echo "=== V1.6b Ten-Node Scaling Test ==="
echo "Date: $(date)"
echo "Nodes: 10"
echo "Total Cores: 1,120"
echo "Purpose: Maximum GPP throughput for Phase 2"
echo ""

mpirun -np 10 ./V1.6b-mpi-openmp 0 100000000 --tag v16b_10nodes

echo ""
echo "=== Expected Results ==="
echo "Throughput: 1.3-1.6B nums/sec"
echo "Scaling: 8-10× vs single-node"
echo "Efficiency: 80-100%"
echo "Per-node: ~150M nums/sec"
echo ""
echo "✅ If throughput > 1.2B: EXCELLENT - Phase 2 complete!"
echo "⚠️  If throughput 800M-1.2B: GOOD - Some scaling loss"
echo "❌ If throughput < 800M: FAIL - Scaling bottleneck"
