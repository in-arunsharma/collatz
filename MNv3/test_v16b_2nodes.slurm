#!/bin/bash
#SBATCH --job-name=v16b_2nodes
#SBATCH --output=collatz_v16b_2nodes_%j.out
#SBATCH --error=collatz_v16b_2nodes_%j.err
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=112
#SBATCH --time=00:10:00
#SBATCH --qos=gp_debug
#SBATCH --account=nct_352

module load intel impi

export OMP_NUM_THREADS=112
export OMP_PROC_BIND=spread
export OMP_PLACES=cores

echo "=== V1.6b Two-Node Scaling Test ==="
echo "Date: $(date)"
echo "Nodes: 2"
echo "Total Cores: 224"
echo "Purpose: Verify MPI bug fix - should scale 2×, not 0.018× like V1.6!"
echo ""

mpirun -np 2 ./V1.6b-mpi-openmp 0 100000000 --tag v16b_2nodes

echo ""
echo "=== Expected Results ==="
echo "Throughput: 270-330M nums/sec"
echo "Scaling: 1.8-2.0× vs single-node"
echo "Efficiency: 90-100%"
echo ""
echo "V1.6 BROKEN result: 2.5M nums/sec (54× slowdown)"
echo ""
echo "✅ If throughput > 250M: SUCCESS - Bug fix worked!"
echo "❌ If throughput < 100M: FAIL - MPI bug still present"
