#!/bin/bash
#SBATCH --job-name=v17_2nodes
#SBATCH --output=collatz_v17_2nodes_%j.out
#SBATCH --error=collatz_v17_2nodes_%j.err
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

echo "=== V1.7 Two-Node Scaling Test ==="
echo "Date: $(date)"
echo "Nodes: 2"
echo ""

mpirun -np 2 ./V1.7-mpi 0 100000000 --threads 112

echo ""
echo "Expected: ~274M nums/sec (2× single-node)"
