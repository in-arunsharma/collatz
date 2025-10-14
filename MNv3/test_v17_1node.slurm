#!/bin/bash
#SBATCH --job-name=v17_1node
#SBATCH --output=collatz_v17_1node_%j.out
#SBATCH --error=collatz_v17_1node_%j.err
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

echo "=== V1.7 Single-Node Test (Based on proven V1.5) ==="
echo "Date: $(date)"
echo "Node: $(hostname)"
echo ""

mpirun -np 1 ./V1.7-mpi 0 100000000 --threads 112

echo ""
echo "Expected: 137M nums/sec (V1.5 baseline)"
