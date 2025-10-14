#!/bin/bash
#SBATCH -J v15mpi_2nodes
#SBATCH -p gpp
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=112
#SBATCH --hint=nomultithread
#SBATCH -t 00:05:00
#SBATCH -o v15mpi_%j.out
#SBATCH -e v15mpi_%j.err

module purge
module load gcc/12 impi/2021.10.0
export I_MPI_CXX=g++

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PLACES=cores
export OMP_PROC_BIND=spread
export I_MPI_PIN=1
export I_MPI_PIN_DOMAIN=omp
export KMP_BLOCKTIME=0

srun --cpu-bind=cores ./V1.5-mpi-lean 0 1000000000 --threads $OMP_NUM_THREADS
