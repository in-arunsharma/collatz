#!/bin/bash
#SBATCH -J collatz71
#SBATCH -p gpp
#SBATCH -N 30
#SBATCH --ntasks-per-node=56
#SBATCH --cpus-per-task=2
#SBATCH --hint=nomultithread
#SBATCH --exclusive
#SBATCH -t 02:00:00
#SBATCH --array=0-249%1
#SBATCH -o collatz_%A_%a.out
#SBATCH -e collatz_%A_%a.err
#SBATCH --qos=gp_debug
#SBATCH --account=nct_352

module purge
module load intel impi/2021.10.0   # same stack as you used to build

# Intel MPI pinning
export I_MPI_PIN=1
export I_MPI_PIN_RESPECT_CPUSET=1
export I_MPI_PIN_DOMAIN=omp
export I_MPI_PIN_ORDER=compact
export SLURM_CPU_BIND=verbose

# OpenMP pinning inside the rank
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OMP_PLACES=cores
export OMP_PROC_BIND=close

# Per-task range control
COUNT=$((1<<42))                    # numbers per array task (raw range length)
: "${ARRAY_BASE_INDEX:=0}"          # set to 250, 500, ... on later waves
TASK_INDEX=$((ARRAY_BASE_INDEX + SLURM_ARRAY_TASK_ID))
START=$((TASK_INDEX * COUNT))
TAG="v15mpi.${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"

# Run (program already adds base 2^71 internally)
srun --cpu-bind=cores ./V1.5-mpi-lean "${START}" "${COUNT}" \
     --threads "${OMP_NUM_THREADS}" --tag "${TAG}"
