#!/bin/bash
#SBATCH -J collatz-mpiomp
#SBATCH -p gpp
#SBATCH -N 4
#SBATCH --ntasks-per-node=56
#SBATCH --cpus-per-task=2
#SBATCH --hint=nomultithread
#SBATCH --exclusive
#SBATCH -t 00:10:00
#SBATCH -o %x-%j.out
#SBATCH -e %x-%j.err
#SBATCH --qos=gp_debug
#SBATCH --account=nct_352

module purge
module load intel impi/2021.10.0   # same stack as you used to build

# Intel MPI pinning: give each rank the whole cpuset and bind OMP inside it
export I_MPI_PIN=1
export I_MPI_PIN_RESPECT_CPUSET=1
export I_MPI_PIN_DOMAIN=omp        # good when using OpenMP threads per rank
export I_MPI_PIN_ORDER=compact     # try 'scatter' if you prefer
# (Enable once to see binding details; then comment out)
# export I_MPI_DEBUG=4
export SLURM_CPU_BIND=verbose
# OpenMP pinning inside the rank
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OMP_PLACES=cores
export OMP_PROC_BIND=close         # try 'spread' as an A/B

# Show the binding Slurm will apply (diagnostic; keep for 1 run)
export SLURM_CPU_BIND=verbose

# Optional: see which libs we link to (helps catch runtime mismatches)
# ldd ./V1.5-mpi-lean | egrep 'iomp|mpi' || true

START=0
COUNT=$((1<<30))
TAG="v15mpi.${SLURM_JOB_ID}"

# Bind rank processes to the cores Slurm assigned
srun --cpu-bind=cores ./V1.5-mpi-lean ${START} ${COUNT} \
     --threads ${OMP_NUM_THREADS} --tag ${TAG}
