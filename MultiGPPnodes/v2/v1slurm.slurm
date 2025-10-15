#!/bin/bash
#SBATCH -J collatz-array
#SBATCH -p gpp
#SBATCH -N 30
#SBATCH --ntasks-per-node=56
#SBATCH --cpus-per-task=2
#SBATCH --hint=nomultithread
#SBATCH --exclusive
#SBATCH -t 02:00:00
#SBATCH -o %x-%A_%a.out
#SBATCH -e %x-%A_%a.err
#SBATCH --qos=gp_debug
#SBATCH --account=nct_352

#SBATCH --array=0-350%1

module purge
module load intel impi/2021.10.0   # same stack as you used to build

# Intel MPI pinning: give each rank the whole cpuset and bind OMP inside it
export I_MPI_PIN=1
export I_MPI_PIN_RESPECT_CPUSET=1
export I_MPI_PIN_DOMAIN=omp        # good when using OpenMP threads per rank
export I_MPI_PIN_ORDER=compact     # try 'scatter' if you prefer

export I_MPI_DEBUG=0
unset COLLATZ_DEBUG

# (Enable once to see binding details; then comment out)
# export I_MPI_DEBUG=4
# OpenMP pinning inside the rank
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OMP_PLACES=cores
export OMP_PROC_BIND=close         # try 'spread' as an A/B
export OMP_DYNAMIC=false

# Show the binding Slurm will apply (diagnostic; keep for 1 run)
unset SLURM_CPU_BIND

# Optional: see which libs we link to (helps catch runtime mismatches)
# ldd ./V1.5-mpi-lean | egrep 'iomp|mpi' || true

CHUNK=$((1<<40))                #
START=$(( SLURM_ARRAY_TASK_ID * CHUNK ))
COUNT=${CHUNK}
TAG="v15mpi.${SLURM_ARRAY_JOB_ID}.a${SLURM_ARRAY_TASK_ID}"
echo "TASK $SLURM_ARRAY_TASK_ID: START=$START COUNT=$COUNT"

# Bind rank processes to the cores Slurm assigned
srun --kill-on-bad-exit=1 --cpu-bind=cores ./V1.5-mpi-lean ${START} ${COUNT} \
     --threads ${OMP_NUM_THREADS} --tag ${TAG} --no-logs

