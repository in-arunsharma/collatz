#!/bin/bash
#SBATCH -J collatz71
#SBATCH -p acc
#SBATCH -N 6
#SBATCH --ntasks-per-node=4         # 1 rank per GPU
#SBATCH --gpus-per-node=4
#SBATCH --cpus-per-task=8
#SBATCH --hint=nomultithread
#SBATCH --exclusive
#SBATCH -t 02:00:00
#SBATCH --array=0-1%1
#SBATCH -o collatz_%A_%a.out
#SBATCH -e collatz_%A_%a.err
#SBATCH --qos=acc_debug
#SBATCH --account=nct_352

# --- Modules: ACC load order matters ---
module purge
module load EB/install
module load intel
module load CUDA/12.6.0            # pick the exact version you built with
module load impi/2021.10.0         # same MPI family you built against

# --- Per-task range ---
COUNT=$((1<<30))
: "${ARRAY_BASE_INDEX:=0}"
TASK_INDEX=$((ARRAY_BASE_INDEX + SLURM_ARRAY_TASK_ID))
START=$((TASK_INDEX * COUNT))
TAG="v15cuda.${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"

# --- Bindings: 1 GPU per task, compact CPU cores ---
export CUDA_DEVICE_ORDER=PCI_BUS_ID
export OMP_NUM_THREADS=1
export SLURM_CPU_BIND=cores

# Intel MPI pinning (safe defaults)
export I_MPI_PIN=1
export I_MPI_PIN_RESPECT_CPUSET=1
export I_MPI_PIN_DOMAIN=core
export I_MPI_PIN_ORDER=compact

# If your site prefers PMI2 with Intel MPI, swap --mpi=pmix -> --mpi=pmi2
srun --mpi=pmix \
     --gpu-bind=single:1 \
     --cpu-bind=cores \
     ./collatz_gpu_mpi_cycles "${START}" "${COUNT}" \
        --small-limit 20 \
        --tag "${TAG}"
        # add --no-print-cycles to suppress per-cycle lines (summary still prints)
