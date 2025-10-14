#!/bin/bash
#SBATCH --job-name=v15mpi_2n
#SBATCH --output=v15mpi_2nodes_%j.out
#SBATCH --error=v15mpi_2nodes_%j.err
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=1      # 1 MPI rank per node (KEY!)
#SBATCH --cpus-per-task=112      # Full GPP node
#SBATCH --hint=nomultithread     # Disable SMT
#SBATCH --time=00:10:00
#SBATCH --qos=gp_debug
#SBATCH --account=nct_352

echo "=== V1.5-mpi-lean: GCC Build (icpc killed performance!) ==="
echo "Date: $(date)"
echo "Nodes: $SLURM_NODELIST"
echo "Job: $SLURM_JOB_ID"
echo ""

# Use GCC with OpenMPI (NOT Intel Classic Compiler!)
module purge
module load gcc/13.2.0 openmpi/4.1.5-gcc

# For Intel MPI with GCC backend, use instead:
# module purge
# module load gcc/13.2.0 impi/2021.10.0
# export I_MPI_CXX=g++
# export I_MPI_PIN=1
# export I_MPI_PIN_DOMAIN=omp

# OpenMP settings
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PLACES=cores
export OMP_PROC_BIND=spread

echo "Configuration:"
echo "  Compiler: GCC (via I_MPI_CXX=g++)"
echo "  MPI ranks: $SLURM_NTASKS (1 per node)"
echo "  CPUs/task: $SLURM_CPUS_PER_TASK"
echo "  OMP_NUM_THREADS: $OMP_NUM_THREADS"
echo "  OMP_PROC_BIND: $OMP_PROC_BIND"
echo "  OMP_PLACES: $OMP_PLACES"
echo "  I_MPI_PIN_DOMAIN: $I_MPI_PIN_DOMAIN"
echo ""

# Use srun with CPU binding (friend's advice)
srun --cpu-bind=cores ./V1.5-mpi-lean 0 100000000 --threads $OMP_NUM_THREADS --cold off

echo ""
echo "Expected: ~274M nums/sec (2 × 137M)"
echo "If per-rank is not ~137M, check binding!"
