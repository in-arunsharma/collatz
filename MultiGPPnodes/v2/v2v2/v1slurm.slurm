#!/bin/bash
#SBATCH -J collatz71
#SBATCH -p gpp
#SBATCH -N 30
#SBATCH --ntasks-per-node=56
#SBATCH --cpus-per-task=2
#SBATCH --hint=nomultithread
#SBATCH --exclusive
#SBATCH -t 02:00:00
#SBATCH --array=0-0%1
#SBATCH -o collatz_%A_%a.out
#SBATCH -e collatz_%A_%a.err
#SBATCH --qos=gp_debug
#SBATCH --account=nct_352

set -euo pipefail

module purge
module load intel impi/2021.10.0

# Intel MPI / OpenMP pinning
export I_MPI_PIN=1
export I_MPI_PIN_RESPECT_CPUSET=1
export I_MPI_PIN_DOMAIN=omp
export I_MPI_PIN_ORDER=compact

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OMP_PLACES=cores
export OMP_PROC_BIND=close

# -------- Range control --------
COUNT=$((1<<36))                 # raw range length per array task
: "${ARRAY_BASE_INDEX:=0}"       # bump for later waves
TASK_INDEX=$((ARRAY_BASE_INDEX + SLURM_ARRAY_TASK_ID))

# program adds 2^71 internally; change START_BASE to resume
START_BASE=0
START=$((START_BASE + TASK_INDEX * COUNT))
TAG="v15mpi.${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
APP=./V1.5-mpi-lean

# -------- Single srun: prints like before, captures perf to files --------
srun --cpu-bind=cores bash -lc '
  APP="'"$APP"'"; START="'"$START"'"; COUNT="'"$COUNT"'"
  THREADS="'"$OMP_NUM_THREADS"'"; TAG="'"$TAG"'"

  # make program output flush promptly
  export STDOUT_LINEBUF=1
  command -v stdbuf >/dev/null 2>&1 && PREFIX="stdbuf -oL -eL" || PREFIX=""

  if command -v perf >/dev/null 2>&1; then
    # one per node: supported events list (helps interpret later)
    if [ "$SLURM_LOCALID" -eq 0 ]; then
      perf list > perf.$(hostname -s).events.txt 2>&1 || true
    fi

    # perf stat: send perf’s own output to file via stderr redirection
    # child program prints to stdout/stderr as usual (goes to SLURM .out/.err)
    exec $PREFIX perf stat -ddd -x, --wall --no-big-num -- \
      "$APP" "$START" "$COUNT" --threads "$THREADS" --tag "$TAG" \
      2> "perf.$(hostname -s).rank${SLURM_PROCID}.csv"
  else
    exec "$APP" "$START" "$COUNT" --threads "$THREADS" --tag "$TAG"
  fi
'

echo "[OK] Job complete: ${TAG}"
