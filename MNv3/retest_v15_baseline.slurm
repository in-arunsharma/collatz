#!/bin/bash
#SBATCH --job-name=v15_retest
#SBATCH --output=v15_retest_%j.out
#SBATCH --error=v15_retest_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=112
#SBATCH --time=00:10:00
#SBATCH --qos=gp_debug
#SBATCH --account=nct_352

echo "=== V1.5-openmp Baseline Retest ==="
echo "Date: $(date)"
echo "Node: $(hostname)"
echo "Job: $SLURM_JOB_ID"
echo ""

# OpenMP configuration (same as original successful run)
export OMP_NUM_THREADS=112
export OMP_PROC_BIND=spread
export OMP_PLACES=cores

echo "Configuration:"
echo "  Threads: $OMP_NUM_THREADS"
echo "  Binding: $OMP_PROC_BIND"
echo "  Places:  $OMP_PLACES"
echo ""

# Same test as original: 100M range, 2^20 memo
echo "Running: ./V1.5-openmp 0 100000000 --threads 112 --small-limit 20"
echo ""

time ./V1.5-openmp 0 100000000 --threads 112 --small-limit 20

echo ""
echo "=== Expected: 137M nums/sec (from job 30668700) ==="
