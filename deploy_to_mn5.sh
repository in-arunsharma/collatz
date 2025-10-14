#!/bin/bash
# Quick deployment script for MareNostrum 5

set -e

echo "========================================="
echo "MNv2 Deployment to MareNostrum 5"
echo "========================================="
echo ""

# Check if archive exists
if [ ! -f "mnv2.tar.gz" ]; then
    echo "[ERROR] mnv2.tar.gz not found!"
    echo "[INFO] Run: tar czf mnv2.tar.gz MNv2/"
    exit 1
fi

echo "Step 1: Transfer archive to MareNostrum..."
echo "Running: scp mnv2.tar.gz nct01225@glogin1.bsc.es:~/"
echo ""
scp mnv2.tar.gz nct01225@glogin1.bsc.es:~/

echo ""
echo "✅ Transfer complete!"
echo ""
echo "========================================="
echo "Next Steps (run on MareNostrum):"
echo "========================================="
echo ""
echo "ssh nct01225@glogin1.bsc.es"
echo ""
echo "# Extract in projects directory"
echo "cd /gpfs/projects/nct_352/nct01225/collatz/"
echo "tar xzf ~/mnv2.tar.gz"
echo "cd MNv2/"
echo ""
echo "# Create output directory"
echo "mkdir -p /gpfs/scratch/nct_352/nct01225/collatz_output"
echo ""
echo "# Build"
echo "./build_gpp.sh"
echo ""
echo "# Quick test (login node, small workload)"
echo "mpirun -np 2 ./collatz_mpi_gpp 0 1000000 quick_test"
echo ""
echo "# Submit to cluster (2 GPP nodes, 1B seeds)"
echo "sbatch slurm_gpp.slurm"
echo ""
echo "# Monitor job"
echo "squeue -u \$USER"
echo "tail -f /gpfs/scratch/nct_352/nct01225/collatz_output/mnv2_gpp_*.out"
echo ""
echo "========================================="
echo "Expected Performance:"
echo "  2 GPP nodes = ~500M nums/sec"
echo "  1B seeds in ~2 seconds"
echo "========================================="
