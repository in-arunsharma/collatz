#!/bin/bash
# Quick local test before deploying to MareNostrum
# Tests with 4 threads on a small range

set -e

echo "========================================="
echo "Local Test - V1.5-openmp"
echo "========================================="
echo ""

# Build if needed
if [ ! -f "V1.5-openmp" ]; then
    echo "Building V1.5-openmp..."
    bash build_openmp.sh
    echo ""
fi

# Test parameters
THREADS=4
OFFSET=0
COUNT=100000  # 100K numbers (should take ~1 second with 4 threads)
TAG="local_test"

echo "Running local test:"
echo "  Threads: $THREADS"
echo "  Range: [$OFFSET, $((OFFSET + COUNT)))"
echo "  Expected time: ~1 second"
echo ""

# Run test
./V1.5-openmp $OFFSET $COUNT --threads $THREADS --small-limit 20 --tag $TAG

echo ""
echo "========================================="
echo "Test complete!"
echo "========================================="
echo ""
echo "If you see:"
echo "  ✓ Throughput >100K nums/sec → Build is correct"
echo "  ✓ Speedup >3× → OpenMP is working"
echo "  ✓ No errors → Ready for MareNostrum"
echo ""
echo "Next steps:"
echo "  1. Review output above"
echo "  2. If successful, deploy to MareNostrum:"
echo "     scp -r MNv3/ nct01225@glogin1.bsc.es:/gpfs/projects/nct_352/nct01225/collatz/"
echo "  3. SSH to MareNostrum and submit job:"
echo "     sbatch MNv3/slurm_openmp_1node.slurm"
