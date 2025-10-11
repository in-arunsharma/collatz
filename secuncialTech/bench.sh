#!/bin/bash
# Quick benchmark script - runs 1M numbers with perf stats

VERSION=${1:-"1.0"}
COUNT=${2:-1000000}

echo "=== Benchmarking V${VERSION} with ${COUNT} numbers ==="
echo ""

cd /home/aruns/Desktop/MN25/secuncialTech

if [ ! -f "V${VERSION}" ]; then
    echo "Error: V${VERSION} not found. Build it first with ./build.sh ${VERSION}"
    exit 1
fi

# Run with perf stat
sudo perf stat -e instructions,cycles,cache-misses,branch-misses,L1-dcache-load-misses,L1-icache-load-misses \
    ./V${VERSION} 0 ${COUNT} 2>&1 | tee perf_v${VERSION}.txt

echo ""
echo "Results saved to perf_v${VERSION}.txt"
echo ""
echo "Key metrics to add to PERFORMANCE.md:"
echo "- Instructions (total)"
echo "- Cycles (total)"  
echo "- IPC (insn per cycle)"
echo "- Cache misses"
echo "- Branch misses"
