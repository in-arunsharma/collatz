#!/bin/bash
# Extract key metrics from perf output for PERFORMANCE.md table

PERF_FILE=${1:-"perf_v1.0.txt"}

if [ ! -f "$PERF_FILE" ]; then
    echo "Usage: ./extract_perf.sh perf_vX.X.txt"
    exit 1
fi

echo "Extracting metrics from $PERF_FILE..."
echo ""

# Extract numbers/sec from program output
NUMSEC=$(grep "Numbers per second:" "$PERF_FILE" | awk '{print $4}')
TIME=$(grep "Total time:" "$PERF_FILE" | awk '{print $3}')

# Extract perf stats (need to handle both atom and core counters)
INSTRUCTIONS=$(grep "instructions" "$PERF_FILE" | grep "cpu_core" | awk '{gsub(/,/, ""); print $1}')
CYCLES=$(grep "cycles" "$PERF_FILE" | grep "cpu_core" | awk '{gsub(/,/, ""); print $1}')
IPC=$(grep "insn per cycle" "$PERF_FILE" | grep "cpu_core" | awk '{print $4}')
CACHE_MISS=$(grep "cache-misses" "$PERF_FILE" | grep "cpu_core" | awk '{gsub(/,/, ""); print $1}')
BRANCH_MISS=$(grep "branch-misses" "$PERF_FILE" | grep "cpu_core" | awk '{gsub(/,/, ""); print $1}')
L1D_MISS=$(grep "L1-dcache-load-misses" "$PERF_FILE" | grep "cpu_core" | awk '{gsub(/,/, ""); print $1}')
L1I_MISS=$(grep "L1-icache-load-misses" "$PERF_FILE" | grep "cpu_core" | awk '{gsub(/,/, ""); print $1}')

# Convert to billions/millions for readability
INS_B=$(echo "scale=1; $INSTRUCTIONS / 1000000000" | bc)
CYC_B=$(echo "scale=1; $CYCLES / 1000000000" | bc)
CM_M=$(echo "scale=1; $CACHE_MISS / 1000000" | bc)
BM_M=$(echo "scale=1; $BRANCH_MISS / 1000000" | bc)
L1D_K=$(echo "scale=0; $L1D_MISS / 1000" | bc)
L1I_M=$(echo "scale=1; $L1I_MISS / 1000000" | bc)

echo "═══════════════════════════════════════════════════════"
echo "PASTE THIS INTO PERFORMANCE.md TABLE:"
echo "═══════════════════════════════════════════════════════"
echo ""
echo "Numbers/sec:  $NUMSEC"
echo "Time (ms):    $TIME"
echo "Instructions: ${INS_B}B"
echo "Cycles:       ${CYC_B}B"
echo "IPC:          $IPC"
echo "Cache Miss:   ${CM_M}M"
echo "Branch Miss:  ${BM_M}M"
echo "L1D Miss:     ${L1D_K}K"
echo "L1I Miss:     ${L1I_M}M"
echo ""
echo "═══════════════════════════════════════════════════════"
