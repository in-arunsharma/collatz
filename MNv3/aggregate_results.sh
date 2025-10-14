#!/bin/bash
# Aggregate results from V1.5 array job

if [ $# -ne 1 ]; then
    echo "Usage: $0 <job_id>"
    echo "Example: $0 30669600"
    exit 1
fi

JOB_ID=$1

echo "=== Aggregating Results for Job $JOB_ID ==="
echo ""

# Find all output files for this job
FILES=(v15_array*_${JOB_ID}_*.out)

if [ ${#FILES[@]} -eq 0 ]; then
    echo "ERROR: No output files found for job $JOB_ID"
    exit 1
fi

TOTAL_TESTED=0
TOTAL_TIME_MS=0
TOTAL_STEPS=0
MAX_STEPS=0
NODES=0

echo "Processing ${#FILES[@]} tasks..."
echo ""

for file in "${FILES[@]}"; do
    if [ -f "$file" ]; then
        NODES=$((NODES + 1))
        
        # Extract statistics using grep/awk
        TESTED=$(grep "^Tested:" "$file" | awk '{print $2}')
        TIME_MS=$(grep "^Time:" "$file" | awk '{print $2}')
        STEPS=$(grep "^Avg steps:" "$file" | awk '{print $3}')
        MAX=$(grep "^Max steps:" "$file" | awk '{print $3}')
        
        if [ -n "$TESTED" ]; then
            TOTAL_TESTED=$((TOTAL_TESTED + TESTED))
            TOTAL_TIME_MS=$((TOTAL_TIME_MS > TIME_MS ? TOTAL_TIME_MS : TIME_MS))
            TOTAL_STEPS=$((TOTAL_STEPS + TESTED * STEPS / 1))
            MAX_STEPS=$((MAX_STEPS > MAX ? MAX_STEPS : MAX))
            
            echo "Task $(basename $file): $TESTED numbers in ${TIME_MS}ms"
        fi
    fi
done

echo ""
echo "=== AGGREGATE RESULTS ==="
echo "Nodes:        $NODES"
echo "Total tested: $TOTAL_TESTED numbers"
echo "Wall time:    ${TOTAL_TIME_MS} ms (max across all tasks)"

if [ $TOTAL_TIME_MS -gt 0 ]; then
    THROUGHPUT=$((TOTAL_TESTED * 1000 / TOTAL_TIME_MS))
    SPEEDUP=$(echo "scale=2; $THROUGHPUT / 2200000" | bc)
    
    echo "Throughput:   ${THROUGHPUT} nums/sec"
    echo "Per-node:     $((THROUGHPUT / NODES)) nums/sec"
    echo "Speedup:      ${SPEEDUP}× (vs 2.2M/s sequential)"
    echo ""
    
    if [ $NODES -eq 2 ]; then
        echo "Expected: ~274M nums/sec (2 × 137M)"
        EFFICIENCY=$(echo "scale=1; $THROUGHPUT / 274000000 * 100" | bc)
        echo "Scaling efficiency: ${EFFICIENCY}%"
    elif [ $NODES -eq 10 ]; then
        echo "Expected: ~1.37B nums/sec (10 × 137M)"
        EFFICIENCY=$(echo "scale=1; $THROUGHPUT / 1370000000 * 100" | bc)
        echo "Scaling efficiency: ${EFFICIENCY}%"
    fi
fi

echo ""
echo "=== Files ==="
ls -lh v15_array*_${JOB_ID}_*
