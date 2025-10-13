#!/bin/bash
# Quick integration helper - links V1.4b core functions into Phase 1

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
V14B_PATH="../parallel/V1.4b-openmp.cpp"
WORKER_FILE="worker_gpp_node.cpp"

echo "========================================="
echo "V1.4b Integration Helper"
echo "========================================="

# Check if V1.4b exists
if [ ! -f "$SCRIPT_DIR/$V14B_PATH" ]; then
    echo "ERROR: V1.4b-openmp.cpp not found at $V14B_PATH"
    echo "Please ensure parallel/V1.4b-openmp.cpp exists"
    exit 1
fi

echo "✓ Found V1.4b-openmp.cpp"

# Backup worker file
cp "$SCRIPT_DIR/$WORKER_FILE" "$SCRIPT_DIR/${WORKER_FILE}.backup"
echo "✓ Backed up $WORKER_FILE to ${WORKER_FILE}.backup"

# Check if already integrated
if grep -q "#include.*V1.4b-openmp.cpp" "$SCRIPT_DIR/$WORKER_FILE"; then
    echo "⚠ V1.4b already integrated in $WORKER_FILE"
    echo "If you want to re-integrate, restore from backup first:"
    echo "  mv ${WORKER_FILE}.backup $WORKER_FILE"
    exit 0
fi

echo ""
echo "Integrating V1.4b core functions..."
echo ""

# Create temporary file with integration
cat > /tmp/worker_integration.tmp << 'EOF'
// ============================================================================
// COLLATZ CORE - V1.4b-openmp Integration
// ============================================================================
// Include the proven sequential algorithm (V1.4b-openmp.cpp)
// This provides: process_number(), initialize_memoization(), etc.

#include "../parallel/V1.4b-openmp.cpp"

// ============================================================================

EOF

# Insert after the existing includes
sed -i '/^#include.*node_config.hpp/a \
\
// V1.4b-openmp core integration\
#include "../parallel/V1.4b-openmp.cpp"' "$SCRIPT_DIR/$WORKER_FILE"

echo "✓ Added #include directive for V1.4b-openmp.cpp"

# Replace the dummy computation
sed -i 's/uint32_t steps = (n % 100);.*$/uint32_t steps = process_number(n);  \/\/ V1.4b core/' "$SCRIPT_DIR/$WORKER_FILE"

echo "✓ Replaced dummy computation with process_number(n)"

# Add initialization call (if not already there)
if ! grep -q "initialize_memoization()" "$SCRIPT_DIR/$WORKER_FILE"; then
    sed -i '/GPPNodeWorker::GPPNodeWorker/,/^}/ {
        /seeds_per_thread =/ a\
    \
    // Initialize V1.4b memoization tables\
    initialize_memoization();
    }' "$SCRIPT_DIR/$WORKER_FILE"
    echo "✓ Added initialize_memoization() call to constructor"
fi

echo ""
echo "========================================="
echo "Integration Complete!"
echo "========================================="
echo ""
echo "Changes made to $WORKER_FILE:"
echo "  1. Added #include for V1.4b-openmp.cpp"
echo "  2. Replaced dummy computation with process_number(n)"
echo "  3. Added initialize_memoization() in constructor"
echo ""
echo "Backup saved: ${WORKER_FILE}.backup"
echo ""
echo "Next steps:"
echo "  1. Build: ./build_phase1_gpp.sh"
echo "  2. Test:  ./test_local_mpi.sh 2"
echo "  3. If successful, deploy to MareNostrum!"
echo ""
echo "========================================="
