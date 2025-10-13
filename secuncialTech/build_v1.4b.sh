#!/bin/bash
# Build script for V1.4b (Hot/Cold Queue Architecture)

set -e  # Exit on error

echo "Building V1.4b (Hot/Cold Queue with 256-bit & Brent)"
echo "=============================================="

# Get git commit hash
GIT_HASH=$(git rev-parse --short HEAD 2>/dev/null || echo "unknown")
echo "Git commit: $GIT_HASH"

# Compilation flags (same as V1.3d/V1.4 for fair comparison)
FLAGS="-O3 -march=native -std=c++17 -Wall -Wextra -flto -fno-exceptions -fno-rtti -funroll-loops -fno-asynchronous-unwind-tables -DNDEBUG"

echo "Flags: $FLAGS"
echo ""

# Compile
g++ $FLAGS -DGIT_HASH=\"$GIT_HASH\" V1.4b.cpp -o V1.4b

echo "✅ Build complete: V1.4b"
echo ""
echo "Features:"
echo "  - Hot path: Pure 128-bit, 100K fuse (~2M nums/sec)"
echo "  - Cold queue 1: Fuse hits → 1M fuse + Brent cycle detection"
echo "  - Cold queue 2: 128-bit overflow → 256-bit + Brent"
echo "  - Batch processing: Every 10K hot seeds"
echo "  - Integrated: Single program, no post-processing"
echo ""
echo "Run with: ./V1.4b 0 1000000 --tag <your_tag>"
