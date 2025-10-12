#!/bin/bash
# Build script for V1.4 with metadata embedding

# Get git hash if available
GIT_HASH=$(git rev-parse --short HEAD 2>/dev/null || echo "unknown")

# Compile flags
FLAGS="-O3 -march=native -std=c++17 -Wall -Wextra -flto -fno-exceptions -fno-rtti -funroll-loops -fno-asynchronous-unwind-tables -DNDEBUG"

echo "Building V1.4 (Production-Ready Sequential)"
echo "Git commit: $GIT_HASH"
echo "Flags: $FLAGS"

# Build without embedding compile flags in macro (avoid escaping issues)
g++ $FLAGS \
    -DGIT_HASH=\"$GIT_HASH\" \
    V1.4.cpp -o V1.4

if [ $? -eq 0 ]; then
    echo "✅ Build complete: V1.4"
    echo ""
    echo "Features:"
    echo "  - Overflow logging (128-bit exceeded)"
    echo "  - Fuse-hit logging (max_steps exceeded)"  
    echo "  - Reproducible metadata (JSON)"
    echo ""
    echo "Run with: ./V1.4 0 1000000 --tag <your_tag>"
else
    echo "❌ Build failed"
    exit 1
fi
