#!/bin/bash
# Build script for V1.3d with enhanced optimization flags

g++ -O3 -march=native -std=c++17 \
    -Wall -Wextra \
    -flto \
    -fno-exceptions \
    -fno-rtti \
    -funroll-loops \
    -fno-asynchronous-unwind-tables \
    -DNDEBUG \
    V1.3d.cpp -o V1.3d

echo "Build complete. Binary: V1.3d"
echo "Enhanced flags: -flto -fno-asynchronous-unwind-tables -DNDEBUG"
