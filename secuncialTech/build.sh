#!/bin/bash
# Build script for Collatz Conjecture versions

# Colors for output
GREEN='\033[0;32m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "${BLUE}=== Collatz Conjecture Build Script ===${NC}"
echo ""

# Get version from argument or default to all
VERSION=${1:-"all"}

build_version() {
    local ver=$1
    local opts=$2
    local desc=$3
    
    if [ -f "V${ver}.cpp" ]; then
        echo -e "${BLUE}Building V${ver} ($desc)...${NC}"
        g++ $opts -std=c++17 -o "V${ver}" "V${ver}.cpp"
        if [ $? -eq 0 ]; then
            echo -e "${GREEN}✓ V${ver} built successfully${NC}"
        else
            echo "✗ V${ver} build failed"
            return 1
        fi
        echo ""
    fi
}

if [ "$VERSION" == "all" ]; then
    # Build all versions
    build_version "1.0" "-O0" "Baseline - No optimizations"
    build_version "2.0" "-O0" "Bitwise operations"
    build_version "3.0" "-O3 -march=native -funroll-loops" "Compiler optimizations"
    build_version "4.0" "-O3 -march=native" "Skip even chains"
    build_version "5.0" "-O3 -march=native" "Optimized 3n+1"
    build_version "6.0" "-O3 -march=native" "Path compression"
    build_version "8.0" "-O3 -march=native -fopenmp" "OpenMP parallel"
else
    # Build specific version
    case $VERSION in
        "1.0")
            build_version "1.0" "-O0" "Baseline - No optimizations"
            ;;
        "2.0")
            build_version "2.0" "-O0" "Bitwise operations"
            ;;
        "3.0")
            build_version "3.0" "-O3 -march=native -funroll-loops" "Compiler optimizations"
            ;;
        "4.0")
            build_version "4.0" "-O3 -march=native" "Skip even chains"
            ;;
        "5.0")
            build_version "5.0" "-O3 -march=native" "Optimized 3n+1"
            ;;
        "6.0")
            build_version "6.0" "-O3 -march=native" "Path compression"
            ;;
        "8.0")
            build_version "8.0" "-O3 -march=native -fopenmp" "OpenMP parallel"
            ;;
        *)
            echo "Unknown version: $VERSION"
            echo "Usage: ./build.sh [version|all]"
            echo "Example: ./build.sh 1.0"
            echo "         ./build.sh all"
            exit 1
            ;;
    esac
fi

echo -e "${GREEN}Build complete!${NC}"
echo ""
echo "Usage examples:"
echo "  ./V1.0 0 10000          # Test 10K numbers from 2^71"
echo "  ./V1.0 0 1000000        # Test 1M numbers from 2^71"
echo "  ./V1.0 1000000 1000000  # Test 1M numbers from 2^71 + 1M offset"
