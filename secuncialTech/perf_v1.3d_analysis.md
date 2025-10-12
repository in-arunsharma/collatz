# V1.3d Performance Analysis

## Benchmark Results (3 runs each)

### V1.3c (Baseline - Precomputed Read-Only Memo)
- Run 1: 2,120,704 nums/sec
- Run 2: 2,084,623 nums/sec  
- Run 3: 2,057,978 nums/sec
- **Average: 2,087,768 nums/sec**

### V1.3d (Final Sequential Optimizations)
- Run 1: 2,253,914 nums/sec
- Run 2: 2,176,402 nums/sec
- Run 3: 2,134,476 nums/sec
- **Average: 2,188,264 nums/sec**

## Performance Gain

**Absolute improvement:** +100,496 nums/sec  
**Relative improvement:** +4.8% faster than V1.3c

## Micro-Optimizations Applied

1. **`__attribute__((always_inline))`** on `compute_collatz_readonly()`
   - Forces inlining of hot kernel function
   - Enables better loop-invariant code motion (LICM)
   - Reduces function call overhead

2. **`__restrict` pointer aliasing hints**
   - `const uint32_t* __restrict memo_ptr` parameter
   - Tells compiler memo array doesn't alias with other pointers
   - Enables more aggressive loop optimizations

3. **Hoisted low-limb extraction**
   ```cpp
   uint64_t lo = (uint64_t)current;
   if ((lo & 1ULL) == 0) {  // Even check on low bits only
   ```
   - Reduces 128-bit operations in common case
   - Faster parity checking

4. **Enhanced compile flags**
   - `-flto`: Link-time optimization (whole-program analysis)
   - `-fno-asynchronous-unwind-tables`: Smaller code, better cache
   - `-DNDEBUG`: Disable debug overhead

5. **Optional progress reporting**
   - Progress printing off by default (use `--progress`)
   - Eliminates fprintf() syscall overhead during benchmarking
   - Cleaner, more accurate timing

## Analysis

The 4.8% gain is respectable for micro-optimizations. Key factors:

1. **Always-inline:** Likely the biggest contributor - enables cross-function optimization
2. **LTO (-flto):** Whole-program analysis caught additional opportunities
3. **No unwind tables:** Reduced code size improves I-cache utilization
4. **Restrict hints:** Helps compiler generate better vectorized/pipelined code

## Comparison to V1.0 Baseline

- V1.0: 1,220,000 nums/sec (baseline)
- V1.3d: 2,188,264 nums/sec
- **Total speedup: 1.79× (79% faster)**

## Next Steps

**V1.3d is READY for parallelization!**

Sequential optimization journey complete:
- ✅ V1.1: CTZ bundling (+28%)
- ✅ V1.3: Early-exit memo (+45%)
- ✅ V1.3a: Path compression (+58%)
- ✅ V1.3c: Precomputed table (+63%)
- ✅ V1.3d: Micro-optimizations (+79%)

**Time to scale horizontally:**
- V1.4: OpenMP (80 cores → 175M+ nums/sec target)
- V1.5: CUDA (4 GPUs → 1B+ nums/sec target)
- V1.6: MPI (1,120 nodes → 2T+ nums/sec target)

## Build Command

```bash
g++ -O3 -march=native -std=c++17 -Wall -Wextra \
    -flto -fno-exceptions -fno-rtti -funroll-loops \
    -fno-asynchronous-unwind-tables -DNDEBUG \
    V1.3d.cpp -o V1.3d
```

## Usage

```bash
# Basic benchmark (no progress spam)
./V1.3d 0 1000000

# With progress reporting
./V1.3d 0 1000000 --progress

# Custom progress interval
./V1.3d 0 1000000 --progress 100000

# Save precomputed table
./V1.3d 0 1000000 --save steps_2p20.bin

# Load precomputed table
./V1.3d 0 1000000 --load steps_2p20.bin
```

---
**Date:** Saturday, October 12, 2025  
**Status:** PRODUCTION-READY - Final sequential version  
**Next:** Begin OpenMP parallelization (V1.4)
