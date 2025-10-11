# Collatz Conjecture - MareNostrum 5 Optimization Project

**Goal:** Push tested limits from 2^71 using MareNostrum 5's 4,480 NVIDIA Hopper GPUs

---

## Quick Start

```bash
cd secuncialTech
./build.sh all                    # Build all versions
sudo perf stat ./V1.0 0 1000000  # Benchmark with perf
```

---

## Optimization Strategy

### V1.0 - Baseline (DONE)
**Bottleneck:** Modulo (%) and division (/) operations are SLOW  
**Perf Data (-O0):** 136M branch misses (5.4%), 1.89 IPC, 533K nums/sec  
**Perf Data (-O3):** TBD (need perf access), 1.21M nums/sec (2.26x improvement!)  
**Performance:** 1.21M numbers/sec with `-O3 -march=native`

**Code:**
```cpp
if (n % 2 == 0)    // Modulo is expensive!
    n = n / 2;     // Division is expensive!
else
    n = 3 * n + 1;
```

**Observation:** Compiler optimizations alone gave 2.26x speedup. Algorithm optimizations should give much more.

---

### V2.0 - Bitwise Operations
**Fix:** Replace expensive ops with bit operations  
**Expected:** 2-3x speedup

**Changes:**
- `n % 2` → `n & 1` (AND operation, 1 cycle vs ~10-40 cycles)
- `n / 2` → `n >> 1` (bit shift, 1 cycle vs ~10-40 cycles)
- Keep `3*n+1` for now (multiplication still reasonably fast)

**Why it works:**  
- Modulo by power of 2 = check least significant bit
- Division by 2 = right shift by 1 position
- CPU can do bitwise ops in 1 cycle, division takes 10-40+ cycles

**Target:** Reduce branch misses, increase IPC

---

### V3.0 - Compiler Optimizations
**Fix:** Let compiler optimize our code  
**Expected:** 2-5x speedup on top of V2.0

**Flags:**
- `-O3`: Aggressive optimizations (loop unrolling, inlining, etc.)
- `-march=native`: Use all CPU features (AVX2, AVX512, etc.)
- `-funroll-loops`: Unroll loops manually

**Why it works:**  
- Compiler can see patterns we can't
- Auto-vectorization (SIMD)
- Better register allocation
- Instruction scheduling

**Target:** Max out single-threaded CPU performance

---

### V4.0 - Skip Even Chains
**Fix:** When even, divide by 2^k at once instead of k times  
**Expected:** 1.5-2x speedup

**Algorithm:**
```cpp
if (n is even) {
    k = count_trailing_zeros(n);  // How many times divisible by 2
    n >>= k;                       // Divide by 2^k in one go
    steps += k;
}
```

**Why it works:**  
- 16 → 8 → 4 → 2 → 1 becomes 16 >> 4 = 1 (one operation)
- Fewer loop iterations
- Better branch prediction

**Use:** `__builtin_ctzll()` or similar intrinsic

---

### V5.0 - Combined Odd Operation
**Fix:** Odd numbers always produce even numbers  
**Expected:** 1.3-1.5x speedup

**Algorithm:**
```cpp
if (n & 1) {
    n = (3*n + 1) >> 1;  // Combine 3n+1 and immediate /2
    steps += 2;           // Count both operations
}
```

**Why it works:**  
- 3*odd + 1 = even (always!)
- Skip one loop iteration per odd number
- One operation instead of two

---

### V6.0 - Path Compression/Memoization
**Fix:** Cache results for numbers we've seen  
**Expected:** 1.2-3x (depends on hit rate)

**Challenge at 2^71 scale:**  
- Can't cache everything (need petabytes)
- Cache only small numbers (< 2^32)
- Use hash table for recent values

**Tradeoff:** Memory vs compute time

---

### V8.0 - OpenMP Multi-threading
**Fix:** Use all CPU cores  
**Expected:** 8-16x (linear with core count)

**Code:**
```cpp
#pragma omp parallel for reduction(+:total_steps)
for (uint128_t n = start; n < end; n++) {
    // Each thread processes different numbers
}
```

**Why it works:**  
- Numbers are independent (embarrassingly parallel)
- No synchronization needed during computation
- MareNostrum has 80 cores per node!

**Compile:** `-fopenmp`

---

### V11.0 - CUDA GPU (Critical for Hackathon!)
**Fix:** Massive parallelization on NVIDIA Hopper GPUs  
**Expected:** 100-1000x per GPU

**Architecture:**
```
1 GPU = thousands of cores
Each thread = one number
4 GPUs/node × 1120 nodes = 4,480 GPUs
```

**Challenges:**
- uint128_t not native in CUDA (need custom type)
- Memory transfers (minimize!)
- Divergent branches (odd/even paths)

**Strategy:**
- Kernel processes batch of numbers
- Each thread independent
- Minimize memory access
- Use shared memory for stats

**This is the game changer for MareNostrum!**

---

### V12.0 - MPI Multi-node
**Fix:** Distribute across nodes  
**Expected:** Linear scaling

**Architecture:**
```
Node 0: tests 2^71 + 0 to 1B
Node 1: tests 2^71 + 1B to 2B
...
Node N: tests 2^71 + N*1B to (N+1)*1B
```

**MPI Operations:**
- `MPI_Init()` / `MPI_Finalize()`
- `MPI_Comm_rank()` - which node am I?
- `MPI_Comm_size()` - how many nodes total?
- `MPI_Reduce()` - aggregate results

---

### V13.0 - Hybrid CPU+GPU+MPI
**Ultimate version:** Combine everything

**Architecture:**
```
1,120 nodes
├─ Each node: 4x Hopper GPUs + 80 CPU cores
├─ GPUs: Main computation
├─ CPUs: Coordination, data prep
└─ MPI: Cross-node communication

Total: 4,480 GPUs working in parallel!
```

**Estimated Performance:**  
If 1 GPU = 1B numbers/sec  
Then 4,480 GPUs = 4.48 TRILLION numbers/sec

**Could test:** Entire range 2^71 to 2^72 in hours!

---

## Performance Tracking

See `secuncialTech/PERFORMANCE.md` for detailed metrics.

**Key Metrics:**
- **Throughput:** Numbers/second (higher = better)
- **IPC:** Instructions per cycle (shows CPU efficiency)
- **Branch Miss:** Mispredictions (causes stalls)
- **Cache Miss:** Memory access delays

---

## MareNostrum 5 Target

**Accelerated Partition:**
- 1,120 nodes
- 4× NVIDIA Hopper H100 GPUs per node
- 2× Intel Sapphire Rapids 8460Y+ (80 cores)
- 512GB DDR5 + 256GB GPU memory per node
- 4× NDR200 network (800Gb/s)

**Why Hopper is Perfect:**
- Massive parallelism (thousands of CUDA cores)
- FP64 performance (for large integer operations)
- Fast HBM2 memory
- PCIe 5.0 for fast CPU-GPU transfer

---

## File Structure

```
MN25/
├── README.md (this file)
├── MareNostrumINFO/
│   └── marenostrum5_specs.md
└── secuncialTech/
    ├── V1.0.cpp, V2.0.cpp, ...
    ├── PERFORMANCE.md (benchmark table)
    ├── build.sh (compilation script)
    └── perf_v*.txt (perf stat outputs)
```

---

## Development Approach

1. **Optimize sequentially first** (V1-V7)
   - Understand bottlenecks with `perf stat`
   - Each version fixes specific bottleneck
   - Measure everything!

2. **Then parallelize** (V8-V13)
   - OpenMP for local multi-core
   - CUDA for GPUs (the real power)
   - MPI for multi-node scaling

3. **Focus on performance, not polish**
   - Measure with perf
   - Identify bottleneck
   - Fix it
   - Repeat

---

## Key Commands

```bash
# Build
cd secuncialTech
g++ -O3 -march=native -std=c++17 -o V3.0 V3.0.cpp

# Benchmark with perf
sudo perf stat -e instructions,cycles,cache-misses,branch-misses ./V3.0 0 1000000

# Detailed profiling
sudo perf record ./V3.0 0 1000000
sudo perf report

# OpenMP version
g++ -O3 -march=native -fopenmp -o V8.0 V8.0.cpp
OMP_NUM_THREADS=16 ./V8.0 0 1000000

# CUDA version (on MareNostrum)
nvcc -O3 -arch=sm_90 -o V11.0 V11.0.cu
./V11.0 0 1000000000
```

---

## Bottleneck Analysis

**Current V1.0 Bottlenecks:**
1. ❌ **Branch misses: 136M (5.4%)** - CPU pipeline stalls
2. ❌ **Modulo operation** - 10-40 cycles per call
3. ❌ **Division operation** - 10-40 cycles per call
4. ❌ **No compiler optimization** - Missing auto-vectorization
5. ❌ **Single threaded** - Using 1 core out of many
6. ❌ **Sequential** - Not using GPU at all

**Fix order:**
- V2: Remove modulo/division → fix #2, #3
- V3: Compiler opts → fix #4
- V4-V5: Algorithm → reduce iterations
- V8: OpenMP → fix #5
- V11: CUDA → fix #6 (biggest impact!)

---

## Success Metrics for Hackathon

**Minimum:** Test range beyond any previous attempt from 2^71  
**Target:** Test 2^71 to 2^72 (double the range)  
**Stretch:** Test 2^71 to 2^73 or beyond

**With 4,480 GPUs at ~1B numbers/sec each:**  
Could process ~4.5 trillion numbers per second!

---

**Focus: Performance. Just performance.**
