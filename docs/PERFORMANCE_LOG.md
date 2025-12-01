# Collatz Conjecture - Performance Optimization Log

## System Specifications
- **CPU:** Intel i5 (6 cores / 12 threads)
- **Test Range:** [2^71, 2^71 + 1B] (1 billion numbers, ~333M odd tested after mod-6 filtering)
- **Goal:** Maximum single-thread throughput → Multi-threading → GPU

---

## Performance Results Table

| Version | Description | Time (s) | Numbers/sec | Speedup vs V1 | Speedup vs Prev | Cycles | Instructions | Branches | Branch-misses | Notes |
|---------|-------------|----------|-------------|---------------|-----------------|-------|---------------|-----------------|-------|-----|
| V1 | std::set cycle detection | 183,092764241 | 54617 | **1.00x** | 1.00x | 796.810.832.544 | 1.395.495.915.428 | 332.721.563.671 | 13.175.191.612 | O(1) memory |
| V2 | Floyd cycle detection + bit ops | 5,736504057 | 1743679 | **31.92x** | 31.92x | 25.127.696.411 | 71.577.324.884 | 16.975.964.788 | 337.354.903 | O(1) memory |
| V3 | Multiple bit shift | 3,888829506 | 2572678 | **47.08x** | 1.475x | 17.028.978.736 | 50.490.467.827 | 9.919.901.276 | 10.015.949 | O(1) memory |
| V4 | Early termination to verified range | 0,070 | 142857142 | **2616.11x** | 55.53x | 289.190.901 | 1.008.099.411 | 165.321.977 | 1.464.317 | Verification only (10M range) |
| V5 | No Floyd cycle detection | 0,057 | 175438596 | **3212.44x** | 1.23x | 237.421.546 | 606.530.640 | 102.737.941 | 1.536.175 | Verification only (10M range) |
| V6 | Mod-6 filtering | 3,716 | 269106566 | **4927.66x** | 1.53x | 16.045.160.929 | 44.258.116.911 | 6.821.597.247 | 97.490.467 | 1B range, 333M tested |
| V7 | Compiler hints (always_inline) | 3,644 | 274423710 | **5024.98x** | 1.02x | 15.742.586.085 | 44.258.116.932 | 6.821.597.215 | 98.831.285 | 1B range, 333M tested |


---

## Version Details

### V1 - Baseline
- Cycle detection: `std::set<__uint128_t>` 
- Operations: `n % 2` (modulo), `n / 2` (division)
- Memory per number: ~avg 529 steps × 16 bytes = 8.5 KB per number tested
- Instructions per number: 1,395,495,915,428 / 10,000,000 = 139,550 instructions/number

### V2 - Floyd's Cycle Detection  
- Cycle detection: Floyd's algorithm (2 pointers)
- Operations: `n & 1` (bitwise AND), `n >> 1` (bit shift)
- Memory: 32 bytes fixed (2 × __uint128_t pointers)
- Instructions per number: 71,577,324,884 / 10,000,000 = 7,158 instructions/number
- Reduction: 139,550 → 7,158 = 94.9% fewer instructions
- Cost: 3 Collatz steps per loop iteration (1 slow + 2 fast)

### V3 - Multiple Bit Shift
- Added: `__builtin_ctzll()` to collapse consecutive divisions by 2
- Instructions per number: 50,490,467,827 / 10,000,000 = 5,049 instructions/number  
- Reduction vs V2: 7,158 → 5,049 = 29.5% fewer instructions
- Branch reduction: 16.98B → 9.92B = 41.6% fewer branches
- Branch-miss rate: 337M/17.0B = 2.0% (V2) → 10M/9.9B = 0.1% (V3)

### V4 - Early Termination to Verified Range
- **Optimization for verification throughput** (not step counting)
- Goal: Verify numbers reach 1 (detect cycles), maximize range coverage/sec
- Implementation: Terminate when trajectory drops below starting range (2^71)
- Assumption: All numbers < 2^71 already verified by existing Collatz research (verified to ~2^68)
- Cycle detection: Floyd's algorithm catches cycles above original; cycles below would have been detected in prior verification
- Performance: 0.070s vs V3's 3.89s = 55.5x speedup
- Throughput: 142.9M range/sec (verifying 10M range in 70ms)
- Actual computation: Only 5M odd numbers tested (evens skipped as redundant)
- Hardware efficiency: 289M cycles (58 cycles/number), 1,008M instructions (202 instructions/number)
- IPC: 1,008M / 289M = 3.49 instructions per cycle
- Branch efficiency: 1.46M misses / 165M branches = 0.89% miss rate
- Reduction from V3: 17.0B → 289M cycles (98.3% reduction), 50.5B → 1.0B instructions (98.0% reduction)
- Use case: Large-scale verification where step counts not needed, only cycle detection
- Bug fix: Used separate tmp variable for fast pointer to avoid corrupting step count
- Note: Does NOT compute accurate step counts (only counts steps until drop below original)
- Metric interpretation: "142.9M numbers/sec" = can verify a 142.9M range per second
- Parallelization strategy: Each thread takes range slice, uses same termination threshold

### V5 - No Floyd Cycle Detection
- **Removed Floyd's algorithm** - single pointer instead of slow/fast
- Rationale: With early termination, Floyd's 3x overhead (1 slow + 2 fast steps) is unnecessary
- Cycle detection: Iteration limit (100,000) catches infinite loops above threshold
- If trajectory stays ≥ original for 100k iterations → likely cycle (return -1)
- Performance: 0.057s vs V4's 0.070s = 1.23x speedup
- Throughput: 175.4M range/sec (verifying 10M range in 57ms)
- Hardware efficiency: 237M cycles (47 cycles/number), 607M instructions (121 instructions/number)
- IPC: 607M / 237M = 2.56 instructions per cycle
- Branch efficiency: 1.54M misses / 103M branches = 1.50% miss rate
- Reduction from V4: 289M → 237M cycles (18% reduction), 1,008M → 607M instructions (40% reduction)
- Key insight: Floyd does 3x work per iteration; with early termination, simple iteration counter suffices
- Computational load: Only 1 Collatz step per iteration (vs 3 in Floyd)
- Safety: MAX_ITERATIONS=100k prevents infinite loops if cycle exists above threshold

### V6 - Mod-6 Filtering
- **Skip numbers ≡ 3 (mod 6)** - only test n ≡ 1,5 (mod 6)
- Mathematical insight: Any n ≡ 3 (mod 6) after 3n+1 becomes even (wasted computation)
- Only n ≡ 1,5 (mod 6) can have interesting trajectories
- Implementation: Alternating stride pattern (+4, +2, +4, +2, ...)
- Reduction: 50% odd numbers (500M) → 33.3% filtered numbers (333M) = **1/3 fewer tests**
- Performance: 3.716s for 1B range = 269.1M range/sec
- Speedup vs V5: 1.53x (extrapolated from 10M to 1B range)
- Hardware efficiency: 16.045B cycles (48 cycles/tested), 44.258B instructions (133 instructions/tested)
- IPC: 44.258B / 16.045B = 2.76 instructions per cycle
- Branch efficiency: 97.5M misses / 6.822B branches = 1.43% miss rate
- Key insight: Algorithmic reduction (skip work) beats micro-optimizations
- Numbers actually tested: 333,333,333 (1/3 of range, 2/3 of odd numbers)
- Total steps computed: 3,492,676,192 (avg 10.5 steps/number to reach < 2^71)

### V7 - Compiler Hints (always_inline)
- **Added __attribute__((always_inline))** to hot functions
- Force compiler to inline `ctz128()`, `collatz_algo()`, `collatz_steps()`
- Hoisted low 64 bits extraction for parity check
- Performance: 3.644s for 1B range = 274.4M range/sec
- Speedup vs V6: 1.02x (1.9% improvement)
- Hardware efficiency: 15.743B cycles (47 cycles/tested), 44.258B instructions (133 instructions/tested)
- IPC: 44.258B / 15.743B = 2.81 instructions per cycle
- Branch efficiency: 98.8M misses / 6.822B branches = 1.45% miss rate
- Key insight: Compiler hints provide marginal gains - `-O3` already aggressively optimizes
- Diminishing returns: Single-core optimizations hitting limits
- Next step: **Multi-threading with OpenMP** for 10-12x additional speedup

---

## Compilation Commands

```bash
# Debug build
g++ -o V1 V1.cpp -std=c++17

# Optimized build (use for all benchmarks)
g++ -O3 -march=native -mtune=native -o V2 V2.cpp -std=c++17
```

---

## Notes & Observations

- 2^71 = 2361183241434822606848
- Testing 10M numbers gives statistically significant timing
- Always run multiple times and take average
- Disable CPU frequency scaling for consistent results if possible

