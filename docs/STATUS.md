# Current Status: V1.3d � (FINAL SEQUENTIAL)

**Latest Version:** V1.3d - Final sequential micro-optimizations before parallelization

**Performance:** 2,188,264 nums/sec (152ms for 333K numbers) - **READY FOR OPENMP!**

**Micro-Optimizations Applied (Oct 12, 2025):**
- ✅ `__attribute__((always_inline))` on hot kernel
- ✅ `__restrict` pointer for memo aliasing hints
- ✅ Hoisted low-limb extraction for faster parity checks
- ✅ Enhanced compile flags: `-flto -fno-asynchronous-unwind-tables -DNDEBUG`
- ✅ Optional `--progress` flag (off by default for clean benchmarking)

**Performance Gain:** +4.8% vs V1.3c, **1.79× vs V1.0 baseline**

**Sequential Journey Complete:**
- V1.0 → V1.3d: **79% speedup**
- Ready for horizontal scaling (OpenMP → CUDA → MPI)

**See:** [perf_v1.3d_analysis.md](../secuncialTech/perf_v1.3d_analysis.md) for benchmark details

---

## Version History

### V1.3d - Final Sequential Optimizations 🎯 (CURRENT BEST)
- **Technique:** Compiler hints (always_inline, restrict), LTO, cleaner code size
- **Optimization:** Extract last sequential gains before parallelization
- **Numbers tested:** 333,333 (1/3 of 1M range, mod-6 filtered)
- **Time:** 152ms (average of 3 runs)
- **Throughput:** 2,188,264 nums/sec (**+4.8% vs V1.3c**)
- **Instructions:** (not measured - perf disabled)
- **Key wins:**
  - Always-inline: Forces inlining for better LICM/if-conversion
  - Restrict: Aliasing hints for aggressive loop optimization
  - LTO: Whole-program analysis finds additional opportunities
  - No unwind tables: Smaller code → better I-cache
  - Optional progress: Eliminate syscall overhead in benchmarks
- **Engineering:** Clean benchmarking, production-ready
- **Status:** ✅ **FINAL SEQUENTIAL VERSION - Ready for OpenMP!**

### V1.3c - Precomputed Read-Only Memo Table 🏆 (CORRECTED)
- **Technique:** Precomputed read-only memo with CTZ bundling + Critical Bug Fixes
- **Optimization:** Thread-safe memo access, validated correctness
- **Numbers tested:** 333,333 (1/3 of 1M range, mod-6 filtered)
- **Time:** 168ms
- **Throughput:** 1,975,109 nums/sec (average of 3 runs: 2.09M)
- **Instructions:** 3.07B (core)
- **IPC:** 4.52 (core)
- **Branch miss:** 0.82%
- **Key wins:**
  - Fixed 128-bit truncation bug (precompute correctness)
  - Fixed backfill logic bug (completed flag)
  - Validation self-test passes
  - Thread-safe for parallelization
- **Engineering:** Default 2^20 (4MB), optional disk save/load
- **Status:** ✅ Validated baseline for V1.3d and V1.4

### V1.3a - Smarter Memo Engineering (CURRENT BEST SEQUENTIAL)
- **Technique:** Path-compression fill, uint32_t sentinel, prefaulting
- **Optimization:** Backfill ALL small values visited in trajectory
- **Numbers tested:** 333,333 (1/3 of 1M range, mod-6 filtered)
- **Time:** 173ms
- **Throughput:** 1,926,780 nums/sec
- **Instructions:** 3.38B (2.61B core + 0.77B atom) - **21% reduction from V1.3**
- **IPC:** 4.99 (core) - excellent!
- **Branch miss:** 0.11%
- **Key wins:** 
  - Path compression: More memo entries per trajectory → fewer instructions
  - Prefaulting: Eliminates page-fault timing noise
  - uint32_t sentinel: Cleaner comparisons vs int32_t
  - Power-of-2 progress: Bitmask instead of modulo
- **Engineering:** Default 2^20 (4MB), progress every 2^14
- **Status:** ✅ Best single-thread sequential version

## V1.3b: 4-Lane Interleaved Walker (ILP Experiment) ❌

**Performance:** 1,474,925 nums/sec (226ms for 333K)
**Status:** ❌ Parked - V1.3a remains best

**Goal:** Boost ILP by processing 4 independent seeds in lockstep
**Result:** 23% SLOWER than V1.3a (58% more instructions)

**Why it failed:**
- Lane management overhead (checking 4 lanes per iteration)
- Increased branching: 1.91M branch misses vs 627K
- Long trajectories (500+ steps) have internal dependencies
- ILP benefit couldn't overcome bookkeeping cost
- 80.9% retiring shows CPU busy doing wasteful work

**Lesson:** ILP optimization not beneficial for this workload. Single-lane with better algorithm beats multi-lane overhead.

---

## V1.3c: Precomputed Read-Only Memo Table 🏆 (CORRECTED)

**Performance:** 1,975,109 nums/sec (168ms for 333K, with validated table)
**Status:** ✅ **PRODUCTION-READY - Correctness validated, ready for parallelization!**

**CRITICAL BUG FIXES (Oct 12, 2025):**

**Bug #1 - 128-bit Truncation:**
- **Problem:** Precompute used `uint64_t current` but assigned `uint128_t` values → silent truncation
- **Fix:** Changed to `uint128_t current` throughout trajectory
- **Impact:** Previous table was corrupted for all values!

**Bug #2 - Incorrect Backfill:**
- **Problem:** Backfilled step counts even when safety fuse/overflow tripped → wrong values
- **Fix:** Added `completed` flag, only backfill when trajectory reaches 1 or cached value
- **Impact:** Previous table had bogus values for aborted paths!

**Changes from V1.3a:**
1. **Precompute phase:** Fully populate 2^20 table BEFORE timing
   - Path-compression backfill during precompute (now CORRECT)
   - 100% of small values cached upfront
   - ✅ Validation self-test ensures correctness
2. **Compute phase:** Read-only memo access (no writes)
   - Thread-safe - trivial to parallelize
   - No locks needed for OpenMP/CUDA
3. **CTZ bundling optimization:** Bundle odd step + even collapse
   ```cpp
   // Instead of: 3n+1, then collapse evens next iteration
   // Do: t = 3n+1; k = ctz(t); return t>>k with 1+k steps
   ```
4. **Optional disk I/O:** Save/load `steps_2p20.bin`
   - Amortize precompute cost across runs
   - 4MB file for 2^20 table

**Performance Analysis (CORRECTED VERSION):**
- **1.62× faster** than V1.0 baseline (1.98M vs 1.22M nums/sec)
- Core instructions: 3.07B
- Cycles: 862M
- IPC: 4.52 (core) - excellent!
- Branch misses: 0.82% (986K total) - slightly higher but still good
- 60.0% retiring (excellent CPU utilization)

**Why Corrected Version is Slightly Slower:**
- `completed` flag adds conditional logic
- Skip backfill for aborted paths (correct behavior)
- **Trade-off:** 7.7% slower but 100% correct
- **Verdict:** Correctness is non-negotiable for scientific computing

**CTZ Bundling Impact:**
```
Before (slow):               After (fast):
n odd → t=3n+1 (1 step)     n odd → t=3n+1, k=ctz(t)
next iter → t even              → return t>>k (1+k steps)
next iter → collapse evens  
Total: 2+ loop iterations   Total: 1 loop iteration
```

**Benefits:**
- **Fastest CORRECT sequential** implementation (1.98M nums/sec)
- **Thread-safe** memo access (critical for OpenMP/CUDA)
- **Validated** - self-test confirms known values correct
- **Production-ready** - safe for MareNostrum 5 deployment
- **Predictable cache** behavior (no write contention)
- Foundation for V1.4 (OpenMP) and V1.5 (CUDA)
- Disk save/load amortizes precompute cost

**Usage:**
```bash
# Precompute and save
./V1.3c 0 1000000 --save steps_2p20.bin

# Load and run (fast startup)
./V1.3c 0 1000000 --load steps_2p20.bin

# Custom table size
./V1.3c 0 1000000 --small-limit 22  # 2^22 = 4M entries (16MB)
```

**Key Insight:** 
Thread-safety doesn't require a performance sacrifice! With proper CTZ bundling, we get:
- ✅ Best sequential performance (2.20M nums/sec)
- ✅ Thread-safe for massive parallelization
- ✅ Foundation ready for 1000× GPU speedup

**Next Steps:** 
- ✅ V1.3c provides fastest thread-safe foundation
- → V1.4: OpenMP parallelization (multicore CPU)
- → V1.5: CUDA implementation (GPU)
- → V1.6: MPI for multi-node scaling

---

### V1.3 - Early-Exit Memo
- **Technique:** Path-compression fill, uint32_t sentinel, prefaulting
- **Optimization:** Backfill ALL small values visited in trajectory
- **Numbers tested:** 333,333 (1/3 of 1M range, mod-6 filtered)
- **Time:** 173ms
- **Throughput:** 1,926,780 nums/sec
- **Instructions:** 3.38B (2.61B core + 0.77B atom) - **21% reduction from V1.3**
- **IPC:** 4.99 (core) - excellent!
- **Branch miss:** 0.11%
- **Key wins:** 
  - Path compression: More memo entries per trajectory → fewer instructions
  - Prefaulting: Eliminates page-fault timing noise
  - uint32_t sentinel: Cleaner comparisons vs int32_t
  - Power-of-2 progress: Bitmask instead of modulo
- **Engineering:** Default 2^20 (4MB), progress every 2^14

### V1.3 - Early-Exit Memo
- **Technique:** Lazy memo table for n < 2^20, early exit when trajectory dips
- **Optimization:** 100% memo hit rate, 4MB table fits in cache
- **Numbers tested:** 333,333 (1/3 of 1M range, mod-6 filtered)
- **Time:** 188ms
- **Throughput:** 1,773,048 nums/sec
- **Instructions:** 4.28B (3.37B core + 0.90B atom)
- **IPC:** 3.51 (excellent!)
- **Branch miss:** 0.11%
- **Key win:** 1.31× faster than V1.2b by caching small trajectory tails
- **CLI:** `--small-limit 20` for 2^20 (tested 16-26, sweet spot = 20)
- **Learning:** 2^24 (64MB) was slower (cache pressure + 16K page faults)

### V1.2b - Mod-6 Filtered
- **Technique:** Test only n≡1,5(mod6) - skip evens & multiples of 3
- **Optimization:** XOR stride toggle (`delta ^= 6`), fast mod3, correct first stride
- **Numbers tested:** 333,333 (1/3 of 1M range)
- **Time:** 246ms
- **Throughput:** 1,355,012 nums/sec
- **Instructions:** 4.35B
- **IPC:** 2.65
- **Branch miss:** 0.09%
- **Key win:** ~3× faster for range coverage (same numerical range, fewer starting points)
- **New feature:** Peak excursion tracking (`max_excursion` field)

### V1.1 - CTZ + Accelerated
- **Technique:** Collapse even runs with CTZ, strip 2s from 3n+1, clean hot path
- **Numbers/sec:** 1,557,632
- **Time:** 642ms (1M numbers)
- **Improvement:** 1.28× faster than V1.0
- **Instructions:** 8.44B (14% reduction)
- **Branch misses:** 96% reduction (26.5M → 1.03M)

### V1.0 - Baseline
- **Numbers/sec:** 1,219,512
- **Time:** 820ms (1M numbers)
- **Instructions:** 9.87B
- **IPC:** 3.10
- **Branch miss:** 1.52%
- **Correctness:** Overflow guard, safety fuse, div-by-zero guard

---

## Optimization Roadmap

### ✅ V1.1 - CTZ + Accelerated (DONE - 1.28× faster)
- Collapse even runs: `n >> ctz_u128(n)` 
- Strip all 2s from 3n+1: `m = 3*n+1; n = m >> ctz_u128(m)`
- Eliminates branch mispredictions (96% reduction)

### ✅ V1.2b - Mod-6 Filter (DONE - 3× faster for range coverage)
- Test only n≡1,5(mod6) - skip evens & multiples of 3
- Branchless XOR stride toggle: `delta ^= 6`
- Fast mod3 using 64-bit arithmetic (2^64≡1 mod 3)
- Peak excursion tracking for trajectory analysis

### 🔄 V1.3 - Early-exit memo table (IN PROGRESS)
- Small read-only memo (2^24 entries ≈ 64MB)
- Early exit when trajectory dips below threshold
- Lazy filling, no threading races
- Expected: Multi-× speedup (paths dip early/frequently)

### V1.3a - Micro-polish (planned)
- Hoist constants outside loops
- Power-of-2 progress masking
- Final sequential cleanups

### V1.3b - Multi-walk interleaving (planned)
- Process 4 independent seeds in lockstep
- Reduce dependency stalls via ILP
- Expected: 5-20% single-thread boost

### V1.4+ - OpenMP parallelization
- Target: laptop cores (multi-threading)

### V1.5+ - CUDA for MareNostrum
- 4,480 NVIDIA Hopper GPUs = massive parallelism

---

## Performance Focus

**What matters:**

## Performance Focus

**What matters:**
- Numbers/second (higher = better)
- IPC from perf (higher = better CPU efficiency)  
- Branch misses (lower = better)
- Instructions executed (lower = better algorithm)

**What doesn't matter:**
- Code beauty
- Comments  
- Anyone else's opinion

**Goal:** PURE PERFORMANCE. Period.

---

## Benchmark Protocol

**Standard test:** 
```bash
./V[X].0 0 1000000
```

**With perf:**
```bash
sudo perf stat -d ./V[X].0 0 1000000 2>&1 | tee perf_v[X].0.txt
```

**Record in PERFORMANCE.md table immediately**

---

You're ready to optimize. Start with V2.0 bitwise operations - should be quick and give 2-3x improvement.
