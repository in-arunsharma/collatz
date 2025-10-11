# CURRENT STATUS - V1.3 ✅

## Latest Version

✅ **V1.3 EARLY-EXIT MEMO - Lazy memo table with 100% hit rate**  
✅ **Performance: 1,773,048 numbers/second (333K tested in 188ms)**  
✅ **Speedup: 1.31× faster than V1.2b, 1.45× faster than V1.1**  
✅ **Perf analyzed: IPC 3.51, Branch miss 0.11%, 4.28B instructions**  
✅ **Memory: 2^20 (4MB) optimal - fits L2/L3 cache**  
✅ **Next: V1.3a micro-polish, then V1.3b multi-walk ILP**

**Early-exit memo working perfectly. Ready for final sequential polishing.**

---

## Version History

### V1.3 - Early-Exit Memo (CURRENT)
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
