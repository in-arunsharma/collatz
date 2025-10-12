# COLLATZ CONJECTURE - PERFORMANCE TRACKING

**Standard Benchmark:** 1,000,000 numbers starting from 2^71  
**Why 1M?** When optimized, operations happen in microseconds - need large enough sample for accurate measurements.

**Benchmark command:**
```bash
./V[X].0 0 1000000
```

**With perf (when enabled):**
```bash
sudo sysctl -w kernel.perf_event_paranoid=-1
perf stat -d ./V[X].0 0 1000000 2>&1 | tee perf_v[X].0.txt
```

---

## PERFORMANCE TABLE

| Ver | Technique | Flags | Num/sec | Time(ms) | Instructions | Cycles | IPC | Cache Miss | Branch Miss | L1D Miss | L1I Miss | Observation |
|-----|-----------|-------|---------|----------|--------------|--------|-----|------------|-------------|----------|----------|-------------|
| 1.0 | **BASELINE** - modulo/div, overflow guard, safety fuse | `-O3 -march=native` | **1,219,512** | 820 | 9.87B | 3.18B | 3.10 | 11.4K LLC | 26.5M (1.52%) | 113K | - | ✅ Good IPC (3.10), low branch miss (1.52%). Compiler did well. Bottleneck: modulo/div ops |
| 1.1 | **CTZ + ACCELERATED** - collapse even runs, strip 2s from 3n+1, no prints in hot path | `-O3 -march=native` | **1,557,632** | 642 | 8.44B | 2.71B | 3.11 | 11.4K LLC | 1.03M (0.065%) | 94K | - | ✅ 1.28x faster. 14% fewer instructions. 96% fewer branch misses. Clean hot path |
| 1.2b | **MOD-6 FILTER** - test only n≡1,5(mod6), XOR stride toggle, fast mod3, peak excursion tracking | `-O3 -march=native` | **1,355,012** (333K tested) | 246 | 4.35B | 1.64B | 2.65 | 9.4K LLC | 603K (0.09%) | 56K | - | ✅ Tests 1/3 of numbers (skips evens & mult-3). ~3× faster for **range coverage**. Per-number cost slightly higher due to stride overhead. Branch misses down to 0.09% |
| 1.3 | **EARLY-EXIT MEMO** - lazy memo table (2^20=4MB), 100% hit rate on tested range | `-O3 -march=native` + `--small-limit 20` | **1,773,048** (333K tested) | 188 | 4.28B | 1.22B | 3.51 | 11.2K LLC | 621K (0.11%) | 57K | - | ✅ **1.31× faster than V1.2b!** Memo table sweet spot: 2^20 (fits L2/L3). Tried 2^24 (64MB) = slower (cache pressure). IPC up to 3.51. 73% retiring (excellent). 100% memo hits |
| **1.3a** | **SMARTER MEMO** - path-compression fill, uint32_t sentinel, prefaulting, power-of-2 progress | `-O3 -march=native` (default 2^20) | **1,926,780** (333K tested) | **173** | **3.38B** | 1.24B | **4.99** | 11.8K LLC | 627K (0.11%) | 54K | - | ✅ **1.09× faster than V1.3!** Path compression = 21% fewer instructions. Prefaulting eliminates page-fault timing noise. IPC 4.99 (core), 59% retiring. **Best lazy-fill!** |
| 1.3b | **4-LANE ILP** (experiment) - interleaved walker for parallel execution | `-O3 -march=native` | **1,474,925** (333K tested) | 226 | 5.35B | 1.48B | 3.62 | 17K LLC | 1.91M (0.19%) | 82K | - | ❌ **23% SLOWER than V1.3a.** Lane overhead (58% more instructions) exceeded ILP benefit. Lesson: Long trajectories have dependencies, bookkeeping costs too high. **Parked.** |
| 1.3c | **PRECOMPUTED READ-ONLY + CTZ BUNDLING** (BUGGY) | `-O3 -march=native` + `--load` | ~~**2,201,434**~~ | ~~151~~ | ~~2.31B~~ | ~~640M~~ | ~~3.61~~ | ~~68K LLC~~ | ~~356K (0.11%)~~ | ~~85K~~ | - | ⚠️ **INCORRECT RESULTS!** 128-bit truncation + wrong backfill → corrupted memo table. **DO NOT USE.** |
| **1.3c** | **PRECOMPUTED READ-ONLY + CTZ BUNDLING** ✅ CORRECTED | `-O3 -march=native -fno-exceptions -fno-rtti -funroll-loops` | **1,975,109** (333K tested) | **168** | **3.07B** | 862M | **4.52** | 76K LLC | 987K (0.82%) | 591K | - | 🏆 **VALIDATED CORRECT! 1.62× vs baseline!** Fixed 128-bit truncation + backfill bugs. Validation self-test passes. 7.7% slower than buggy version but 100% correct. **Thread-safe for OpenMP/CUDA!** |
| **1.3d** | **FINAL SEQUENTIAL** - always_inline, restrict, -flto, optional progress | `-O3 -march=native -flto -fno-exceptions -fno-rtti -funroll-loops -fno-asynchronous-unwind-tables -DNDEBUG` | **2,188,264** (333K tested) | **152** | - | - | - | - | - | - | - | 🎯 **+4.8% vs V1.3c! 1.79× vs baseline!** Micro-optimizations: always_inline hot kernel, restrict aliasing hints, LTO, cleaner code size. **READY FOR PARALLELIZATION!** Sequential journey complete. |
| | | | | | | | | | | | | |

---

**Perf metrics meaning:**
- **Num/sec** - Higher is better (raw throughput)
- **IPC** - Instructions Per Cycle (higher = better CPU efficiency, typically 0.5-4.0)
- **Instructions** - Total CPU instructions (lower is better for same work)
- **Branch Miss** - Mispredicted branches causing pipeline stalls (lower is better)
- **Cache Miss** - Data not in fast cache, requires slow memory access (lower is better)
