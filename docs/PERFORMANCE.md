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
| | | | | | | | | | | | | |

---

**Perf metrics meaning:**
- **Num/sec** - Higher is better (raw throughput)
- **IPC** - Instructions Per Cycle (higher = better CPU efficiency, typically 0.5-4.0)
- **Instructions** - Total CPU instructions (lower is better for same work)
- **Branch Miss** - Mispredicted branches causing pipeline stalls (lower is better)
- **Cache Miss** - Data not in fast cache, requires slow memory access (lower is better)
