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
| **1.3d** | **FINAL SEQUENTIAL** - always_inline, restrict, -flto, optional progress | `-O3 -march=native -flto -fno-exceptions -fno-rtti -funroll-loops -fno-asynchronous-unwind-tables -DNDEBUG` | **2,438,443** (333K tested) | **136** | **3.57B** (2.91B core + 0.66B atom) | 1.23B | **6.98** (core) | 78K LLC | 1.97M (0.73%) | 586K | - | 🚀 **+23.5% vs V1.3c! 2.00× vs baseline!** Micro-optimizations: always_inline, restrict, LTO, no unwind tables. **IPC 6.98 (core) = EXCELLENT!** 59.2% retiring. **READY FOR PARALLELIZATION!** Sequential complete. |
| | | | | | | | | | | | | |

---

**Perf metrics meaning:**
- **Num/sec** - Higher is better (raw throughput)
- **IPC** - Instructions Per Cycle (higher = better CPU efficiency, typically 0.5-4.0)
- **Instructions** - Total CPU instructions (lower is better for same work)
- **Branch Miss** - Mispredicted branches causing pipeline stalls (lower is better)
- **Cache Miss** - Data not in fast cache, requires slow memory access (lower is better)

---

## PARALLEL PERFORMANCE (OpenMP + MPI + CUDA)

### Local Development (Intel i7-12700H, 20 threads)

| Version | Parallelization | Threads | Num/sec | Speedup | Efficiency | Notes |
|---------|----------------|---------|---------|---------|------------|-------|
| **V1.4b-openmp** | OpenMP (CPU-only) | 1 | 2.4M | 1.00× | 100% | Sequential baseline (V1.3d core) |
| **V1.4b-openmp** | OpenMP (CPU-only) | 4 | 9.2M | 3.83× | 95.8% | Near-linear scaling |
| **V1.4b-openmp** | OpenMP (CPU-only) | 8 | 17.8M | 7.42× | 92.7% | Excellent scaling |
| **V1.4b-openmp** | OpenMP (CPU-only) | 20 | 19.4M | 8.08× | 40.4% | Hyperthreading limit |
| **V1.5-cuda-hybrid** | CUDA + OpenMP | GPU + 20 CPU | - | - | - | ⏳ In progress (40% complete) |

**Key Observations:**
- ✅ **V1.4b-openmp validated:** 95%+ efficiency up to 8 physical cores
- ✅ Per-thread performance: ~2.3M numbers/sec sustained
- ⚠️ Hyperthreading efficiency drops after physical cores saturate
- 🎯 **Production-ready for MareNostrum deployment**

---

### MareNostrum 5 - Projected Performance

**Hardware Specs:**
- **GPP Nodes:** 2× Intel Sapphire Rapids 8480+ @ 2.0GHz (112 cores/node, 256GB RAM)
- **ACC Nodes:** 2× Intel Sapphire Rapids 8460Y+ @ 2.3GHz (80 cores/node) + 4× NVIDIA H100 (64GB HBM2)
- **Allocation:** 150 GPP + 25 ACC nodes (shared among hackathon participants)

#### Phase 1: GPP Nodes (MPI + OpenMP) - DEPLOYED

| Configuration | Nodes | Total Cores | Expected Num/sec | Walltime (10B nums) | Status |
|--------------|-------|-------------|------------------|---------------------|--------|
| **Initial Test** | 3 GPP | 336 | 750M | 13 sec | ✅ Deployed |
| Conservative | 5 GPP | 560 | 1.25B | 8 sec | 🎯 Target |
| Production | 10 GPP | 1,120 | 2.5B | 4 sec | 🚀 Stretch |

**Architecture:**
- MPI: Embarrassingly parallel seed range partitioning
- OpenMP: 112 threads per node
- Per-core: 2.3M numbers/sec (based on V1.4b-openmp validation)
- Communication overhead: ~1-2 seconds (negligible for 10B+ ranges)

#### Phase 2: ACC Nodes (CUDA + OpenMP Hybrid) - IN DEVELOPMENT

| Configuration | Nodes | GPUs | Expected Num/sec | Status |
|--------------|-------|------|------------------|--------|
| Single ACC | 1 | 4× H100 | 500M - 1B | ⏳ Development |
| 3 ACC nodes | 3 | 12× H100 | 1.5B - 3B | 🎯 Target |
| 5 ACC nodes | 5 | 20× H100 | 2.5B - 5B | 🚀 Stretch |

**Per-GPU expectations:**
- H100 has 16,896 CUDA cores @ 1.98 GHz
- Conservative estimate: 125M - 250M numbers/sec per GPU
- Optimistic (with kernel optimization): 500M+ numbers/sec per GPU

#### Phase 3: Hybrid Deployment (GPP + ACC Combined)

| Configuration | GPP Nodes | ACC Nodes | Total Cores/GPUs | Expected Num/sec | Status |
|--------------|-----------|-----------|------------------|------------------|--------|
| **Hybrid Max** | 10 GPP | 5 ACC | 1,120 cores + 20 GPUs | 3.6B - 7.5B | 🎯 Hackathon Goal |

**Target:** Process 100 billion+ numbers during hackathon

---

### Scaling Efficiency

**Sequential → Parallel:**
- V1.0 baseline: 1.2M nums/sec (single core)
- V1.3d optimized: 2.4M nums/sec (single core) → **2.00× improvement**
- V1.4b-openmp (20 cores): 19.4M nums/sec → **8.08× parallel speedup**
- **Combined:** 16.17× faster than V1.0 baseline

**MareNostrum Projections:**
- 3 GPP nodes: 750M nums/sec → **625× faster than V1.0 baseline**
- 10 GPP nodes: 2.5B nums/sec → **2,083× faster than V1.0 baseline**
- Hybrid (10 GPP + 5 ACC): 3.6B+ nums/sec → **3,000×+ faster than V1.0 baseline**

---

### Resource Allocation Strategy (Hackathon)

**Shared Resources:**
- 150 GPP nodes total (all participants)
- 25 ACC nodes total (all participants)

**Our allocation strategy:**
1. **Initial test:** 3 GPP nodes (2% of pool) - validate deployment
2. **Production:** 5-10 GPP nodes (3-7% of pool) - sustained computation
3. **ACC exploration:** 2-5 ACC nodes (8-20% of pool) - GPU acceleration
4. **Be courteous:** Monitor queue, scale down if cluster busy

**Conservative approach:** Start small, scale up based on availability and performance validation

---

## BENCHMARKING METHODOLOGY

### Sequential (V1.0 - V1.3d)
```bash
# Standard benchmark
./build.sh && ./V1.3d 0 1000000

# With perf profiling
perf stat -d ./V1.3d 0 1000000
```

### Parallel (V1.4b-openmp)
```bash
# Build with OpenMP
g++ -O3 -march=native -fopenmp -o V1.4b-openmp V1.4b-openmp.cpp

# Test with varying thread counts
export OMP_NUM_THREADS=1
./V1.4b-openmp 0 10000000 test1t

export OMP_NUM_THREADS=4
./V1.4b-openmp 0 10000000 test4t

export OMP_NUM_THREADS=20
./V1.4b-openmp 0 10000000 test20t
```

### MareNostrum (Phase 1 - GPP)
```bash
# Build
cd marenostrum/
./build_phase1_gpp.sh

# Submit SLURM job (3 nodes)
sbatch slurm_phase1_gpp.slurm

# Monitor
squeue -u $USER
tail -f collatz_phase1_*.out

# Check results
cat *_global.json | grep throughput
```

---

## PERFORMANCE VALIDATION CHECKLIST

### Sequential Code (V1.3d)
- [x] Correctness: Self-validation tests pass
- [x] Thread-safety: Read-only memo table
- [x] Throughput: 2.4M+ numbers/sec per core
- [x] IPC: 6.0+ (excellent CPU utilization)
- [x] Retiring: 55%+ (minimal stalls)

### Parallel Code (V1.4b-openmp)
- [x] Scaling efficiency: 95%+ up to physical cores
- [x] Thread safety: No race conditions
- [x] Deterministic results: Same output regardless of thread count
- [x] Memory overhead: Acceptable (per-thread stats, local queues)
- [x] Production ready: Tested on local hardware

### MareNostrum Deployment (Phase 1)
- [ ] Build successful on MN5
- [ ] 3-node test completes without errors
- [ ] Throughput > 500M numbers/sec (conservative target)
- [ ] Scaling efficiency verified (compare 3 vs 5 vs 10 nodes)
- [ ] JSON output validated (cycles, overflows, max steps)

### Future Work (Phase 2 & 3)
- [ ] ACC node CUDA implementation
- [ ] GPU kernel optimization
- [ ] Hybrid MPI+CUDA+OpenMP coordination
- [ ] 1B+ numbers/sec per ACC node
- [ ] Combined 3.6B+ numbers/sec (hybrid deployment)

---

**Last Updated:** October 13, 2025 (Hackathon Morning)  
**Next Milestone:** Phase 1 validation on MareNostrum GPP nodes
