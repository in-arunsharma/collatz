# Parallel Versions for MareNostrum 5

**Foundation:** V1.3c (validated, 2.07M nums/sec single-thread)

---

## Collatz Conjecture Parallelization

This directory contains parallelized versions of the Collatz conjecture verification code, optimized for MareNostrum 5's architecture.

## V1.4b-openmp: OpenMP Multi-Core Parallelization

**Status:** ✅ **PRODUCTION READY** (Tested locally: 1/2/4/14/20 threads)

### Performance

- **Sequential baseline (V1.4b):** 2.38M nums/sec
- **1 thread:** 2.32M nums/sec (parity with sequential)
- **4 threads:** 5.40M nums/sec (2.33× speedup)
- **14 threads (P-cores):** 14.1M nums/sec (6.1× speedup) - ⭐ **Recommended for hybrid CPU**
- **20 threads (P+E):** 19.4M nums/sec (8.3× speedup) - lower efficiency due to E-cores
- **Target (80 cores MareNostrum):** ~150-192M nums/sec

**Tip for hybrid CPUs (Intel 12th gen+):** Use physical P-core count (14) instead of thread count (20) for best efficiency.

### Latest Features (Oct 13, 2025)

✅ **Output directory support:** All logs/metadata now in `output/` folder (no mixing with code)
✅ **Per-thread logs:** `output/*_<tag>_t<tid>.txt`
✅ **Hybrid CPU optimization:** Better scheduling for P-core + E-core architectures

### Building

```bash
cd parallel/
./build_v1.4b-openmp.sh
```

### Local Testing

```bash
# Optimal for hybrid CPU (physical P-cores only)
OMP_NUM_THREADS=14 ./V1.4b-openmp 0 10000000 --tag test --output-dir output

# All threads (includes E-cores, lower efficiency but higher absolute throughput)
OMP_NUM_THREADS=20 ./V1.4b-openmp 0 10000000 --tag test --output-dir output

# Check output
ls -lh output/
cat output/run_metadata_test.json
```

### MareNostrum 5 Deployment

#### Single-Node Job (80 cores)

```bash
# Submit single-node job
sbatch marenostrum_single_node.slurm

# Monitor
tail -f collatz_omp_<JOBID>.out

# Results will be in output/ directory on compute node
```

**Expected throughput:** 150-192M nums/sec

#### Multi-Node Array Job (100 nodes × 80 cores)

```bash
# Submit array job
sbatch marenostrum_array_job.slurm

# Check status
squeue -u $USER

# Aggregate results
cat output/run_metadata_mn5_job*_r*.json | jq -s 'map(.results.tested) | add'
```

**Expected throughput:** 15-19B nums/sec

---

## V1.5-cuda-hybrid: CUDA + OpenMP Hybrid ⚡ **NEW!**

**Status:** 🚧 **IN DEVELOPMENT** (Kernel compiles, full pipeline pending)

### Architecture

**GPU (CUDA):** Hot path only
- 128-bit arithmetic (manual uint2 implementation)
- Fuse=100K, no Brent, no logging
- Pure compute, minimal branching
- Output: `(steps, peak, status)` per seed

**CPU (OpenMP):** Cold queue processing
- Brent cycle detection
- 256-bit overflow handling
- Logging and metadata

**Pipeline:** Double-buffered batches
- GPU processes batch N
- CPU handles batch N-1 cold queues
- Overlapped compute for maximum throughput

### Local GPU Test (RTX 3060 Laptop)

```bash
./build_v1.5-cuda-hybrid.sh
./V1.5-cuda-hybrid

# Output:
# [CUDA] Device: NVIDIA GeForce RTX 3060 Laptop GPU
# [CUDA] Compute: 8.6
# [CUDA] Memory: 6.09 GB
# [CUDA] SMs: 30
```

**Expected throughput (when complete):**
- RTX 3060 Laptop: 100-500M nums/sec (10-50× CPU)
- H100 (MareNostrum 5): 5-50B nums/sec per GPU
- 4,480 GPUs: Trillions/sec (friend's recommendation)

### Implementation Checklist

- ✅ CUDA kernel skeleton (128-bit hot path)
- ✅ Build script (nvcc + OpenMP)
- ✅ GPU detection and info
- ⏳ Memo table upload to GPU
- ⏳ Batch pipeline with double buffering
- ⏳ CPU cold queue integration
- ⏳ Profiling (nvprof/nsight)

**Friend's Recommendation:**
> "Use CUDA for hot path on ACC nodes, OpenMP for CPU-only GPP nodes and all cold-queue work. Keep GPUs laser-focused on the cheap 128-bit hot kernel; let CPUs mop up rare/expensive paths."

---

## Development Timeline

- **V1.0-V1.3d:** Sequential optimizations (2.00× speedup)
- **V1.4:** Production logging
- **V1.4b:** Hot/cold queue architecture (2.38M nums/sec sequential)
- **V1.4b-openmp:** ✅ OpenMP parallelization (14.1M @ 14 cores, 19.4M @ 20 cores)
- **V1.5-cuda-hybrid:** 🚧 CUDA + OpenMP (in development, kernel ready)
- **V1.6-mpi:** Planned for post-hackathon multi-node

**Hackathon:** Monday, Oct 14, 2025 (TOMORROW!)

---

## Files

```
parallel/
├── V1.4b-openmp.cpp          # ✅ OpenMP (production ready)
├── build_v1.4b-openmp.sh     # ✅ Build script
├── V1.5-cuda-hybrid.cu       # 🚧 CUDA + OpenMP (in dev)
├── build_v1.5-cuda-hybrid.sh # 🚧 CUDA build
├── marenostrum_single_node.slurm  # SLURM: 80 cores
├── marenostrum_array_job.slurm    # SLURM: multi-node
├── output/                    # All logs & metadata
└── README.md                  # This file
```

---

*Last updated: Late Saturday night (Oct 13, 2025) - Final prep for Monday's hackathon!*


### Architecture

**Hot/Cold Queue Design:**
- **Hot path:** Pure 128-bit computation, 100K iteration fuse, no cycle detection
  * Per-thread statistics (ThreadStats struct)
  * Deterministic seed mapping (index → seed)
  * Static scheduling for reproducibility
  
- **Cold Queue 1 (fuse hits):** Extended 1M iteration fuse + Brent cycle detection
- **Cold Queue 2 (128-bit overflow):** 256-bit arithmetic + Brent cycle detection
- **Post-parallel processing:** Single-thread cold queue verification after hot scan

**Thread Safety:**
- Per-thread log files: `*_<tag>_t<tid>.txt`
- std::atomic counters for cold queue stats
- std::once_flag for initialization
- Read-only memo table (precomputed before parallel region)

### Building

```bash
cd parallel/
./build_v1.4b-openmp.sh
```

**Requirements:**
- g++ 7.5+ with OpenMP support
- Linux (uses `clock_gettime` for timing)

### Local Testing

```bash
# Single-threaded (verify matches V1.4b)
OMP_NUM_THREADS=1 ./V1.4b-openmp 0 1000000 --tag test1t

# Multi-threaded scaling test
OMP_NUM_THREADS=4 ./V1.4b-openmp 0 1000000 --tag test4t

# Determinism verification (compare max_steps across runs)
OMP_NUM_THREADS=1 ./V1.4b-openmp 0 1000000 --tag det1
OMP_NUM_THREADS=4 ./V1.4b-openmp 0 1000000 --tag det4
# Should report same max_steps value (e.g., 1242)
```

### MareNostrum 5 Deployment

#### Single-Node Job (80 cores)

```bash
# Submit single-node job
sbatch marenostrum_single_node.slurm

# Monitor progress
tail -f collatz_omp_<JOBID>.out

# Check results
ls -lh run_metadata_mn5_job*.json
ls -lh *_mn5_job*_t*.txt
```

**Expected throughput:** 150-192M nums/sec (2.4M × 80 cores × 0.8 efficiency)

#### Multi-Node Array Job (100 nodes × 80 cores)

```bash
# Submit 100-task array job (100B seeds total)
sbatch marenostrum_array_job.slurm

# Check array job status
squeue -u $USER

# Monitor specific task
tail -f collatz_omp_<JOBID>_<TASKID>.out

# Aggregate results
cat run_metadata_mn5_job*_r*.json | jq -s 'map(.results.tested) | add'
```

**Expected throughput:** 15-19B nums/sec (192M × 100 nodes)

#### Environment Variables

**Critical for NUMA performance:**
```bash
export OMP_NUM_THREADS=80         # Match --cpus-per-task
export OMP_PROC_BIND=close        # Bind threads to nearby cores
export OMP_PLACES=cores           # One thread per core (avoid HT)
export OMP_SCHEDULE=static        # Deterministic scheduling
```

**Optional tuning:**
```bash
export OMP_DYNAMIC=false          # Disable dynamic thread adjustment
export OMP_WAIT_POLICY=active     # Active spinning (low latency)
export OMP_MAX_ACTIVE_LEVELS=1    # No nested parallelism
```

### Verification

**Determinism:** All runs with same `[start, count)` should report:
- Same `max_steps` value
- Same total tested count
- Potentially different `hardest_n` (but with same max_steps)

**Correctness:**
- Compare cold queue triggers across thread counts
- Verify no cycles found (unless genuine!)
- Check 256-bit overflow counts (should be identical)

### Output Files

Per-thread logs (where `<tag>` = `--tag` parameter, `<tid>` = thread ID):
- `overflow_seeds_<tag>_t<tid>.txt` - 128-bit overflow cases
- `fuse_seeds_<tag>_t<tid>.txt` - Iteration fuse hits
- `cycle_seeds_<tag>_t<tid>.txt` - **GENUINE CYCLES** (if found!)
- `256bit_overflow_<tag>_t<tid>.txt` - 256-bit overflow cases

Global metadata:
- `run_metadata_<tag>.json` - Full reproducibility metadata

### Troubleshooting

**Low scaling efficiency (<50%):**
- Check NUMA binding: `numactl --hardware`
- Verify thread placement: `export OMP_DISPLAY_AFFINITY=true`
- Increase workload: Use `COUNT` > 10M seeds per thread

**Non-deterministic results:**
- Disable dynamic scheduling: `export OMP_SCHEDULE=static`
- Check for race conditions in logs (shouldn't happen with per-thread files)
- Verify same start/count across runs

**Memory issues:**
- Memo table: 4 MB (2^20 entries, shared read-only)
- Per-thread overhead: ~1 KB ThreadStats + cold queue vectors
- Total: <100 MB for 80 threads (negligible)

### Next Steps (CUDA GPU Parallelization - V1.5)

- **Target:** 4,480 NVIDIA Hopper GPUs on MareNostrum 5
- **Architecture:** GPU kernel for hot path (billions/sec per GPU)
- **Challenges:** 
  * GPU-compatible 256-bit arithmetic (already implemented!)
  * Host-side cold queue processing (CPU handles rare cases)
  * Multi-GPU coordination (NCCL for all-reduce)
- **Expected throughput:** Trillions of seeds/sec

**Current blocking issues:** None - V1.4b-openmp is production-ready for Monday's hackathon!

---

## Development Timeline

- **V1.0-V1.3d:** Sequential optimizations (2.00× speedup over naive)
- **V1.4:** Production logging infrastructure
- **V1.4b:** Hot/cold queue architecture (2.38M nums/sec)
- **V1.4b-openmp:** ✅ OpenMP parallelization (5.40M nums/sec @ 4 cores)
- **V1.5:** CUDA GPU parallelization (planned)
- **V1.6:** MPI multi-node (planned for post-hackathon)

**Hackathon deadline:** Monday, Oct 14, 2025 - MareNostrum 5 deployment

---

*Last updated: Saturday evening (pre-hackathon testing)*


---

## V1.3c Foundation - Production-Ready Sequential

**Performance:** 2,072,988 nums/sec (160ms for 333K numbers)

**Key Features:**
- ✅ **Thread-safe:** Read-only memo during compute (no locks needed)
- ✅ **Validated:** Self-test passes all known values
- ✅ **Overflow-safe:** MAX_SAFE guard prevents 128-bit overflow
- ✅ **CTZ bundling:** Bundles odd step + even collapse (1+k steps)
- ✅ **Mod-6 filtering:** Tests only n ≡ 1,5 (mod 6)
- ✅ **Path compression:** Backfills all visited small values
- ✅ **Portable:** Clean printf formats, bounds checking
- ✅ **Configurable:** --small-limit, --max-steps, --save/--load

**Critical Bug Fixes Applied:**
1. **128-bit truncation** - Fixed uint64_t → uint128_t in precompute
2. **Incorrect backfill** - Only backfill when trajectory completes
3. **Safety fuse on cache hits** - Honor max_steps even on memo hit

**Robustness:**
- Table size limits: 8 ≤ SMALL_LIMIT_BITS ≤ 28
- Validation guards against tiny tables
- Cumulative step tracking for exact backfill
- Read-only memo pointer hint for compiler optimization

---

## MareNostrum 5 Target Architecture

**ACC Partition Specs:**
- **Nodes:** 1,120
- **CPUs per node:** 2× Intel Sapphire Rapids 8460Y+ (40 cores each = 80 total)
- **GPUs per node:** 4× NVIDIA Hopper H100 (80GB each)
- **Total GPUs:** 4,480
- **Memory per node:** 512GB DDR5 + 256GB GPU (4×64GB)
- **Network:** 4× NDR200 per node (800 Gb/s total)
- **Peak Performance:** 260 PFlops

**Our Strategy:**
1. **V1.4 (OpenMP):** Baseline CPU parallelization (80 cores/node)
2. **V1.5 (CUDA):** GPU acceleration (4 GPUs/node, 4,480 total)
3. **V1.6 (MPI):** Multi-node scaling (1,120 nodes)

**Goal:** Push tested limits of Collatz conjecture from 2^71 using 4,480 GPUs

---

## Hackathon Timeline

**Date:** Monday, October 14, 2025  
**Time Remaining:** ~48 hours

**Saturday (Oct 12) - TODAY:**
- ✅ V1.3c validated and polished
- 🔄 V1.4 OpenMP implementation
- 🔄 V1.4 testing and benchmarking

**Sunday (Oct 13):**
- V1.5 CUDA kernel implementation
- V1.5 GPU benchmarking
- V1.6 MPI wrapper (basic)

**Monday (Oct 14) - DEPLOYMENT:**
- Deploy to MareNostrum 5
- Run large-scale tests
- Present results at hackathon

---

## Development Guidelines

**Correctness First:**
- Always validate memo table after precompute/load
- Test with known values before large runs
- Use MAX_SAFE overflow guard everywhere
- Honor safety fuse in all code paths

**Performance Second:**
- Profile before optimizing
- Measure parallel efficiency (speedup / cores)
- Watch for false sharing in OpenMP
- Minimize GPU memory transfers in CUDA

**Portability:**
- Use portable printf formats (%llu with casts)
- Test on different architectures
- Document architecture-specific optimizations
- Keep sequential V1.3c as reference

---

## Next Steps

1. **V1.4 OpenMP Implementation:**
   - Add `#pragma omp parallel for` to main loop
   - Thread-local accumulators for stats
   - Reduction for totals
   - Verify linear speedup (80 cores → 80× faster)

2. **Testing Strategy:**
   - Start with 2 threads, verify 2× speedup
   - Scale to 4, 8, 16, 40, 80 threads
   - Measure parallel efficiency at each level
   - Profile for bottlenecks (cache contention, false sharing)

3. **Benchmarking:**
   - Compare V1.4 vs V1.3c (single-thread)
   - Measure speedup curve (threads vs throughput)
   - Validate results match sequential version
   - Document performance characteristics

---

**Status:** V1.3c foundation ready. Beginning V1.4 OpenMP implementation now.
