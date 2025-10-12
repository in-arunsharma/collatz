# Parallel Versions for MareNostrum 5

**Foundation:** V1.3c (validated, 2.07M nums/sec single-thread)

---

## Parallelization Roadmap

### V1.4 - OpenMP Multi-Threading (CPU)
**Target:** 80 cores per MareNostrum 5 node  
**Strategy:** Partition range across threads, each with read-only memo access  
**Expected:** 80× speedup → 150M+ nums/sec per node  
**Status:** In progress

### V1.5 - CUDA GPU Parallelization
**Target:** 4× NVIDIA Hopper H100 GPUs per node  
**Strategy:** Copy memo table to GPU constant memory, massive parallelism  
**Expected:** 300-500× speedup per GPU → 500M-1B nums/sec per GPU  
**Status:** Planned

### V1.6 - MPI Multi-Node Distribution
**Target:** 1,120 nodes × 4 GPUs = 4,480 total GPUs  
**Strategy:** Range partitioning across nodes, gather results  
**Expected:** 2+ trillion nums/sec on full MareNostrum 5  
**Status:** Planned

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
