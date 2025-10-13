# Phase 1 GPP Integration Guide

## Overview
Phase 1 wraps the proven V1.4b-openmp sequential algorithm with MPI for multi-node scaling on MareNostrum 5 GPP nodes.

## Current Status
✅ MPI infrastructure complete (coordinator, workers, main program)  
✅ Build scripts and SLURM job submission ready  
⏳ **Need to link V1.4b-openmp core computation functions**  

## Integration Steps

### 1. Extract Core Functions from V1.4b-openmp.cpp

The following functions need to be made available to `worker_gpp_node.cpp`:

```cpp
// From parallel/V1.4b-openmp.cpp - extract these into a shared header

// Core Collatz computation
inline uint32_t collatz_steps(uint64_t n);

// Memoization structures
struct __attribute__((aligned(64))) HotValue {
    uint32_t steps;
    uint32_t hit_count;
};

extern HotValue* hot_memo;
extern std::atomic<uint32_t> hot_memo_size;

// Cycle detection
extern uint64_t* cycle_seeds;
extern std::atomic<uint32_t> cycle_count;

// Overflow tracking
extern uint64_t* overflow_seeds;
extern std::atomic<uint32_t> overflow_count;

// Initialize memoization tables
void initialize_memoization();

// Process a single number (returns steps, updates memos)
uint32_t process_number(uint64_t n);
```

### 2. Option A: Create Shared Header (Recommended)

Create `marenostrum/collatz_core.hpp`:

```cpp
#ifndef COLLATZ_CORE_HPP
#define COLLATZ_CORE_HPP

#include <cstdint>
#include <atomic>

// Include all core computation functions from V1.4b-openmp
// This keeps the proven sequential logic intact

// ... copy relevant structs and functions here ...

#endif
```

Then in `worker_gpp_node.cpp`:
```cpp
#include "collatz_core.hpp"  // Add this include
```

### 3. Option B: Direct Include (Quick & Dirty)

In `worker_gpp_node.cpp`, simply include the V1.4b-openmp implementation:

```cpp
// At top of worker_gpp_node.cpp
#define V14B_IMPLEMENTATION
#include "../parallel/V1.4b-openmp.cpp"  // Direct include
```

**Pros**: Fastest integration, no code duplication  
**Cons**: Less clean separation

### 4. Update worker_gpp_node.cpp Process Loop

Replace the placeholder in `GPPNodeWorker::process_work()`:

**Current placeholder:**
```cpp
// TODO: Replace with actual V1.4b-openmp core computation
uint32_t steps = (n % 100);  // Dummy computation
```

**Replace with:**
```cpp
// Call the proven V1.4b-openmp processing function
uint32_t steps = process_number(n);

// Update thread-local stats
if (steps > stats.max_steps) {
    stats.max_steps = steps;
    stats.max_steps_seed = n;
}
```

### 5. Initialize Memoization in Constructor

In `GPPNodeWorker::GPPNodeWorker()` constructor:

**Add after config initialization:**
```cpp
// Initialize memoization tables (from V1.4b-openmp)
initialize_memoization();
```

### 6. Quick Integration Script

Create `marenostrum/integrate_v14b.sh`:

```bash
#!/bin/bash
# Quick integration: symlink V1.4b core into marenostrum/

cd "$(dirname "$0")"

# Option 1: Symlink the proven code
ln -sf ../parallel/V1.4b-openmp.cpp collatz_core_impl.cpp

# Option 2: Or copy it if you want independence
# cp ../parallel/V1.4b-openmp.cpp collatz_core_impl.cpp

echo "✓ V1.4b-openmp core integrated"
echo "Now rebuild: ./build_phase1_gpp.sh"
```

## Testing Integration

After integration:

```bash
cd marenostrum/

# 1. Make scripts executable
chmod +x *.sh

# 2. Integrate V1.4b core
./integrate_v14b.sh  # (if you create it)

# 3. Build
./build_phase1_gpp.sh

# 4. Test locally with small MPI job
./test_local_mpi.sh 2  # 2 MPI ranks

# 5. If test passes, deploy to MareNostrum
```

## Expected Integration Time

- **Option A (Shared header)**: ~15 minutes (cleaner)
- **Option B (Direct include)**: ~5 minutes (faster)

## Verification Checklist

After integration, verify:

- [ ] Builds without errors
- [ ] Local MPI test runs (./test_local_mpi.sh)
- [ ] Output matches V1.4b-openmp behavior (cycle detection, overflows)
- [ ] Performance is similar to V1.4b per-core (19M nums/sec ÷ threads)
- [ ] JSON output files generated correctly

## Performance Expectations

**Local test (your laptop):**
- 4 MPI ranks × 4 threads = 16 parallel threads
- ~10-15M numbers/sec (depending on hardware)

**MareNostrum 5 GPP (5 nodes):**
- 5 MPI ranks × 112 threads = 560 parallel threads  
- Expected: **1.3B numbers/sec** (560 × 2.3M/sec)
- Walltime for 10B numbers: ~8 seconds

## Next Phase

Once Phase 1 is validated:
- **Phase 2**: Add ACC node support (CUDA + OpenMP hybrid)
- **Phase 3**: Combined deployment (5 ACC + 10 GPP = 3.6B nums/sec)

---

*This integration preserves the proven V1.4b-openmp sequential algorithm while scaling to multi-node MPI.*
