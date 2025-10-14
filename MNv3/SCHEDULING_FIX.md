# CRITICAL BUG FOUND: Static Scheduling Bottleneck

## The Problem

V1.5-openmp achieves **137M nums/sec** with:
```cpp
#pragma omp for schedule(dynamic, 512)
```

V1.5-mpi-lean gets only **2.48M nums/sec** with:
```cpp
#pragma omp for schedule(static, 8192)
```

## Root Cause

Collatz sequences have **highly variable computation time**:
- Easy numbers: ~50 steps (fast)
- Hard numbers: ~1400 steps (28× slower)

With **static scheduling**:
- Thread 0 gets indices [0, 8192) 
- Thread 1 gets indices [8192, 16384)
- etc.

If one chunk has many hard numbers, that thread becomes a **bottleneck**. Other threads finish and sit idle.

With **dynamic scheduling**:
- Threads grab work dynamically
- When a thread finishes its chunk, it immediately gets more work
- **Perfect load balancing** even with uneven work distribution

## The Fix

Change line 769 in V1.5-mpi-lean.cpp from:
```cpp
#pragma omp for schedule(static, 8192)
```

To:
```cpp
#pragma omp for schedule(dynamic, 512)
```

## Expected Result

This single line change should restore **137M nums/sec** performance!

## Why Your Friend's Code Had This Issue

Your friend optimized for:
1. Eliminating vector allocation (✅ good)
2. On-the-fly seed generation (✅ good)
3. Static scheduling for "less overhead" (❌ WRONG for uneven workload)

Static scheduling works great for **uniform** workloads (e.g., matrix operations).
For **highly variable** workloads like Collatz, dynamic scheduling is ESSENTIAL.

## Action Items

1. Change scheduling to `dynamic,512`
2. Rebuild
3. Retest - should see ~137M nums/sec
4. Then scale to 2+ nodes
