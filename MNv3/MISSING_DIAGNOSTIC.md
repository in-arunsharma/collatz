## 🔍 Missing Diagnostic Output - Added More Debug

### What We Found

✅ OpenMP **IS** linked (`libgomp.so.1`)
✅ Environment shows `threads/rank=112`
❌ But diagnostic message `[RANK 0] OpenMP started with N threads` **IS MISSING**

This means either:
1. The parallel region never executes
2. The diagnostic print is being skipped
3. stderr is being lost somehow

### New Diagnostic Build

I've added MORE diagnostics:
- **BEFORE** parallel region starts
- **INSIDE** parallel region with explicit flush

### Commands to Run

```bash
cd /gpfs/projects/nct_352/nct01225/collatz/MNv3/

# Rebuild with new diagnostics
./build_v15_mpi_lean.sh

# Run sanity test again
sbatch test_v15mpi_1node_sanity.slurm

# Wait for completion
squeue -u nct01225

# Check stderr - should now see DEBUG messages
cat v15mpi_1node_sanity_*.err | tail -20
```

### What to Look For

**Expected stderr output:**
```
[VALIDATE] Self-test passed
[DEBUG] About to enter OpenMP parallel region
[DEBUG] omp_get_max_threads() = 112
[DEBUG] total_candidates = 66666664
[RANK 0] OpenMP started with 112 threads (max=112)
```

**If you see:**
- `[DEBUG] omp_get_max_threads() = 1` → Environment variable issue!
- No `[RANK 0] OpenMP started...` → Parallel region not executing!
- `[RANK 0] OpenMP started with 1 threads` → OpenMP running but only 1 thread!

### The Mystery

The diagnostic I added should have printed, but it didn't appear in stderr. This new version has explicit `fflush(stderr)` and prints BEFORE the parallel region to help us understand what's happening.
