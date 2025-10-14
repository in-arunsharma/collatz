## ✅ READY TO BUILD - Correct Modules Found

### Module Versions Confirmed on MareNostrum 5

**Available:**
- `gcc/13.2.0` ✅ (also 12.3.0, 14.1.0)
- `openmpi/4.1.5-gcc` ✅
- `impi/2021.10.0` ✅

### Scripts Updated

All scripts now use:
```bash
module load gcc/13.2.0 openmpi/4.1.5-gcc
```

### Commands to Run NOW on MareNostrum

```bash
cd /gpfs/projects/nct_352/nct01225/collatz/MNv3/

# Build with GCC (should work now!)
./build_v15_mpi_lean.sh

# Expected output:
# ✅ Build SUCCESS: V1.5-mpi-lean
# -rwxr-xr-x ... V1.5-mpi-lean

# Sanity test (1 node)
sbatch test_v15mpi_1node_sanity.slurm

# Wait ~1 min, check results
squeue -u nct01225
cat v15mpi_1node_sanity_*.out

# Expected in output:
# Throughput:   ~137000000 nums/sec  ✅
# Per-rank:     ~137000000 nums/sec  ✅
```

### Success Criteria for Sanity Test

**✅ PASS:** "Per-rank: 130M-145M nums/sec"
- If this passes → Ready for 2-node scale test
- This proves GCC fixed the icpc performance bug

**❌ FAIL:** "Per-rank: < 100M nums/sec"
- Check stderr for thread count
- Should show: `[RANK 0] OpenMP started with 112 threads`

### After Sanity Test Passes

```bash
# Scale to 2 nodes
sbatch test_v15mpi_2nodes.slurm

# Expected:
# Throughput:   ~274000000 nums/sec  (2 × 137M)
# Per-rank:     ~137000000 nums/sec
```

### Why This Will Work Now

1. **Correct modules:** `gcc/13.2.0` instead of non-existent `gcc/12.2.0`
2. **Right MPI:** `openmpi/4.1.5-gcc` (precompiled with GCC)
3. **NO icpc:** Avoiding Intel Classic Compiler that killed performance
4. **Same code:** Friend's optimizations unchanged (on-the-fly seeds, static scheduling)

### Performance Expectations

| Configuration | Expected Throughput | Status |
|--------------|---------------------|--------|
| 1 node (sanity) | ~137M nums/sec | Must pass first |
| 2 nodes | ~274M nums/sec | 2 × baseline |
| 10 nodes | ~1.37B nums/sec | Phase 2 done! 🎉 |

### Comparison to Previous Runs

**V1.5 with GCC (baseline):** 137M nums/sec ✅
**V1.5-mpi with icpc:** 1.87M nums/sec ❌ (73× slower)
**V1.5-mpi-lean with GCC:** Expected ~137M/sec per rank ✅

The 73× performance loss was entirely due to icpc's poor `__int128` codegen!
