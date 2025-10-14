## 🔍 URGENT: Performance Still Broken (2.5M instead of 137M)

### Current Status
- ✅ GCC build succeeded (56KB binary)
- ❌ Performance: 2.5M nums/sec (should be 137M) - **54× too slow**
- ⚠️  This is DIFFERENT from icpc (was 1.87M), so new issue!

### Most Likely Cause
**OpenMP not linked or not running with 112 threads**

### IMMEDIATE Action on MareNostrum

```bash
cd /gpfs/projects/nct_352/nct01225/collatz/MNv3/

# Run diagnostic script
./diagnose_performance.sh
```

### What to Look For

**Scenario 1: OpenMP NOT linked**
```
⚠️  NO OpenMP library linked!
⚠️  NO OpenMP symbols found!
```
**Solution:** Rebuild with correct OpenMP flags

**Scenario 2: OpenMP linked but only 1 thread**
```
stderr shows: [RANK 0] OpenMP started with 1 threads (max=1)
```
**Solution:** Environment variable issue

**Scenario 3: OpenMP working with 112 threads**
```
stderr shows: [RANK 0] OpenMP started with 112 threads (max=112)
```
**Solution:** Something else is slow (need deeper investigation)

### If OpenMP Not Linked (Scenario 1)

The build used `-fopenmp` but maybe it needs `-lgomp` explicitly:

```bash
# Try explicit OpenMP library linking
mpicxx -O3 -march=native -flto -fno-exceptions -fno-rtti -funroll-loops \
       -DNDEBUG -fopenmp -lgomp \
       -o V1.5-mpi-lean V1.5-mpi-lean.cpp

# Or try verbose build to see what's happening
./build_v15_verbose.sh
```

### Key Numbers to Remember

| Metric | V1.5 Baseline | Current | Should Be |
|--------|---------------|---------|-----------|
| Throughput | 137M/sec | 2.5M/sec | 137M/sec |
| Time for 66M | ~500ms | 26,656ms | ~500ms |
| Per-number | 0.0073 μs | 0.4 μs | 0.0073 μs |

The 54× slowdown suggests single-threaded execution (112 threads → 1 thread = 112× theoretical, but with overhead ~50-60×).

### Next Steps

1. Run `./diagnose_performance.sh` 
2. Check stderr: `cat v15mpi_1node_sanity_*.err`
3. Report what the diagnostic script shows
4. We'll fix based on which scenario it is
