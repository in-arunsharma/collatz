## 🚀 V1.5-mpi-lean: Quick Commands (GCC Build)

### CRITICAL ISSUE: Still Getting 2.5M/sec with GCC!

**Status:** Build with GCC succeeded, but performance still 54× too slow
**Likely cause:** OpenMP threads not spawning OR memory access issue

### Check stderr NOW:
```bash
cat v15mpi_1node_sanity_*.err
```

Look for:
- `[RANK 0] OpenMP started with N threads` - should be 112
- Any warnings about OpenMP or threading
- Any MPI binding messages

### On MareNostrum - Diagnostic Commands

```bash
cd /gpfs/projects/nct_352/nct01225/collatz/MNv3/

# Check stderr for thread diagnostics
cat v15mpi_1node_sanity_*.err | grep -i "thread\|omp\|rank"

# Step 2: 1-node sanity test (CRITICAL!)
sbatch test_v15mpi_1node_sanity.slurm
squeue -u nct01225
# Wait ~1 min
cat v15mpi_1node_sanity_*.out

# Check output:
# ✅ PASS: "Per-rank: ~137M nums/sec" → Go to Step 3
# ❌ FAIL: "Per-rank: < 100M nums/sec" → Check stderr for binding issues

# Step 3: Scale to 2 nodes (if sanity passed)
sbatch test_v15mpi_2nodes.slurm
squeue -u nct01225
cat v15mpi_2nodes_*.out

# Expected: "Throughput: ~274M nums/sec"
```

### What Changed
- **Compiler:** icpc → GCC (via `I_MPI_CXX=g++`)
- **Pinning:** Added `I_MPI_PIN_DOMAIN=omp` + `KMP_AFFINITY`
- **Code:** ZERO changes (same friend's optimizations)

### Success Criteria
| Test | Expected | Status |
|------|----------|--------|
| 1-node sanity | ~137M/sec per-rank | Must pass first! |
| 2-node scale | ~274M/sec total | 2 × 137M |
| 10-node scale | ~1.37B/sec total | Phase 2 COMPLETE! |

### If Sanity Test Still Fails
Check stderr for:
```
[RANK 0] OpenMP started with N threads
```
Should show N=112 threads per rank.

If N=1 or wrong number → Binding issue, check `I_MPI_DEBUG=4` output.

### After Success
1. Set `I_MPI_DEBUG=0` in SLURM scripts (reduces stderr spam)
2. Scale to 10 nodes
3. Move to Phase 3 (CUDA) 🎉
