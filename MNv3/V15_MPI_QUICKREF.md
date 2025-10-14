## V1.5-mpi-lean: Quick Reference

### 🎯 What Changed (Friend's Optimizations)

**1. On-the-fly seed generation** (NO vector)
- OLD: Built `std::vector<uint128_t> numbers_to_test` (~528MB for 33M numbers)
- NEW: `n = start + ((idx>>1)*6) + ((idx&1)?4:0)` - computed in loop

**2. Static scheduling** (uniform work)
- OLD: `schedule(dynamic, 512)` - high overhead
- NEW: `schedule(static, 8192)` - zero overhead

**3. MPI work distribution**
- Simple offset splitting per rank
- Each rank: `my_offset = global_start + rank * block`

**4. Optional cold queue disable**
- NEW: `--cold off` flag for pure throughput runs

### 📋 Commands on MareNostrum

**Build:**
```bash
cd /gpfs/projects/nct_352/nct01225/collatz/MNv3/
./build_v15_mpi_lean.sh
```

**Test 2 nodes:**
```bash
sbatch test_v15mpi_2nodes.slurm
squeue -u nct01225
# Wait ~1min
cat v15mpi_2nodes_*.out
```

**Expected Results:**
- Per-rank: ~137M nums/sec (same as V1.5 baseline)
- Total: ~274M nums/sec (2 × 137M)
- Wall time: ~120-130ms

### 🔧 Binding Rules (CRITICAL!)

**Quote from friend:**
> "If your per-rank 'Per-rank' number isn't ≈137M on a node (± a bit), your binding is off. Fix binding before touching code."

**Correct binding (already in SLURM script):**
- `--ntasks-per-node=1` - 1 MPI rank per node
- `--cpus-per-task=112` - All cores to OpenMP
- `export OMP_PLACES=cores`
- `export OMP_PROC_BIND=spread`
- `srun --cpu-bind=cores` - Bind to cores not threads

### 🚀 Scaling to 10 Nodes

If 2-node test shows ~270M+ nums/sec:

```bash
cp test_v15mpi_2nodes.slurm test_v15mpi_10nodes.slurm
# Edit: --nodes=10
sbatch test_v15mpi_10nodes.slurm
```

**Expected:** ~1.37B nums/sec (10 × 137M)

### ✅ Success Criteria

**2-node test:**
- ✅ Per-rank ≥ 130M nums/sec (95% of baseline)
- ✅ Total ≥ 260M nums/sec (95% scaling efficiency)
- ❌ If per-rank < 100M → Binding problem

**10-node test:**
- ✅ Per-rank ≥ 130M nums/sec
- ✅ Total ≥ 1.3B nums/sec
- ✅ If successful → Phase 2 COMPLETE!

### 🔑 Key Files

- `V1.5-mpi-lean.cpp` - Source (friend's optimizations)
- `build_v15_mpi_lean.sh` - Build with aggressive flags
- `test_v15mpi_2nodes.slurm` - 2-node test job
- `V1.5-openmp.cpp` - Original (UNCHANGED - still works)

### 🎓 What Friend Fixed

**Problem 1:** V1.6/V1.7/V1.8 built huge vectors → cache thrashing
**Solution:** Compute seeds on-the-fly (zero allocation)

**Problem 2:** Dynamic scheduling overhead
**Solution:** Static scheduling with large chunks (8192)

**Problem 3:** Complex work distribution
**Solution:** Simple offset arithmetic per rank

**All V1.5 kernels UNCHANGED** - proven fast code stays fast!
