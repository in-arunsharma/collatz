## 🐛 ROOT CAUSE FOUND: Intel Classic Compiler (icpc)

### The Problem

**Your fast V1.5:** 137M nums/sec with GCC
**Your V1.5-mpi:** 1.87M nums/sec with icpc
**Performance loss:** 73× slower!

### Why icpc Failed

Intel Classic Compiler (icpc) has **terrible codegen for `__int128` arithmetic**:
- Lowers 128-bit ops to helper function calls
- Misses critical inlining opportunities
- No vectorization-friendly transforms
- Deprecated by Intel (EOL in 2023)

**GCC/Clang:** Generate fast inline SIMD code for `__int128`
**icpc:** Generates slow library calls

### The Fix

**Use GCC as the MPI backend compiler:**

```bash
module purge
module load gcc/12.2.0 impi/2021.10.0
export I_MPI_CXX=g++  # Force GCC backend (NOT icpc!)
mpicxx -O3 -march=native -flto -fopenmp -fno-exceptions -fno-rtti -funroll-loops -DNDEBUG \
       -o V1.5-mpi-lean V1.5-mpi-lean.cpp
```

### Enhanced Pinning (Intel MPI)

Add these environment variables before `srun`:

```bash
# OpenMP settings
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PLACES=cores
export OMP_PROC_BIND=spread
export KMP_AFFINITY=granularity=fine,scatter
export KMP_BLOCKTIME=0

# Intel MPI pinning (CRITICAL!)
export I_MPI_PIN=1
export I_MPI_PIN_DOMAIN=omp       # Give whole OMP team to rank
export I_MPI_DEBUG=4              # Verify binding (first run only)
```

### Testing Strategy

**Step 1: Sanity Check (1 node)**
```bash
sbatch test_v15mpi_1node_sanity.slurm
```
**Expected:** ~137M nums/sec (same as V1.5 baseline)
**If fails:** Compiler or binding issue - fix before scaling

**Step 2: Scale to 2 Nodes**
```bash
sbatch test_v15mpi_2nodes.slurm
```
**Expected:** ~274M nums/sec (2 × 137M)

**Step 3: Scale to 10 Nodes**
```bash
# Edit test_v15mpi_2nodes.slurm: --nodes=10
sbatch test_v15mpi_10nodes.slurm
```
**Expected:** ~1.37B nums/sec (10 × 137M)

### Files Updated

✅ `build_v15_mpi_lean.sh` - Now uses GCC via `I_MPI_CXX=g++`
✅ `test_v15mpi_2nodes.slurm` - Enhanced pinning settings
✅ `test_v15mpi_1node_sanity.slurm` - NEW: Sanity test before scaling

### What Stays the Same

**All computation kernels are UNCHANGED from V1.5:**
- Same `compute_collatz_readonly()` hot path
- Same memo table structure
- Same Brent cycle detection
- Same on-the-fly seed generation (friend's optimization)
- Same static scheduling

**Only change:** GCC compiler instead of icpc

### Commands on MareNostrum

```bash
ssh nct01225@glogin2.bsc.es
cd /gpfs/projects/nct_352/nct01225/collatz/MNv3/

# Upload new files
# (from local: scp build_v15_mpi_lean.sh test_v15mpi_*.slurm nct01225@glogin2.bsc.es:...)

# Rebuild with GCC
./build_v15_mpi_lean.sh

# Sanity test (1 node) - MUST pass before scaling!
sbatch test_v15mpi_1node_sanity.slurm
squeue -u nct01225
cat v15mpi_1node_sanity_*.out

# If sanity passes (~137M/sec), scale to 2 nodes
sbatch test_v15mpi_2nodes.slurm
```

### Success Criteria

**1-node sanity:** Per-rank ≥ 130M nums/sec ✅
**2-node scale:** Total ≥ 260M nums/sec ✅
**10-node scale:** Total ≥ 1.3B nums/sec ✅ → Phase 2 COMPLETE!

### Why This Will Work

1. **GCC codegen:** Fast `__int128` inline code (proven in V1.5)
2. **Enhanced pinning:** `I_MPI_PIN_DOMAIN=omp` prevents thread migration
3. **Friend's optimizations:** On-the-fly seeds, static scheduling (no changes needed)
4. **Proven kernel:** Same hot path that achieved 137M/sec

**Quote from expert:**
> "If your per-rank number isn't ≈137M on a node (± a bit), your binding is off. Fix binding before touching code."

With GCC + proper pinning, you should see ~137M per rank immediately.
