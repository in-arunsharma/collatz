# V1.6b Testing Plan - October 14, 2025

## Quick Summary
Test the **bug-fixed** V1.6b-mpi-openmp.cpp to verify:
1. Single-node matches Phase 1 baseline (137M+ nums/sec)
2. Multi-node scales linearly (not 54× regression!)
3. L2 cache optimization (2^19 table) gives 10-20% boost

---

## Step 1: Deploy V1.6b to MareNostrum

```bash
# From local machine
cd /home/aruns/Desktop/MN25/MNv3

scp V1.6b-mpi-openmp.cpp V1.6_BUGFIX_REPORT.md \
    nct01225@glogin1.bsc.es:/gpfs/projects/nct_352/nct01225/collatz/MNv3/
```

---

## Step 2: Build on MareNostrum

```bash
# SSH to MareNostrum
ssh nct01225@glogin1.bsc.es
cd /gpfs/projects/nct_352/nct01225/collatz/MNv3/

# Load modules
module load gcc openmpi

# Build V1.6b
mpicxx -O3 -march=native -fopenmp -o V1.6b-mpi-openmp V1.6b-mpi-openmp.cpp

# Verify binary
ls -lh V1.6b-mpi-openmp
```

---

## Step 3: Single-Node Validation (CRITICAL)

**Purpose**: Verify V1.6b matches V1.5 baseline (137M nums/sec)

Create `test_v16b_1node.slurm`:
```bash
cat > test_v16b_1node.slurm << 'EOF'
#!/bin/bash
#SBATCH --job-name=v16b_1node
#SBATCH --output=collatz_v16b_1node_%j.out
#SBATCH --error=collatz_v16b_1node_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=112
#SBATCH --time=00:10:00
#SBATCH --qos=debug
#SBATCH --partition=gpp

module load gcc openmpi

export OMP_NUM_THREADS=112
export OMP_PROC_BIND=spread
export OMP_PLACES=cores

echo "=== V1.6b Single-Node Test ==="
echo "Date: $(date)"
echo "Node: $(hostname)"
echo "Cores: 112"
echo ""

mpirun -np 1 ./V1.6b-mpi-openmp 0 100000000 --tag v16b_1node

echo ""
echo "=== Expected: ~150M nums/sec (improved with L2 cache) ==="
echo "=== Baseline: 137M nums/sec (Phase 1 result) ==="
EOF

sbatch test_v16b_1node.slurm
```

**Success Criteria**:
- ✅ Throughput: 137-165M nums/sec (137M baseline, up to +20% from L2 optimization)
- ✅ Tested: ~33.3M numbers (after mod-6 filtering)
- ✅ Time: ~220ms compute (excluding precompute)

**If this fails**: V1.6b has a new bug, don't proceed to multi-node

---

## Step 4: Two-Node Scaling Test (THE BIG TEST)

**Purpose**: Verify bug fix - should get 2× single-node, NOT 0.018× like V1.6!

Create `test_v16b_2nodes.slurm`:
```bash
cat > test_v16b_2nodes.slurm << 'EOF'
#!/bin/bash
#SBATCH --job-name=v16b_2nodes
#SBATCH --output=collatz_v16b_2nodes_%j.out
#SBATCH --error=collatz_v16b_2nodes_%j.err
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=112
#SBATCH --time=00:10:00
#SBATCH --qos=debug
#SBATCH --partition=gpp

module load gcc openmpi

export OMP_NUM_THREADS=112
export OMP_PROC_BIND=spread
export OMP_PLACES=cores

echo "=== V1.6b Two-Node Scaling Test ==="
echo "Date: $(date)"
echo "Nodes: 2"
echo "Total Cores: 224"
echo ""

mpirun -np 2 ./V1.6b-mpi-openmp 0 100000000 --tag v16b_2nodes

echo ""
echo "=== Expected: ~300M nums/sec (2× single-node) ==="
echo "=== V1.6 broken result: 2.5M nums/sec ==="
EOF

sbatch test_v16b_2nodes.slurm
```

**Success Criteria**:
- ✅ Throughput: 270-330M nums/sec (near-linear scaling)
- ✅ Tested: ~33.3M numbers total (same as single-node, split between ranks)
- ✅ Time: ~110ms compute (half of single-node)
- ✅ Speedup: 1.8-2.0× vs single-node (90-100% scaling efficiency)

**Expected Issues (if any)**:
- Slightly less than 2× due to MPI overhead
- 90-95% efficiency is excellent for this workload

---

## Step 5: Ten-Node Scaling (FULL POWER)

**Purpose**: Test maximum GPP scaling for hackathon

Create `test_v16b_10nodes.slurm`:
```bash
cat > test_v16b_10nodes.slurm << 'EOF'
#!/bin/bash
#SBATCH --job-name=v16b_10nodes
#SBATCH --output=collatz_v16b_10nodes_%j.out
#SBATCH --error=collatz_v16b_10nodes_%j.err
#SBATCH --nodes=10
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=112
#SBATCH --time=00:15:00
#SBATCH --qos=acc_debug
#SBATCH --partition=gpp

module load gcc openmpi

export OMP_NUM_THREADS=112
export OMP_PROC_BIND=spread
export OMP_PLACES=cores

echo "=== V1.6b Ten-Node Scaling Test ==="
echo "Date: $(date)"
echo "Nodes: 10"
echo "Total Cores: 1,120"
echo ""

mpirun -np 10 ./V1.6b-mpi-openmp 0 100000000 --tag v16b_10nodes

echo ""
echo "=== Expected: ~1.5B nums/sec (10× single-node) ==="
EOF

sbatch test_v16b_10nodes.slurm
```

**Success Criteria**:
- ✅ Throughput: 1.3-1.6B nums/sec (near-linear scaling to 10 nodes)
- ✅ Speedup: 8-10× vs single-node (80-100% scaling efficiency)
- ✅ Per-node: ~150M nums/sec (consistent with single-node)

---

## Step 6: Monitor Jobs

```bash
# Check queue
squeue -u nct01225

# Watch specific job (replace JOBID)
watch -n 2 "squeue -j JOBID"

# Check results as they complete
ls -lht collatz_v16b_*.out collatz_v16b_*.err

# Quick view
tail -50 collatz_v16b_1node_*.out
tail -50 collatz_v16b_2nodes_*.out
tail -50 collatz_v16b_10nodes_*.out
```

---

## Step 7: Compare Results

**What to look for**:

### ✅ **SUCCESS Pattern**:
```
1 node:  150M nums/sec  (baseline)
2 nodes: 300M nums/sec  (2.0× scaling)
10 nodes: 1.5B nums/sec (10.0× scaling)
```

### ❌ **FAILURE Pattern** (like V1.6):
```
1 node:  150M nums/sec
2 nodes: 2.5M nums/sec  (0.017× - WRONG!)
```

### Key Metrics:
1. **Throughput**: Should scale linearly (within 90%)
2. **Numbers tested**: Should be ~33.3M for all tests (same workload)
3. **Compute time**: Should decrease proportionally (2 nodes = half time)
4. **Edge cases**: Should all be 0 (overflow, fuse, cycles)

---

## Step 8: Analysis Commands

```bash
# Extract key metrics
grep "Throughput:" collatz_v16b_*.out
grep "Tested:" collatz_v16b_*.out
grep "Time:" collatz_v16b_*.out
grep "Speedup:" collatz_v16b_*.out

# Calculate scaling efficiency
# (2-node throughput / 1-node throughput) should be ~2.0
# (10-node throughput / 1-node throughput) should be ~10.0

# Check for errors
cat collatz_v16b_*.err
```

---

## Step 9: Decision Tree

### If 1-node test FAILS (< 120M nums/sec):
- Problem: V1.6b introduced new bug or L2 optimization backfired
- Action: Compare with V1.5-openmp (known-good baseline)
- Fix: May need to revert to 2^20 table or debug V1.6b

### If 1-node PASSES, 2-node FAILS (< 250M nums/sec):
- Problem: MPI work distribution still broken
- Action: Review V1.6b work splitting logic
- Fix: Debug the global list construction

### If 1-node and 2-node PASS, 10-node FAILS:
- Problem: Network bottleneck or SLURM issue
- Action: Test 4 nodes, 6 nodes to find scaling limit
- Fix: May need different MPI configuration

### If ALL PASS:
- **SUCCESS!** Phase 2 complete
- Action: Document results, move to Phase 3 (CUDA)
- Next: Start GPU prototype on ACC partition

---

## Expected Timeline

- **Step 1-2** (Deploy + Build): 5 minutes
- **Step 3** (1-node test): 10 minutes (queue + run)
- **Step 4** (2-node test): 10 minutes
- **Step 5** (10-node test): 15 minutes
- **Total**: ~40 minutes to full validation

---

## Quick Commands Summary

```bash
# === ON LOCAL MACHINE ===
cd /home/aruns/Desktop/MN25/MNv3
scp V1.6b-mpi-openmp.cpp V1.6_BUGFIX_REPORT.md \
    nct01225@glogin1.bsc.es:/gpfs/projects/nct_352/nct01225/collatz/MNv3/

# === ON MARENOSTRUM ===
ssh nct01225@glogin1.bsc.es
cd /gpfs/projects/nct_352/nct01225/collatz/MNv3/

# Build
module load gcc openmpi
mpicxx -O3 -march=native -fopenmp -o V1.6b-mpi-openmp V1.6b-mpi-openmp.cpp

# Create test jobs (copy from Step 3, 4, 5 above)
# Then submit:
sbatch test_v16b_1node.slurm    # Wait for this to finish first!
sbatch test_v16b_2nodes.slurm   # Only if 1-node passes
sbatch test_v16b_10nodes.slurm  # Only if 2-node passes

# Monitor
squeue -u nct01225
tail -f collatz_v16b_1node_*.out

# Results
grep "Throughput:" collatz_v16b_*.out
```

---

## Critical: Test in Sequence!

**DO NOT** submit all jobs at once. Test in order:
1. ✅ 1-node MUST pass before testing 2-node
2. ✅ 2-node MUST pass before testing 10-node

This saves your hackathon time allocation if there's still a bug.

---

## Success Criteria Summary

| Test | Expected Throughput | Scaling Efficiency | Critical? |
|------|-------------------|-------------------|-----------|
| 1 node | 137-165M nums/sec | N/A (baseline) | YES - Blocks all other tests |
| 2 nodes | 270-330M nums/sec | 90-100% | YES - Proves bug fix worked |
| 10 nodes | 1.3-1.6B nums/sec | 80-100% | GOAL - Maximum GPP throughput |

---

## After Testing: Next Steps

### If Phase 2 succeeds:
1. Document final GPP results
2. Update Phase 2 status to ✅
3. Begin Phase 3: CUDA GPU prototype (V1.7-cuda)
4. Target: Single Hopper GPU on ACC partition

### If Phase 2 needs iteration:
1. Analyze failure mode from test results
2. Create V1.6c with additional fixes
3. Re-test with same methodology

---

## Questions to Answer from Tests

1. **Does the bug fix work?** (2-node should be 2× single-node)
2. **Does L2 cache help?** (Compare to V1.5's 137M baseline)
3. **What's the scaling limit?** (10 nodes should be near-linear)
4. **Are we ready for GPUs?** (If GPP maxes at 1.5B, GPU should target 5-50B)

---

**Good luck! The bug fix looks solid - expecting Phase 2 success! 🚀**
