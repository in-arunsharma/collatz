# 🚀 HACKATHON MORNING - DEPLOYMENT GUIDE

## 📍 Your MareNostrum Login
```
ssh nct01225@glogin1.bsc.es
```

## 🎯 Resource Allocation
**Shared pool for all participants:**
- 150 GPP nodes (Intel Sapphire Rapids, 112 cores each)
- 25 ACC nodes (Intel Sapphire Rapids + 4× NVIDIA H100)

**Our strategy:** Start with 3 nodes, scale up as needed

---

## ✅ SELF-CONTAINED SETUP

**Everything is now in the `marenostrum/` folder!**
- Core algorithm: `collatz_core_v14b.cpp` (V1.4b-openmp)
- MPI infrastructure: coordinator, workers, main program
- Build scripts: local + MareNostrum
- SLURM job script (pre-configured for 3 nodes)
- Documentation: all guides

**No external dependencies!** Just copy one folder to MareNostrum.

---

## 🏃 QUICK START (3 commands to test locally)

```bash
cd ~/Desktop/MN25/marenostrum/

# 1. Integrate V1.4b core
./integrate_v14b.sh

# 2. Build
./build_phase1_gpp.sh

# 3. Test locally
./test_local_mpi.sh 2
```

**Expected:** ✓ Test PASSED

---

## 📦 DEPLOY TO MARENOSTRUM

### Step 1: Package (everything in one folder!)
```bash
cd ~/Desktop/MN25
tar czf collatz_mn5.tar.gz marenostrum/
```

### Step 2: Transfer
```bash
scp collatz_mn5.tar.gz nct01225@glogin1.bsc.es:~/
```

### Step 3: SSH and extract
```bash
ssh nct01225@glogin1.bsc.es

# On MareNostrum:
tar xzf collatz_mn5.tar.gz
cd marenostrum/
ls -lh  # Verify all files are there
```

### Step 4: Build on MareNostrum
```bash
./build_phase1_gpp.sh
```

**Expected output:**
```
Loading Intel OneAPI modules...
Loading OpenMPI modules...
Compiling for MareNostrum 5...
Build successful: collatz_mpi_gpp
```

### Step 5: Submit test job (3 nodes)
```bash
sbatch slurm_phase1_gpp.slurm
```

**Expected:**
```
Submitted batch job 123456
```

### Step 6: Monitor
```bash
# Check queue
squeue -u $USER

# Watch output
tail -f collatz_phase1_*.out

# Check for errors
tail -f collatz_phase1_*.err
```

---

## 📊 INITIAL TEST (3 nodes, 1 billion numbers)

**Configuration (already set in slurm script):**
- Nodes: 3
- Seed range: 0 to 1,000,000,000
- Walltime: 30 minutes (will finish in ~2 seconds)

**Expected performance:**
- Throughput: ~750M numbers/sec
- Runtime: ~1.5 seconds
- Total cores: 336 (3 nodes × 112 cores)

**Success criteria:**
- ✅ Job completes without errors
- ✅ Throughput > 500M/sec
- ✅ All 3 ranks produce JSON output
- ✅ Cycle seeds detected (~50-75)
- ✅ Overflow seeds detected (~10-15)

---

## 🔧 SCALING UP (after initial test passes)

### Option A: More GPP nodes (if available)
Edit `slurm_phase1_gpp.slurm`:
```bash
#SBATCH --nodes=5        # 5 nodes = 1.25B nums/sec
#SBATCH --ntasks=5
END_SEED=10000000000     # 10 billion for production
```

### Option B: ACC nodes (Phase 2 - CUDA)
- Wait for Phase 2 implementation
- ACC nodes have 4× H100 GPUs each
- Expected: 500M-1B nums/sec per ACC node

### Conservative Resource Usage
```
Initial test:    3 GPP nodes (2% of pool)
Production:      5-10 GPP nodes (3-7% of pool)
Be courteous:    Monitor queue, don't hog resources
```

---

## 📋 VERIFICATION CHECKLIST

After job completes:

```bash
# 1. Check output files exist
ls -lh collatz_phase1_*.out
ls -lh *.json

# 2. Verify all ranks completed
grep "Rank.*completed" collatz_phase1_*.out

# 3. Check throughput
cat *_global.json | grep throughput

# 4. View aggregated results
cat *_global.json
```

**Example successful output:**
```json
{
  "global_results": {
    "total_numbers": 1000000000,
    "num_workers": 3,
    "elapsed_seconds": 1.35,
    "throughput_per_sec": 740740740,
    "max_steps_global": 525,
    "cycles_found_total": 57,
    "overflows_found_total": 12
  }
}
```

---

## 🚨 TROUBLESHOOTING

### Build fails on MareNostrum
```bash
# Check modules
module list

# Expected:
# intel/2023.2
# openmpi/4.1.5

# If missing, load manually:
module load intel/2023.2
module load openmpi/4.1.5
./build_phase1_gpp.sh
```

### Job stays pending (PD)
```bash
# Check queue
squeue -p gpp  # See all GPP partition jobs

# If very busy, try:
#SBATCH --nodes=2  # Fewer nodes
#SBATCH --time=00:15:00  # Shorter time
```

### Low performance
```bash
# Check error file
cat collatz_phase1_*.err

# Verify OpenMP threads
grep OMP_NUM_THREADS collatz_phase1_*.out
# Should show: OMP_NUM_THREADS: 112
```

### No output files
```bash
# Check if job ran
cat collatz_phase1_*.out

# Look for errors
cat collatz_phase1_*.err

# Verify executable exists
ls -lh collatz_mpi_gpp
```

---

## ⏱️ HACKATHON TIMELINE

**09:00 - Local Testing (15 min)**
- Run integration script
- Build locally
- Test with 2 MPI ranks
- Package for transfer

**09:15 - MareNostrum Deployment (15 min)**
- Transfer tar.gz
- SSH to cluster
- Extract and build
- Submit test job (3 nodes, 1B numbers)

**09:30 - Validation (15 min)**
- Monitor job completion
- Check results (JSON files)
- Verify throughput > 500M/sec
- Analyze cycle/overflow detection

**09:45 - Scale Up (30 min)**
- Increase to 5-10 nodes if resources available
- Production run with 10B+ numbers
- Achieve 1B+ numbers/sec

**10:15 - Phase 1 Complete! ✅**
- Move to Phase 2 (ACC nodes with CUDA)

---

## 📖 DOCUMENTATION

All docs are in `marenostrum/` folder:

- **HACKATHON_SETUP.md** ← You are here!
- **QUICKSTART.md** - Quick reference
- **EXECUTION_GUIDE.md** - Detailed walkthrough
- **INTEGRATION_GUIDE.md** - V1.4b integration details
- **CHECKLIST.txt** - Printable checklist

---

## 🎯 SUCCESS METRICS

**Phase 1 Goals:**
- ✅ Deploy on 3+ GPP nodes
- ✅ Achieve 500M+ numbers/sec
- ✅ Validate cycle/overflow detection
- ✅ Stable multi-node execution
- ✅ Clean JSON output

**Stretch Goals:**
- 🚀 Scale to 10 GPP nodes (2.5B nums/sec)
- 🚀 Process 100B+ numbers total
- 🚀 Find rare cycle seeds
- 🚀 Move to ACC nodes (Phase 2)

---

**LET'S GO! START WITH THE 3 QUICK START COMMANDS! 💪**

```bash
cd ~/Desktop/MN25/marenostrum/
./integrate_v14b.sh
./build_phase1_gpp.sh
./test_local_mpi.sh 2
```
