# 🚀 DEPLOY TO MARENOSTRUM - RIGHT NOW!

## Your Login Info
```
ssh nct01225@glogin1.bsc.es
```

---

## 📦 Step 1: Package Code (30 seconds)

```bash
cd ~/Desktop/MN25
tar czf collatz_mn5.tar.gz marenostrum/
ls -lh collatz_mn5.tar.gz  # Should be ~180KB
```

---

## 📤 Step 2: Transfer to MareNostrum (1 minute)

```bash
scp collatz_mn5.tar.gz nct01225@glogin1.bsc.es:~/
```

**Expected:** Progress bar, completes in ~10-30 seconds

---

## 🔐 Step 3: SSH to MareNostrum

```bash
ssh nct01225@glogin1.bsc.es
```

**You're now on the cluster!**

---

## 📂 Step 4: Setup in SCRATCH Filesystem (1 minute)

**IMPORTANT:** Don't run from `/gpfs/home` - use `/gpfs/scratch` instead!

```bash
# Check your group name (e.g., bsc01, bsc32, etc.)
groups

# Create workspace in SCRATCH (REPLACE bscXX with your group!)
mkdir -p /gpfs/scratch/bscXX/nct01225/collatz
cd /gpfs/scratch/bscXX/nct01225/collatz

# Extract code here
tar xzf ~/collatz_mn5.tar.gz

# Navigate to marenostrum folder
cd marenostrum/

# Make scripts executable
chmod +x *.sh

# Verify files
ls -lh

# Optional: Check disk quota
bsc_quota
```

**Expected:** See all your files (16 files, ~180KB total)

**Why scratch?** Optimized for job I/O, much faster than home directory!

---

## 🔨 Step 5: Build on MareNostrum (1 minute)

```bash
./build_phase1_gpp.sh
```

**Expected output:**
```
Loading Intel OneAPI modules...
Loading OpenMPI modules...
Compiling...
✅ Build successful: collatz_mpi_gpp
```

**If it fails:** Check error messages, might need to load modules manually

---

## 🎯 Step 6: Submit First Validation Job (INSTANT)

**Current configuration (already set):**
- **Nodes:** 2 (conservative, only 1.3% of pool)
- **Test range:** 100 million numbers (ultra-fast validation)
- **Walltime:** 15 minutes (will finish in ~1 second)

```bash
sbatch slurm_phase1_gpp.slurm
```

**Expected:**
```
Submitted batch job 123456
```

**WRITE DOWN THE JOB ID!** ___________

---

## 👀 Step 7: Monitor Job

### Check queue:
```bash
squeue -u $USER
```

**States:**
- `PD` = Pending (waiting in queue)
- `R` = Running
- `CG` = Completing

### Watch output in real-time:
```bash
# Replace 123456 with your actual job ID
tail -f collatz_phase1_123456.out
```

**Press Ctrl+C to stop watching**

### Check for errors:
```bash
tail -f collatz_phase1_123456.err
```

---

## ✅ Step 8: Check Results (after job completes)

```bash
# List output files
ls -lh *.json *.out *.err

# View main output
cat collatz_phase1_*.out

# Check for success
grep "SUCCESS" collatz_phase1_*.out
```

**Success indicators:**
- ✅ Output file says "SUCCESS!"
- ✅ No errors in .err file
- ✅ JSON files created
- ✅ Throughput reported

---

## 📊 Expected Results (2 nodes, 100M numbers)

**With placeholder computation:**
- Throughput: ~0.5-2 M numbers/sec
- Runtime: ~50-200 seconds
- Just validating infrastructure

**What validates:**
- ✅ SLURM job submission works
- ✅ MPI ranks spawn on 2 different nodes
- ✅ OpenMP threads detected (112 per node)
- ✅ Work distribution functions
- ✅ Results aggregation works
- ✅ Output files generated

---

## 🔄 Step 9: Scale Up (after validation passes)

### Option A: More nodes (3 → 5 → 10)
```bash
nano slurm_phase1_gpp.slurm

# Change:
#SBATCH --nodes=5
#SBATCH --ntasks=5
END_SEED=1000000000  # 1 billion

# Submit again
sbatch slurm_phase1_gpp.slurm
```

### Option B: Longer run
```bash
nano slurm_phase1_gpp.slurm

# Change:
END_SEED=10000000000  # 10 billion
#SBATCH --time=01:00:00  # 1 hour

# Submit again
sbatch slurm_phase1_gpp.slurm
```

---

## 🚨 Troubleshooting

### Build fails:
```bash
# Load modules manually
module load intel/2023.2
module load openmpi/4.1.5
./build_phase1_gpp.sh
```

### Job stays pending (PD) for long:
```bash
# Check queue status
squeue -p gpp  # See all jobs on GPP partition

# If very busy, try 1 node:
nano slurm_phase1_gpp.slurm
# Change to: --nodes=1 --ntasks=1
```

### Job fails immediately:
```bash
# Check error file
cat collatz_phase1_*.err

# Check job details
scontrol show job <JOBID>
```

### Can't SSH:
- Check VPN if required
- Verify username: `nct01225`
- Verify hostname: `glogin1.bsc.es`

---

## 📝 Recommended Progression

**Job 1: Infrastructure Validation** (NOW!)
- 2 nodes, 100M numbers, 15 min walltime
- Goal: Verify everything works
- Expected: ~1-2 minute runtime

**Job 2: Small Production** (if Job 1 succeeds)
- 3 nodes, 1B numbers, 30 min walltime
- Goal: Get real performance baseline
- Expected: ~2-5 minute runtime

**Job 3: Scale Up** (if Job 2 succeeds)
- 5 nodes, 10B numbers, 1 hour walltime
- Goal: Target 1B+ nums/sec throughput
- Expected: ~10-15 minute runtime

**Job 4: Full Run** (if all validate)
- 10 nodes, 100B numbers, 2 hour walltime
- Goal: Maximum throughput on GPP only
- Expected: ~40-60 minute runtime

---

## 🎯 Success Criteria for FIRST JOB

You'll know it worked if:
1. ✅ Job completes without errors
2. ✅ Both nodes are detected as GPP type (or treated as such)
3. ✅ 2 MPI ranks × 112 threads = 224 total threads reported
4. ✅ All ranks complete computation
5. ✅ Results aggregated successfully
6. ✅ JSON files generated

**Throughput doesn't matter yet** - you're validating infrastructure!

---

## ⏱️ Timeline for First Deployment

```
Now       → +1 min     Package & transfer
+1 min    → +2 min     SSH & extract
+2 min    → +3 min     Build
+3 min    → +3 min     Submit job
+3 min    → +5 min     Job runs (queue + execution)
+5 min    → +10 min    Analyze results
────────────────────────────────────────
~10 min   FIRST VALIDATION COMPLETE! ✅
```

---

## 🚀 GO NOW!

```bash
# 1. Package code
cd ~/Desktop/MN25
tar czf collatz_mn5.tar.gz marenostrum/

# 2. Transfer
scp collatz_mn5.tar.gz nct01225@glogin1.bsc.es:~/

# 3. SSH
ssh nct01225@glogin1.bsc.es

# 4. Setup in SCRATCH (IMPORTANT!)
groups  # Note your group name
mkdir -p /gpfs/scratch/<YOUR_GROUP>/nct01225/collatz
cd /gpfs/scratch/<YOUR_GROUP>/nct01225/collatz
tar xzf ~/collatz_mn5.tar.gz
cd marenostrum/

# 5. Build and submit
./build_phase1_gpp.sh
sbatch slurm_phase1_gpp.slurm
```

**You're 10 minutes away from running on a supercomputer! 💪**

**See `MARENOSTRUM_FILESYSTEM.md` for detailed filesystem guide!**
