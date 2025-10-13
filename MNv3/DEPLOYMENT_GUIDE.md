# MareNostrum 5 Deployment Quick Guide

## 🎯 Strategy: Start with GPP, then Scale

### Why GPP First?
1. **Faster deployment** (30 min vs 2-3 hours for GPU)
2. **Easier debugging** (CPU is familiar territory)
3. **Guaranteed resources** (less queue competition)
4. **Good performance** (100-200M nums/sec with 112 cores)
5. **Foundation for MPI** (test before GPU phase)

### Performance Targets
```
Sequential V1.4b:    2.2M nums/sec         (baseline)
Phase 1 (1 GPP node):  150M nums/sec       (70× speedup)
Phase 2 (10 GPP nodes): 1.5B nums/sec      (700× speedup)
Phase 3 (1 GPU):       300M nums/sec       (variable)
Phase 4 (100 GPUs):    30B nums/sec        (maximum)
```

## 📋 Deployment Steps

### Step 1: Local Test (5 minutes)
```bash
cd /home/aruns/Desktop/MN25/MNv3
bash test_local.sh
```

**Expected output:**
- Build successful
- Test runs in ~1 second
- Throughput >100K nums/sec
- Speedup 3-4× (on 4 threads)

### Step 2: Deploy to MareNostrum (5 minutes)
```bash
# Copy MNv3 folder to MareNostrum
cd /home/aruns/Desktop/MN25
scp -r MNv3/ nct01225@glogin1.bsc.es:/gpfs/projects/nct_352/nct01225/collatz/

# SSH to MareNostrum
ssh nct01225@glogin1.bsc.es
cd /gpfs/projects/nct_352/nct01225/collatz/MNv3
```

### Step 3: Test on Login Node (2 minutes)
```bash
# Load Intel compiler (if available)
module load intel/2023.2.0  # or whatever version is available

# Build
bash build_openmp.sh

# Quick test (4 threads, 1M numbers)
./V1.5-openmp 0 1000000 --threads 4 --tag login_test
```

**Expected:**
- Build completes without errors
- Test runs in <10 seconds
- Throughput >1M nums/sec
- Speedup ~3-4×

### Step 4: Submit 1-Node Job (10 minutes)
```bash
# Submit to queue
sbatch slurm_openmp_1node.slurm

# Check status
squeue -u nct01225

# Watch output (once running)
watch tail -30 collatz_omp_*.out
```

**Expected:**
- Job starts within 5 minutes (debug queue)
- Completes in ~1 minute (100M numbers)
- Throughput 100-200M nums/sec
- Speedup 50-90×

### Step 5: Validate Results
```bash
# Check output
cat collatz_omp_*.out

# Verify:
# ✓ All 100M numbers tested
# ✓ Throughput >100M nums/sec
# ✓ No crashes or errors
# ✓ Edge cases logged (overflow/fuse files)
```

## 🚀 Next Phase Decision

### Option A: Scale GPP with MPI (Recommended for Hackathon)
**When to choose:**
- Phase 1 achieved >100M nums/sec ✓
- You want guaranteed scaling
- Time is limited (hackathon deadline)

**Action:**
- Create V1.6-mpi.cpp (MPI multi-node)
- Test with 2 nodes, then 10 nodes
- Expected: 1-2B nums/sec

### Option B: Jump to GPU/ACC (Maximum Performance)
**When to choose:**
- You need >1B nums/sec immediately
- You have GPU experience
- You have 2-3 hours for development

**Action:**
- Create V1.7-cuda.cu (CUDA single GPU)
- Port hot path to GPU kernel
- Expected: 300M+ nums/sec per GPU

### Option C: Conservative Approach
**When to choose:**
- Phase 1 results are uncertain
- You want solid foundation
- Hackathon has multiple days

**Action:**
- Finish MPI multi-node first (Phase 2)
- Validate with 10 nodes (1.5B nums/sec)
- Then attempt GPU (Phase 3)

## 📊 Performance Analysis

### Bottlenecks to Watch

**Memory Bandwidth:**
```
- L3 cache: ~200 GB/s (shared by 56 cores)
- DDR5: ~200 GB/s per socket
- Memo table: 4MB (fits in L3)
- Expected: 60-80% scaling efficiency
```

**Load Imbalance:**
```
- Some numbers take longer (high step counts)
- Dynamic scheduling helps
- Monitor: max_steps per thread
```

**Cold Queue Overhead:**
```
- Rare edge cases (overflow, fuse hits)
- Thread-local queues minimize impact
- Should be <1% of total time
```

### Optimization Flags

**For Intel Compiler (icpx):**
```bash
-O3                # Maximum optimization
-xHost             # Target specific CPU (Sapphire Rapids)
-qopenmp           # Enable OpenMP
-ipo               # Interprocedural optimization
-no-prec-div       # Fast division
-fp-model fast=2   # Fast floating point
```

**For GCC:**
```bash
-O3                # Maximum optimization
-march=native      # Target CPU
-mtune=native      # Tune for CPU
-fopenmp           # Enable OpenMP
-ffast-math        # Fast math
-funroll-loops     # Unroll loops
```

## 🐛 Troubleshooting

### Job Won't Start
```bash
# Check queue status
squeue -u nct01225

# Check account
sacctmgr show assoc where user=nct01225

# Try different QOS
# Edit slurm file: --qos=gp_bsccs instead of gp_debug
```

### Low Performance
```bash
# Check thread binding
export OMP_DISPLAY_ENV=TRUE
./V1.5-openmp 0 1000 --threads 4

# Try different scheduling
export OMP_SCHEDULE="dynamic,1024"  # Larger chunks

# Check NUMA
numactl --cpunodebind=0 --membind=0 ./V1.5-openmp 0 1000000 --threads 56
```

### Build Errors
```bash
# Check compiler version
g++ --version  # Need 9.0+
icpx --version  # Need 2021.1+

# Try minimal build
g++ -O2 -fopenmp V1.5-openmp.cpp -o test

# Check for missing libraries
ldd V1.5-openmp
```

## 📈 Scaling Projections

### Single Node (112 cores)
```
Ideal:     112× speedup = 246M nums/sec
Expected:   70× speedup = 154M nums/sec
Minimum:    50× speedup = 110M nums/sec
```

### 10 Nodes (1,120 cores)
```
Ideal:     1,120× speedup = 2.46B nums/sec
Expected:    700× speedup = 1.54B nums/sec
Minimum:     500× speedup = 1.10B nums/sec
```

### 100 GPUs (Hopper)
```
Ideal:     100× 500M = 50B nums/sec
Expected:  100× 300M = 30B nums/sec
Minimum:   100× 200M = 20B nums/sec
```

## ✅ Success Criteria

### Phase 1 Success (OpenMP)
- [x] Builds without errors
- [x] Runs on login node (4 threads)
- [ ] Runs on compute node (112 threads)
- [ ] Achieves >100M nums/sec
- [ ] No crashes or data corruption
- [ ] Edge cases logged correctly

### Ready for Phase 2 (MPI)
- [ ] Phase 1 all criteria met
- [ ] Throughput >150M nums/sec per node
- [ ] Results validated (compare with V1.4b)
- [ ] No mysterious errors

### Ready for Phase 3 (GPU)
- [ ] Comfortable with CUDA programming
- [ ] Have 2-3 hours for development
- [ ] GPU resources available (ACC partition)
- [ ] Or: MPI not meeting performance target

## 🎓 Learning Points

### Key Insights
1. **Embarrassingly parallel** = Perfect for HPC
2. **Read-only data** = Zero synchronization overhead
3. **Dynamic scheduling** = Automatic load balancing
4. **Thread-local** = Minimal critical sections

### Common Mistakes to Avoid
1. ❌ Using static scheduling (load imbalance)
2. ❌ Shared cold queues (mutex contention)
3. ❌ Too small chunks (scheduling overhead)
4. ❌ Forgetting OMP_PROC_BIND (NUMA issues)

### Best Practices
1. ✅ Test locally first
2. ✅ Start with debug queue (faster)
3. ✅ Monitor resource usage
4. ✅ Validate results against sequential
5. ✅ Log edge cases for analysis

## 📞 Need Help?

### Resources
- **MareNostrum Docs**: https://www.bsc.es/user-support/mn5.php
- **Intel Compiler**: https://www.intel.com/content/www/us/en/developer/tools/oneapi/dpc-compiler.html
- **OpenMP Guide**: https://www.openmp.org/specifications/

### Contact
- **Hackathon Mentors**: Ask during sessions
- **BSC Support**: support@bsc.es
- **Your Team**: Discuss strategy together

---

**GOOD LUCK! 🚀**

Remember: Start simple (OpenMP), validate, then scale!
