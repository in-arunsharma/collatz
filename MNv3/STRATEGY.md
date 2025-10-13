# 🎯 MareNostrum 5 Collatz Scaling Strategy - EXECUTIVE SUMMARY

**Date:** October 14, 2025  
**Hackathon:** MareNostrum 5  
**Baseline:** V1.4b @ 2.2M nums/sec (sequential)  
**Goal:** Maximum throughput for Collatz conjecture verification

---

## ✅ RECOMMENDED STRATEGY: Start with GPP, then decide

### Phase 1: OpenMP Single-Node GPP ⭐ **START HERE** (30 min)
- **Status:** ✅ **READY TO DEPLOY**
- **Target:** 1 node × 112 cores
- **Expected:** 150-200M nums/sec (70-90× speedup)
- **Risk:** LOW ✅
- **Files:** `MNv3/` folder ready
- **Action:** Deploy now, validate in 30 minutes

### Phase 2: MPI Multi-Node GPP (2 hours)
- **Status:** ⏳ Create after Phase 1 success
- **Target:** 10 nodes × 112 cores = 1,120 cores
- **Expected:** 1-2 BILLION nums/sec (700× speedup)
- **Risk:** LOW ✅
- **When:** If Phase 1 achieves >100M nums/sec

### Phase 3: CUDA Single-GPU ACC (3-4 hours)
- **Status:** ⏳ Alternative to Phase 2
- **Target:** 1 NVIDIA Hopper GPU
- **Expected:** 200-500M nums/sec per GPU
- **Risk:** MEDIUM ⚠️ (new CUDA code)
- **When:** If you want maximum single-node performance

### Phase 4: Multi-GPU + MPI Hybrid (1 day)
- **Status:** 🎯 Final target for maximum throughput
- **Target:** 50-100 nodes = 200-400 GPUs
- **Expected:** 10-50 BILLION nums/sec
- **Risk:** MEDIUM-HIGH ⚠️
- **When:** After Phase 3 validation

---

## 📊 PERFORMANCE COMPARISON

| Phase | Hardware | Throughput | Speedup | Time to Deploy |
|-------|----------|------------|---------|----------------|
| V1.4b | 1 CPU core | 2.2M/s | 1× | - |
| **Phase 1** | **1 GPP node** | **150M/s** | **70×** | **30 min** ⭐ |
| Phase 2 | 10 GPP nodes | 1.5B/s | 700× | 2 hours |
| Phase 3 | 1 GPU | 300M/s | 140× | 3 hours |
| Phase 4 | 100 GPUs | 30B/s | 14,000× | 1 day |

---

## 🚀 IMMEDIATE ACTION PLAN

### Right Now (5 minutes)
```bash
cd /home/aruns/Desktop/MN25
scp -r MNv3/ nct01225@glogin1.bsc.es:/gpfs/projects/nct_352/nct01225/collatz/
```

### On MareNostrum (10 minutes)
```bash
ssh nct01225@glogin1.bsc.es
cd /gpfs/projects/nct_352/nct01225/collatz/MNv3

# Test build
bash build_openmp.sh

# Quick validation
./V1.5-openmp 0 1000000 --threads 4 --tag login_test

# Submit 1-node job
sbatch slurm_openmp_1node.slurm
```

### Monitor (15 minutes)
```bash
# Check queue
squeue -u nct01225

# Watch output
watch tail -30 collatz_omp_*.out
```

---

## 🎯 DECISION TREE

```
Phase 1 Result
│
├─ >150M nums/sec ✅
│  ├─ Option A: Phase 2 (MPI) → 1.5B nums/sec [SAFE]
│  └─ Option B: Phase 3 (GPU) → 300M nums/sec [AMBITIOUS]
│
├─ 100-150M nums/sec ⚠️
│  ├─ Option A: Optimize Phase 1, then MPI
│  └─ Option B: Jump to GPU (better single-node performance)
│
└─ <100M nums/sec ❌
   └─ Debug: Check NUMA, thread binding, profiling
```

---

## 💡 WHY START WITH GPP?

### ✅ Advantages
1. **Fast deployment** (30 min vs 3 hours for GPU)
2. **Low risk** (OpenMP is simple, stable)
3. **Validates parallelization** (test before GPU complexity)
4. **Guaranteed resources** (less queue competition)
5. **Foundation for MPI** (same code, add MPI wrapper)
6. **Good performance** (150M nums/sec is respectable)

### ⚠️ Why Not GPU First?
1. **Higher development time** (3-4 hours)
2. **More complex** (CUDA kernels, memory management)
3. **Higher risk** (bugs harder to debug)
4. **Queue competition** (ACC partition may be busy)
5. **Overkill for hackathon?** (GPP might be sufficient)

### 🎯 When to Choose GPU?
- You need >1B nums/sec immediately
- You have CUDA experience
- You have 3+ hours available
- GPP results are disappointing (<100M nums/sec)

---

## 📁 FILES READY FOR DEPLOYMENT

```
MNv3/
├── V1.5-openmp.cpp              OpenMP parallelized version
├── build_openmp.sh              Build script (g++ or icpx)
├── slurm_openmp_1node.slurm     SLURM job (1 node, 112 cores)
├── test_local.sh                Local testing (4 threads)
├── README.md                    Quick reference
├── DEPLOYMENT_GUIDE.md          Detailed step-by-step
└── PHASE1_COMPLETE.md           Status report
```

**Status:** ✅ All files created and tested locally

---

## 🔬 TECHNICAL HIGHLIGHTS

### V1.5-openmp Features
- **Thread-safe:** Read-only memo table after precompute
- **Load-balanced:** Dynamic scheduling (512-number chunks)
- **Minimal overhead:** Thread-local statistics and cold queues
- **Edge cases:** Preserved overflow and cycle detection
- **Scalable:** No synchronization in hot path

### Optimization Flags
```bash
# GCC
-O3 -march=native -mtune=native -fopenmp -ffast-math

# Intel (recommended on MareNostrum)
-O3 -xHost -qopenmp -ipo -no-prec-div -fp-model fast=2
```

### SLURM Configuration
```bash
--nodes=1                   # Single node
--cpus-per-task=112         # All 112 cores
--time=02:00:00             # 2 hours (plenty)
--qos=gp_debug              # Fast queue
```

---

## ⏱️ TIMELINE TO SUCCESS

```
T+0:00    Deploy MNv3 to MareNostrum
T+0:10    Build and test on login node
T+0:15    Submit 1-node job
T+0:20    Job starts (debug queue)
T+0:21    Job completes (100M numbers in ~1 min)
T+0:30    Analyze results, decide next phase
---
Phase 1:  30 minutes total

T+0:30    Start Phase 2 (MPI) development
T+1:30    Submit 2-node test
T+2:00    Submit 10-node production run
T+2:30    Achieve 1.5B nums/sec
---
Phase 2:  2.5 hours total from start

OR

T+0:30    Start Phase 3 (CUDA) development
T+2:00    Test single GPU kernel
T+3:00    Submit 1-GPU job
T+3:30    Achieve 300M nums/sec per GPU
T+4:00    Start multi-GPU (Phase 4)
---
Phase 3:  4 hours total from start
```

---

## ✅ SUCCESS CRITERIA

### Phase 1 (Minimum)
- [ ] Job completes without errors
- [ ] Throughput >100M nums/sec
- [ ] All edge cases logged correctly
- [ ] Results validate against V1.4b

### Phase 1 (Target)
- [ ] Throughput >150M nums/sec
- [ ] Speedup >70×
- [ ] Scaling efficiency >60%
- [ ] Ready for Phase 2

### Overall Hackathon (Minimum)
- [ ] Verify 10 billion numbers
- [ ] Throughput >500M nums/sec sustained
- [ ] Complete scientific report

### Overall Hackathon (Target)
- [ ] Verify 100 billion numbers
- [ ] Throughput >1B nums/sec sustained
- [ ] Demonstrate GPU scaling

---

## 🆘 TROUBLESHOOTING QUICK REFERENCE

### Job won't start
```bash
squeue -u nct01225           # Check queue
sacctmgr show assoc user=$USER  # Check account
# Try: --qos=gp_bsccs instead of gp_debug
```

### Low performance
```bash
export OMP_PROC_BIND=close   # Thread pinning
export OMP_PLACES=cores      # NUMA awareness
numactl --cpunodebind=0 ./V1.5-openmp ...  # Force NUMA node
```

### Build errors
```bash
module load intel/2023.2.0   # Use Intel compiler
g++ --version                # Check version (need 9.0+)
```

---

## 📞 SUPPORT

- **MareNostrum Docs:** https://www.bsc.es/user-support/mn5.php
- **Hackathon Mentors:** Ask during sessions
- **BSC Support:** support@bsc.es

---

## 🎓 KEY INSIGHTS

1. **Start simple, scale incrementally** (OpenMP → MPI → GPU)
2. **Validate each phase** (don't skip ahead without proof)
3. **Embarrassingly parallel = HPC paradise** (zero communication overhead)
4. **Read-only data = zero synchronization** (perfect scaling)
5. **Dynamic scheduling = automatic load balancing** (handles imbalanced work)

---

## 🚀 FINAL RECOMMENDATION

**FOR HACKATHON SUCCESS:**

1. ✅ **Deploy Phase 1 NOW** (30 minutes to validation)
2. ⏰ **Wait for results** (likely >150M nums/sec)
3. 🚀 **Proceed to Phase 2** (MPI multi-node → 1.5B nums/sec)
4. 🎯 **Done!** (1.5B nums/sec is excellent for hackathon)

**FOR MAXIMUM PERFORMANCE:**

1. ✅ **Deploy Phase 1** (validate parallelization works)
2. 🎮 **Jump to Phase 3** (GPU development in parallel)
3. 💪 **Scale to Phase 4** (100 GPUs → 30B nums/sec)
4. 🏆 **Win hackathon!** (Highest throughput possible)

---

**STATUS: 🟢 READY TO DEPLOY**

**Next command:**
```bash
cd /home/aruns/Desktop/MN25
scp -r MNv3/ nct01225@glogin1.bsc.es:/gpfs/projects/nct_352/nct01225/collatz/
```

**Good luck! 🚀**
