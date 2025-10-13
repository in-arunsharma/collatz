# Phase 1 - Quick Start Summary

## 🎯 Goal
Deploy proven V1.4b-openmp code on **5 MareNostrum GPP nodes** using MPI for multi-node scaling.

**Expected Performance:** 1.25B numbers/sec (5 nodes × 250M/sec)

---

## ✅ What's Ready

### Files Created (all in `marenostrum/` folder):

| File | Purpose | Status |
|------|---------|--------|
| `marenostrum_config.hpp` | Configuration constants | ✅ Complete |
| `node_config.hpp` | Hardware auto-detection | ✅ Complete |
| `mpi_coordinator.hpp` | MPI work distribution | ✅ Complete |
| `worker_gpp_node.cpp` | GPP node worker class | ⚠️ Needs V1.4b integration |
| `collatz_mpi_gpp.cpp` | Main MPI program | ✅ Complete |
| `build_phase1_gpp.sh` | Build script | ✅ Ready |
| `slurm_phase1_gpp.slurm` | SLURM job submission | ✅ Ready |
| `test_local_mpi.sh` | Local MPI testing | ✅ Ready |
| `INTEGRATION_GUIDE.md` | How to link V1.4b core | ✅ Complete |
| `EXECUTION_GUIDE.md` | Full deployment guide | ✅ Complete |
| `README.md` | Architecture overview | ✅ Complete |

---

## 🔧 What You Need to Do (Monday AM)

### 1. Integrate V1.4b Core (~10 minutes)

**Quick option:** Edit `worker_gpp_node.cpp`

Add at top of file:
```cpp
// Include the proven V1.4b-openmp algorithm
#include "../parallel/V1.4b-openmp.cpp"
```

Replace placeholder in `process_work()` function:
```cpp
// OLD (line ~45):
uint32_t steps = (n % 100);  // Dummy computation

// NEW:
uint32_t steps = process_number(n);  // Use V1.4b function
```

**See `INTEGRATION_GUIDE.md` for details**

### 2. Build & Test Locally (~5 minutes)

```bash
cd marenostrum/

# Build
./build_phase1_gpp.sh

# Test with 2 MPI ranks
./test_local_mpi.sh 2
```

**Expected test output:**
```
✓ Test PASSED
Ready for MareNostrum deployment!
```

### 3. Deploy to MareNostrum (~15 minutes)

```bash
# From your local machine
cd ~/Desktop/MN25
tar czf collatz_mn5.tar.gz marenostrum/ parallel/V1.4b-openmp.cpp
scp collatz_mn5.tar.gz <username>@mn5.bsc.es:~/

# SSH to MareNostrum
ssh <username>@mn5.bsc.es
tar xzf collatz_mn5.tar.gz
cd marenostrum/
./build_phase1_gpp.sh
sbatch slurm_phase1_gpp.slurm

# Monitor job
squeue -u $USER
tail -f collatz_phase1_*.out
```

---

## 📊 How It Works

### Architecture
```
[Master Node]      [Worker Node 1]   [Worker Node 2]   [Worker Node 3]   [Worker Node 4]
MPI Rank 0         MPI Rank 1        MPI Rank 2        MPI Rank 3        MPI Rank 4
112 threads        112 threads       112 threads       112 threads       112 threads
    ↓                  ↓                 ↓                 ↓                 ↓
Seeds 0-2B         Seeds 2B-4B       Seeds 4B-6B       Seeds 6B-8B       Seeds 8B-10B
```

**Total:** 5 nodes × 112 threads = **560 parallel threads**

### Work Distribution (Embarrassingly Parallel)
- Master (rank 0) divides seed range equally among workers
- Each worker processes its chunk independently (no communication during computation)
- Workers send results back to master at the end
- Master aggregates and writes final JSON

---

## 📋 Pre-Flight Checklist

Before submitting to MareNostrum:

- [ ] V1.4b core integrated into `worker_gpp_node.cpp`
- [ ] Local build successful (`./build_phase1_gpp.sh`)
- [ ] Local MPI test passed (`./test_local_mpi.sh 2`)
- [ ] Code transferred to MareNostrum
- [ ] Built on MareNostrum successfully
- [ ] SLURM script parameters adjusted (seed range, output prefix)
- [ ] Job submitted (`sbatch slurm_phase1_gpp.slurm`)

---

## 🚀 Expected Results

### Test Run (10B numbers, 5 GPP nodes)
- **Throughput:** ~1.25B numbers/sec
- **Walltime:** ~8 seconds computation + ~2 seconds overhead
- **Output:** JSON files with cycle seeds, overflows, max steps

### Scaling (if you get more nodes)
- 10 GPP nodes → **2.5B numbers/sec**
- 15 GPP nodes → **3.75B numbers/sec**

---

## 📖 Documentation

- **INTEGRATION_GUIDE.md** - How to connect V1.4b-openmp core
- **EXECUTION_GUIDE.md** - Full deployment walkthrough
- **README.md** - Architecture and design decisions

---

## 🆘 Quick Troubleshooting

**Build fails?**
```bash
# Check MPI compiler
which mpic++
mpic++ --version
```

**Local test fails?**
```bash
# Try fewer ranks
./test_local_mpi.sh 2

# Check MPI installation
mpirun --version
```

**SLURM job pending?**
```bash
# Check queue position
squeue -u $USER

# Try fewer nodes first
#SBATCH --nodes=2  # Edit slurm script
```

**Low performance?**
```bash
# Check OpenMP threads
export OMP_DISPLAY_ENV=TRUE
```

---

## ⏱️ Monday Timeline

| Time | Task | Duration |
|------|------|----------|
| 9:00 | Integrate V1.4b core | 10 min |
| 9:10 | Build & test locally | 5 min |
| 9:15 | Transfer to MareNostrum | 5 min |
| 9:20 | Build on MareNostrum | 5 min |
| 9:25 | Submit test job (1B numbers) | 2 min run |
| 9:30 | Analyze results | 10 min |
| 9:40 | Submit production run | Ready! |

**Phase 1 completion: ~10:00 AM**

Then move to Phase 2 (ACC nodes with CUDA) 🚀

---

**You're 30 minutes away from 1.25 billion numbers/sec!** 💪
