# MNv2 - Clean MareNostrum 5 Implementation

## Quick Start (Hackathon Team)

### 1. Configure Your Resources
Edit `config.hpp`:
```cpp
constexpr int GPP_NODES = 2;  // Your GPP allocation
constexpr int ACC_NODES = 0;  // Your ACC allocation (Phase 2)
```

### 2. Build Phase 1 (GPP nodes, CPU-only)
```bash
./build_gpp.sh
```

### 3. Test Locally
```bash
# Test with 4 MPI ranks (1 master + 3 workers)
mpirun -np 4 ./collatz_mpi_gpp 0 1000000 test_run
```

### 4. Deploy to MareNostrum
```bash
# Package
tar czf mnv2.tar.gz MNv2/

# Transfer
scp mnv2.tar.gz nct01225@glogin1.bsc.es:~/

# SSH and extract
ssh nct01225@glogin1.bsc.es
cd /gpfs/projects/nct_352/nct01225/collatz/  # Projects dir (backed up!)
tar xzf ~/mnv2.tar.gz
cd MNv2/

# Create output directory for logs (in scratch - temporary)
mkdir -p /gpfs/scratch/nct_352/nct01225/collatz_output

# Submit job
sbatch slurm_gpp.slurm
```

---

## MareNostrum 5 Filesystem Guide

**Important:** Use the correct filesystem for each purpose!

| Filesystem | Purpose | Backed Up | What to Store |
|------------|---------|-----------|---------------|
| `/gpfs/projects/nct_352/$USER/` | Project data | ✅ Yes (3-4 days) | **Source code, executables** ⭐ |
| `/gpfs/scratch/nct_352/$USER/` | Temporary | ❌ No | Job output logs, temp data |
| `/gpfs/home/$USER/` | Home | ✅ Yes (~1 day) | Config files, small scripts |
| `$TMPDIR` (local SSD) | Job-local | ❌ No | Fast I/O during job (auto-cleaned) |

**Our deployment:**
- **Build/Run from:** `/gpfs/projects/nct_352/nct01225/collatz/MNv2/`
- **Logs go to:** `/gpfs/scratch/nct_352/nct01225/collatz_output/`
- **Temp files:** `$TMPDIR` (automatically set by SLURM)

---

## Architecture

```
MPI Master (Rank 0)
   ├── Partition work across N workers
   ├── Receive results
   └── Aggregate statistics
   
MPI Workers (Rank 1..N)
   ├── Receive work assignment
   ├── Process with OpenMP threads (GPP: 112 cores, ACC: 80 cores)
   │   └── V1.4b-openmp core: hot path + cold queues
   ├── Send results to master
   └── Each worker is independent (no inter-worker communication)
```

---

## Files

- **config.hpp** - Central configuration (node counts, cores, memo size, fuses)
- **collatz_mpi_gpp.cpp** - Phase 1: MPI+OpenMP for GPP nodes (COMPLETE)
- **build_gpp.sh** - Build script for GPP phase
- **slurm_gpp.slurm** - SLURM job script for MareNostrum 5
- **test_local.sh** - Local MPI testing helper

---

## Phase Roadmap

### Phase 1: GPP Nodes ✅ (Current)
- **Goal:** Get working ASAP with proven V1.4b-openmp
- **Nodes:** 2-10 GPP nodes (112 cores each)
- **Expected:** ~250M nums/sec per node = 500M-2.5B total
- **Status:** READY TO DEPLOY

### Phase 2: ACC Nodes (Next)
- **Goal:** Add GPU acceleration
- **Nodes:** 1-3 ACC nodes (4× H100 + 80 cores each)
- **Expected:** ~500M-1B nums/sec per node (GPU-bound)
- **File:** `collatz_mpi_acc.cu` (CUDA kernel + CPU cold queues)

### Phase 3: Hybrid (Final)
- **Goal:** Mix GPP + ACC nodes intelligently
- **Auto-detection:** Workers detect their node type and optimize accordingly
- **Expected:** 3-5B nums/sec (combined throughput)

---

## Performance Targets

| Configuration | Nodes | Cores | GPUs | Expected Throughput |
|---------------|-------|-------|------|---------------------|
| 2 GPP         | 2     | 224   | 0    | 500M nums/sec       |
| 5 GPP         | 5     | 560   | 0    | 1.2B nums/sec       |
| 10 GPP        | 10    | 1,120 | 0    | 2.5B nums/sec       |
| 1 ACC         | 1     | 80    | 4    | 500M-1B nums/sec    |
| 5 GPP + 2 ACC | 7     | 720   | 8    | 3-4B nums/sec       |

**⚠️ Workload Size Critical for Accurate Benchmarks!**

MPI has fixed overhead (~50ms). For realistic performance:
- ❌ **<10M seeds:** Overhead dominates → appears slow
- ✅ **100M-1B seeds:** Overhead <5% → shows true speed
- ✅ **10B+ seeds:** Production-scale on MareNostrum

**Example (local 3 MPI ranks, 4 threads):**
- 100K seeds → 0.6M nums/sec (misleading!)
- 10M seeds → 15M nums/sec (better)
- 100M seeds → 18M nums/sec (realistic)

---

## Troubleshooting

**Build fails:**
- Check module load: `module load intel/2023.2 impi/2021.10`
- Verify compiler: `which mpicxx`

**Job fails immediately:**
- Check account/QOS: `squeue -u $USER`
- Verify node availability: `sinfo -p gpp`

**Low performance:**
- Check thread count: `echo $OMP_NUM_THREADS` (should be 112 for GPP, 80 for ACC)
- Verify affinity: `echo $OMP_PROC_BIND` (should be "close")
- Check memo table size in config.hpp (default 2^20 = 4MB is good)

**MPI errors:**
- Test locally first: `mpirun -np 4 ./collatz_mpi_gpp 0 10000 test`
- Check network: `ibstat` (should show InfiniBand UP)

---

## Key Differences from Old `marenostrum/` Implementation

| Old (marenostrum/) | New (MNv2/) |
|--------------------|-------------|
| ❌ Hollow placeholders | ✅ Complete working code |
| ❌ sleep(2) simulation | ✅ Real V1.4b-openmp computation |
| ❌ Scattered config | ✅ Central config.hpp |
| ❌ Complex file structure | ✅ Self-contained single file |
| ❌ No local testing | ✅ Easy local MPI testing |
| ❌ Hardcoded values | ✅ Flexible configuration |

---

## License

Same as parent project (see ../LICENSE)
