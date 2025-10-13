# MNv3 - MareNostrum 5 Scalable Deployment

## ✅ Phase 1: OpenMP Single-Node (GPP) - **COMPLETE!**

**Result:** 137M nums/sec (62× speedup) on 1 node × 112 cores

## 🚀 Phase 2: MPI Multi-Node (GPP) - **READY TO DEPLOY!**

**Expected:** 1.37B nums/sec (623× speedup) on 10 nodes × 1,120 cores

---

## Phase 1: OpenMP Single-Node (GPP) ✅ **COMPLETE**

### Quick Start
```bash
# Build
bash build_openmp.sh

# Local test (4 threads)
./V1.5-openmp 0 1000000 --threads 4 --small-limit 20 --tag local

# Submit to MareNostrum GPP (1 node, 112 cores)
sbatch slurm_openmp_1node.slurm
```

### Expected Performance
- **Sequential (V1.4b)**: 2.2M nums/sec
- **OpenMP 56 threads**: 100-150M nums/sec (45-70× speedup)
- **OpenMP 112 threads**: 150-200M nums/sec (70-90× speedup)

### Architecture
- **Target**: GPP nodes (Intel Sapphire Rapids 8480+)
- **Cores**: 2× 56 cores = 112 cores per node
- **Memory**: 256GB DDR5
- **Memo table**: 2^20 = 4MB (fits in L3 cache)

### Key Features
1. **Thread-safe**: Read-only memo table after precompute
2. **Dynamic scheduling**: Automatic load balancing
3. **Thread-local stats**: Minimal synchronization overhead
4. **Batch cold queue**: Per-thread edge case handling
5. **Crash-safe logs**: All edge cases logged

## Phase 2: MPI Multi-Node (GPP) ⭐ **READY NOW**

### Status: READY TO DEPLOY
- **Target**: 10 nodes = 1,120 cores
- **Expected**: 1.37 BILLION nums/sec (based on Phase 1: 137M/node)
- **Strategy**: Embarrassingly parallel (zero inter-node communication)
- **Files**: V1.6-mpi-openmp.cpp, build_mpi.sh, slurm scripts ready

### Quick Start
```bash
cd /gpfs/projects/nct_352/nct01225/collatz/MNv3
module load impi/2021.9.0
bash build_mpi.sh
sbatch slurm_mpi_2nodes.slurm      # Validate with 2 nodes first
sbatch slurm_mpi_10nodes.slurm     # Then run 10-node production
```

## Phase 3: CUDA Single-GPU (ACC) - Coming Soon

### Plan
- **Target**: 1 NVIDIA Hopper GPU
- **Expected**: 50-500M nums/sec per GPU
- **Strategy**: Warp-level parallelism, coalesced memory

## Phase 4: Multi-GPU + MPI Hybrid (ACC) - Final Target

### Plan
- **Target**: 50-100 nodes = 200-400 Hopper GPUs
- **Expected**: 10-50 BILLION nums/sec
- **Strategy**: Maximum throughput for hackathon

## File Structure
```
MNv3/
├── README.md                    # This file
├── V1.5-openmp.cpp              # Phase 1: OpenMP single-node
├── build_openmp.sh              # Build script
├── slurm_openmp_1node.slurm     # SLURM job (1 node, 112 cores)
├── V1.6-mpi.cpp                 # Phase 2: MPI multi-node (TODO)
├── V1.7-cuda.cu                 # Phase 3: CUDA single-GPU (TODO)
└── V1.8-cuda-mpi.cu             # Phase 4: Multi-GPU hybrid (TODO)
```

## Workflow

### Step 1: Local Testing
```bash
# On your laptop/local machine
bash build_openmp.sh
./V1.5-openmp 0 100000 --threads 4 --tag local_test
```

### Step 2: MareNostrum Login Node Test
```bash
# SSH to MareNostrum
ssh nct01225@glogin1.bsc.es
cd /gpfs/projects/nct_352/nct01225/collatz/MNv3

# Quick test (4 threads, 1M numbers)
./V1.5-openmp 0 1000000 --threads 4 --tag login_test
```

### Step 3: Submit 1-Node Job
```bash
# Submit to queue
sbatch slurm_openmp_1node.slurm

# Check status
squeue -u $USER

# Watch output
tail -f collatz_omp_*.out
```

### Step 4: Scale Up
- Once 1-node works, proceed to Phase 2 (MPI)
- Or jump to Phase 3 (GPU) if you want maximum throughput

## Performance Monitoring

### Expected Output
```
=== V1.5-openmp Results (OpenMP Parallel) ===
Threads:      112
Tested:       100000000 numbers
Time:         600000 ms
Throughput:   166666666 nums/sec
Speedup:      75.8× (vs 2.2M/s sequential)
Avg steps:    350.5
Max steps:    2254 at n=...
Peak value:   ...
```

### Scaling Efficiency
- **Linear scaling**: 112 threads = 112× speedup (ideal)
- **Expected scaling**: 60-80× speedup (realistic)
- **Bottlenecks**: Memory bandwidth, cache misses

## Troubleshooting

### Build Errors
```bash
# Check compiler
g++ --version  # Should be GCC 9+ or Intel 2021+

# Try without OpenMP first
g++ -O3 V1.5-openmp.cpp -o test_nomp  # Will fail, but shows compile errors

# Check module availability
module avail intel
module load intel/2023.2.0
```

### Runtime Errors
```bash
# Check OpenMP availability
echo $OMP_NUM_THREADS  # Should be set

# Test with fewer threads
./V1.5-openmp 0 1000 --threads 2 --tag debug

# Check memory
free -h  # Should have >4MB available for memo table
```

### Performance Issues
```bash
# Enable verbose output
export OMP_DISPLAY_ENV=TRUE
./V1.5-openmp 0 1000000 --threads 4

# Profile with perf (if available)
perf stat -d ./V1.5-openmp 0 1000000 --threads 56

# Check NUMA balance
numactl --hardware
```

## Next Steps

After validating Phase 1:

1. **If throughput meets target** (>100M nums/sec):
   - ✅ Proceed to Phase 2 (MPI) for multi-node scaling
   
2. **If throughput exceeds expectations** (>200M nums/sec):
   - ✅ Consider staying on GPP with Phase 2 (cheaper, easier)
   
3. **If you need MAXIMUM throughput** (>1B nums/sec):
   - ✅ Jump to Phase 3 (CUDA) for GPU acceleration

## Contact
- **Hackathon**: MareNostrum 5 Hackathon (October 14, 2025)
- **Account**: nct_352
- **User**: nct01225
- **Support**: Check MareNostrum documentation or ask mentors

## License
Based on V1.4b Collatz Conjecture implementation
