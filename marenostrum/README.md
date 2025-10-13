# MareNostrum 5 Deployment
## Multi-Node MPI + OpenMP + CUDA Hybrid Implementation

### Architecture Overview

```
MPI Layer (Inter-Node)
    ↓
┌─────────────────┬─────────────────┐
│   GPP Nodes     │   ACC Nodes     │
│   (CPU-only)    │   (GPU+CPU)     │
└─────────────────┴─────────────────┘
    ↓                   ↓
OpenMP (112 threads)   CUDA (4 GPUs) + OpenMP (80 threads)
```

### Hardware Specifications

#### GPP Nodes (General Purpose Partition)
- **CPUs:** 2× Intel Sapphire Rapids 8480+ @ 2.0 GHz
- **Cores:** 56 cores/CPU × 2 = **112 cores per node**
- **Memory:** 256GB DDR5 (some nodes have 1TB)
- **Network:** NDR200 (100 Gb/s per node)
- **Use case:** Pure CPU parallelization with V1.4b-openmp

#### ACC Nodes (Accelerated Partition)
- **CPUs:** 2× Intel Sapphire Rapids 8460Y+ @ 2.3 GHz
- **Cores:** 40 cores/CPU × 2 = **80 cores per node**
- **GPUs:** 4× NVIDIA Hopper H100
- **GPU Memory:** 64GB HBM2 per GPU (256GB total)
- **System Memory:** 512GB DDR5
- **Network:** 4× NDR200 (800 Gb/s total)
- **Use case:** GPU hot path + CPU cold queue processing

### Implementation Phases

#### Phase 1: GPP-Only (Monday Morning) - PRIORITY
**Goal:** Get working ASAP with tested V1.4b-openmp code
- **Nodes:** 5-10 GPP nodes
- **Strategy:** MPI work distribution + OpenMP per-node
- **Expected:** 1-2 billion nums/sec
- **Files:** `V1.6-mpi-cpu.cpp` (wrapper around V1.4b-openmp)

#### Phase 2: ACC-Only (Monday Afternoon)
**Goal:** Add GPU acceleration
- **Nodes:** 2-5 ACC nodes
- **Strategy:** CUDA hot path + OpenMP cold queues
- **Expected:** 400M-1B nums/sec (GPU-bound)
- **Files:** `V1.6-mpi-gpu.cpp` (CUDA + OpenMP hybrid)

#### Phase 3: Hybrid (Monday Evening)
**Goal:** Combine everything
- **Nodes:** 5 ACC + 10 GPP = 15 nodes total
- **Strategy:** Auto-detect node type, optimize workload
- **Expected:** 3-4 billion nums/sec
- **Files:** `V1.6-mpi-hybrid.cpp` (unified)

### File Structure

```
marenostrum/
├── README.md                      # This file
├── node_config.hpp                # Hardware detection & config
├── mpi_coordinator.hpp            # MPI work distribution
├── worker_gpp.cpp                 # GPP node worker (OpenMP)
├── worker_acc.cpp                 # ACC node worker (CUDA+OpenMP)
├── V1.6-mpi-cpu.cpp              # Phase 1: GPP-only executable
├── V1.6-mpi-gpu.cpp              # Phase 2: ACC-only executable
├── V1.6-mpi-hybrid.cpp           # Phase 3: Unified executable
├── build_gpp.sh                   # Build for GPP nodes
├── build_acc.sh                   # Build for ACC nodes
├── build_hybrid.sh                # Build unified
├── slurm_phase1_gpp.slurm        # Phase 1 job script
├── slurm_phase2_acc.slurm        # Phase 2 job script
├── slurm_phase3_hybrid.slurm     # Phase 3 job script
└── test_local_mpi.sh             # Local MPI testing
```

### Dependencies

**On MareNostrum 5:**
```bash
module load intel/oneapi       # Intel compilers + MPI
module load openmpi            # MPI implementation
module load cuda/12.x          # CUDA toolkit (for ACC nodes)
```

### Build Instructions

**Phase 1 (GPP nodes):**
```bash
./build_gpp.sh
# Creates: V1.6-mpi-cpu
```

**Phase 2 (ACC nodes):**
```bash
./build_acc.sh
# Creates: V1.6-mpi-gpu
```

**Phase 3 (Hybrid):**
```bash
./build_hybrid.sh
# Creates: V1.6-mpi-hybrid
```

### Deployment

**Phase 1:**
```bash
sbatch slurm_phase1_gpp.slurm
```

**Phase 2:**
```bash
sbatch slurm_phase2_acc.slurm
```

**Phase 3:**
```bash
sbatch slurm_phase3_hybrid.slurm
```

### Performance Targets

| Phase | Nodes | Resources | Expected Throughput |
|-------|-------|-----------|---------------------|
| 1 (GPP) | 10 | 1,120 CPU cores | 1.3-2.6 billion/sec |
| 2 (ACC) | 5 | 20 GPUs + 400 cores | 1-2 billion/sec |
| 3 (Hybrid) | 15 | 1,520 cores + 20 GPUs | 3-4 billion/sec |

### Design Principles

1. **Modular:** Each phase is standalone executable
2. **Testable:** Can test GPP/ACC independently
3. **Safe:** Phase 1 uses PROVEN V1.4b-openmp code
4. **Scalable:** Easy to add more nodes
5. **Monitorable:** Per-node progress tracking via MPI

### Work Distribution Strategy

**Embarrassingly Parallel:**
- Each MPI rank gets non-overlapping seed range
- Minimal inter-node communication
- Independent result files per rank
- Master rank (0) aggregates final results

**Seed Partitioning:**
```
Total range: [2^71 + offset, 2^71 + offset + N]
Ranks: R (total MPI ranks)

Rank i gets: [start + i*(N/R), start + (i+1)*(N/R))
```

**No conflicts, no synchronization needed!**

### Next Steps

1. ✅ Create node_config.hpp (hardware detection)
2. ✅ Create mpi_coordinator.hpp (work distribution)
3. ✅ Create worker_gpp.cpp (reuse V1.4b-openmp)
4. ⏳ Create worker_acc.cpp (CUDA + OpenMP)
5. ⏳ Create build scripts
6. ⏳ Create SLURM scripts
7. ⏳ Test locally with MPI
8. ⏳ Deploy to MareNostrum!
