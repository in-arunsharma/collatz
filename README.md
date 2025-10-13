# Collatz Conjecture - MareNostrum 5 Hackathon

**Goal:** Push tested limits from 2^71 using optimization techniques and MareNostrum 5 HPC

**Status:** ✅ **PRODUCTION READY** - V1.4b-openmp parallelization complete, tested locally, ready for MareNostrum deployment!

---

## Quick Start

### Sequential Version (V1.4b)
```bash
cd secuncialTech
./build_v1.4b.sh
./V1.4b 0 1000000 --tag test  # 2.38M nums/sec
```

### OpenMP Parallel Version (V1.4b-openmp) - **NEW!**
```bash
cd parallel
./build_v1.4b-openmp.sh
OMP_NUM_THREADS=4 ./V1.4b-openmp 0 1000000 --tag test  # 5.40M nums/sec @ 4 cores
```

---

## Current Status

**Performance Evolution:**
- **V1.0 Baseline:** 1.22M nums/sec
- **V1.3d Sequential:** 2.44M nums/sec (2.00× speedup)
- **V1.4b Hot/Cold Queue:** 2.38M nums/sec (sequential, production logging)
- **V1.4b-openmp (4 cores):** 5.40M nums/sec (2.33× speedup, 58% efficiency)
- **Target (80 cores):** ~192M nums/sec for MareNostrum 5 hackathon

**Latest Achievement:** ✅ OpenMP parallelization with:
- Per-thread state (no false sharing)
- Deterministic seed mapping (reproducible results)
- Thread-safe logging (per-thread files)
- Post-parallel cold queue processing
- NUMA-aware affinity settings

See [`docs/PERFORMANCE.md`](docs/PERFORMANCE.md) for complete performance tracking table.

---

## Project Structure

```
MN25/
├── secuncialTech/          # Sequential implementations
│   ├── V1.4b.cpp          # Production: Hot/Cold queue (2.38M nums/sec)
│   ├── build_v1.4b.sh     # Build script
│   └── perf_v1.4b.txt     # Performance analysis
├── parallel/               # Parallel implementations - **NEW!**
│   ├── V1.4b-openmp.cpp   # ✅ OpenMP parallelization (5.40M @ 4 cores)
│   ├── build_v1.4b-openmp.sh
│   ├── marenostrum_single_node.slurm  # 80-core SLURM job
│   ├── marenostrum_array_job.slurm    # Multi-node array job
│   └── README.md          # Deployment guide
├── docs/                   # Documentation
│   ├── PERFORMANCE.md     # Performance tracking table
│   ├── V1.4b_ANALYSIS.md  # Hot/cold queue architecture
│   └── READY_FOR_PARALLELIZATION.md  # Friend's review
└── MareNostrumINFO/       # MareNostrum 5 specs (80 cores, 4480 GPUs)
```

---

## Documentation

- **[parallel/README.md](parallel/README.md)** - ✅ **OpenMP deployment guide (MareNostrum 5)**
- **[docs/PERFORMANCE.md](docs/PERFORMANCE.md)** - Performance tracking table
- **[docs/READY_FOR_PARALLELIZATION.md](docs/READY_FOR_PARALLELIZATION.md)** - Friend's cluster-ready checklist
- **[docs/V1.4b_ANALYSIS.md](docs/V1.4b_ANALYSIS.md)** - Hot/cold queue architecture

---

**Focus:** PURE PERFORMANCE
