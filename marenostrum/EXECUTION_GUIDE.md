# Phase 1 Execution Guide

## Quick Start (Monday Morning)

```bash
# 1. Navigate to marenostrum directory
cd ~/Desktop/MN25/marenostrum/

# 2. Make scripts executable
chmod +x *.sh

# 3. Integrate V1.4b core (see INTEGRATION_GUIDE.md)
#    Quick option: directly include V1.4b-openmp.cpp in worker_gpp_node.cpp

# 4. Build locally
./build_phase1_gpp.sh

# 5. Test with local MPI (2-4 ranks)
./test_local_mpi.sh 2

# 6. If test passes → Deploy to MareNostrum
```

---

## Local Development & Testing

### Prerequisites

**On your local machine:**
```bash
# Ubuntu/Debian
sudo apt install build-essential openmpi-bin libopenmpi-dev

# Verify MPI
mpirun --version
which mpic++
```

### Build Process

```bash
cd marenostrum/
./build_phase1_gpp.sh
```

**What it does:**
- Detects local vs MareNostrum environment
- Uses `mpic++` with OpenMP support
- Optimization: `-O3 -march=native -flto -fopenmp`
- Output: `collatz_mpi_gpp` executable

### Local Testing

**Test with 2 MPI ranks:**
```bash
./test_local_mpi.sh 2
```

**Test with 4 MPI ranks:**
```bash
./test_local_mpi.sh 4
```

**Manual testing:**
```bash
export OMP_NUM_THREADS=4
mpirun -np 2 ./collatz_mpi_gpp 0 1000000 test_run
```

**Expected output:**
```
[Rank 0] Master: Distributing work to 1 workers
[Rank 0] Work assigned: seeds 0 to 499999
[Rank 1] Worker: Received work range [500000, 1000000)
...
[Rank 0] Results collected from 1 workers
```

---

## MareNostrum 5 Deployment

### Step 1: Transfer Code

**From your local machine:**
```bash
# Compress the code
cd ~/Desktop/MN25
tar czf collatz_mn5.tar.gz marenostrum/ parallel/V1.4b-openmp.cpp

# Transfer to MareNostrum
scp collatz_mn5.tar.gz <username>@mn5.bsc.es:~/

# Example:
# scp collatz_mn5.tar.gz bsc12345@mn5.bsc.es:~/
```

### Step 2: Setup on MareNostrum

**SSH to MareNostrum:**
```bash
ssh <username>@mn5.bsc.es
```

**Extract and prepare:**
```bash
cd ~
tar xzf collatz_mn5.tar.gz
cd marenostrum/
chmod +x *.sh
```

### Step 3: Build on MareNostrum

```bash
# Build will auto-detect MareNostrum environment
./build_phase1_gpp.sh
```

**Modules loaded automatically:**
- `intel/2023.2` (Intel C++ compiler with OneAPI)
- `openmpi/4.1.5` (MPI implementation)

### Step 4: Submit Job

**Edit SLURM script if needed:**
```bash
nano slurm_phase1_gpp.slurm

# Adjust these variables:
START_SEED=0
END_SEED=10000000000  # 10 billion (adjust as needed)
OUTPUT_PREFIX="phase1_gpp_run1"
```

**Submit to queue:**
```bash
sbatch slurm_phase1_gpp.slurm
```

**Check job status:**
```bash
squeue -u $USER
```

**Monitor output in real-time:**
```bash
# Find job ID from squeue, then:
tail -f collatz_phase1_<JOBID>.out
```

---

## Understanding the Execution

### MPI + OpenMP Hierarchy

**Phase 1 (GPP nodes):**
```
MareNostrum Job
  ├─ Node 1 (MPI rank 0 - Master)
  │   └─ 112 OpenMP threads
  ├─ Node 2 (MPI rank 1 - Worker)
  │   └─ 112 OpenMP threads
  ├─ Node 3 (MPI rank 2 - Worker)
  │   └─ 112 OpenMP threads
  ├─ Node 4 (MPI rank 3 - Worker)
  │   └─ 112 OpenMP threads
  └─ Node 5 (MPI rank 4 - Worker)
      └─ 112 OpenMP threads

Total parallelism: 5 nodes × 112 threads = 560 threads
```

### Work Distribution

**Example: 10 billion numbers, 5 nodes**

```
Rank 0 (Master): Seeds 0 to 1,999,999,999       (2B numbers)
Rank 1:          Seeds 2,000,000,000 to 3,999,999,999
Rank 2:          Seeds 4,000,000,000 to 5,999,999,999
Rank 3:          Seeds 6,000,000,000 to 7,999,999,999
Rank 4:          Seeds 8,000,000,000 to 9,999,999,999
```

Each rank processes its range with 112 OpenMP threads in parallel.

### Performance Expectations

**Per-node (GPP):**
- 112 cores @ 2.0 GHz
- Expected: ~250M numbers/sec per node
- Walltime for 2B numbers: ~8 seconds

**5-node total:**
- Expected: **1.25B numbers/sec**
- Walltime for 10B numbers: ~8 seconds
- Plus overhead: MPI communication (~1-2 sec), I/O (~1 sec)

---

## Output Files

### JSON Metadata

Each run produces a JSON file: `{OUTPUT_PREFIX}_rank{N}.json`

**Example: `phase1_gpp_run1_rank0.json`**
```json
{
  "config": {
    "start_seed": 0,
    "end_seed": 1999999999,
    "num_threads": 112
  },
  "results": {
    "total_numbers": 2000000000,
    "elapsed_seconds": 8.234,
    "throughput_per_sec": 242876543,
    "max_steps": 596,
    "max_steps_seed": 63728127,
    "cycles_found": 23,
    "overflows_found": 5
  }
}
```

### Master Aggregated Results

Rank 0 (master) produces: `{OUTPUT_PREFIX}_global.json`

```json
{
  "global_results": {
    "total_numbers": 10000000000,
    "num_workers": 5,
    "elapsed_seconds": 8.5,
    "throughput_per_sec": 1176470588,
    "max_steps_global": 596,
    "max_steps_seed_global": 63728127,
    "cycles_found_total": 115,
    "overflows_found_total": 25
  }
}
```

---

## Troubleshooting

### Build Issues

**Problem:** `mpic++ not found`
```bash
# On local machine
sudo apt install libopenmpi-dev

# On MareNostrum (shouldn't happen, but if it does)
module load openmpi/4.1.5
```

**Problem:** OpenMP not enabled
```bash
# Check compiler flags
mpic++ --version
mpic++ -fopenmp test.cpp  # Should compile without errors
```

### Runtime Issues

**Problem:** MPI ranks not starting
```bash
# Check SLURM output
cat collatz_phase1_<JOBID>.err

# Verify executable exists
ls -lh collatz_mpi_gpp
```

**Problem:** Low performance
```bash
# Check OpenMP threads are spawning
export OMP_DISPLAY_ENV=TRUE
mpirun -np 2 ./collatz_mpi_gpp 0 1000 test

# Check CPU binding
export OMP_PROC_BIND=close
export OMP_PLACES=cores
```

**Problem:** Job pending in queue
```bash
# Check queue
squeue -u $USER

# Check job details
scontrol show job <JOBID>

# MareNostrum may have limits, adjust:
#SBATCH --nodes=2  # Try fewer nodes first
```

---

## Quick Reference Commands

### Local Testing
```bash
cd marenostrum/
./build_phase1_gpp.sh          # Build
./test_local_mpi.sh 2          # Test with 2 MPI ranks
mpirun -np 4 ./collatz_mpi_gpp 0 1000000 test  # Manual test
```

### MareNostrum Deployment
```bash
# Transfer
scp -r marenostrum/ <user>@mn5.bsc.es:~/collatz/

# SSH
ssh <user>@mn5.bsc.es

# Build & submit
cd ~/collatz/marenostrum/
./build_phase1_gpp.sh
sbatch slurm_phase1_gpp.slurm

# Monitor
squeue -u $USER
tail -f collatz_phase1_<JOBID>.out
```

### After Job Completes
```bash
# Check results
ls -lh *.json

# Download results to local machine (from local terminal)
scp <user>@mn5.bsc.es:~/collatz/marenostrum/*.json ./results/
```

---

## Monday Hackathon Timeline

**9:00 AM** - Integrate V1.4b core, build, local test  
**9:30 AM** - Transfer to MareNostrum, build there  
**10:00 AM** - Submit first test job (small range)  
**10:30 AM** - Analyze results, adjust parameters  
**11:00 AM** - Submit production run (large range)  
**12:00 PM** - Phase 1 complete, start Phase 2 (ACC nodes)  

---

**You're ready to scale to 1.25 billion numbers/sec on MareNostrum 5! 🚀**
