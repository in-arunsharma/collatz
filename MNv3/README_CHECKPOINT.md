# V1.5-openmp-checkpoint: Overnight Run with Progress Tracking

## What This Does

- Uses **V1.5-openmp** (proven 137M nums/sec)
- Adds **checkpoint logging** every 1 billion numbers (~7 seconds)
- Logs progress to `checkpoint_overnight_job<JOBID>.log`
- If job crashes, you can see exactly how far it got

## Quick Start

```bash
cd /gpfs/projects/nct_352/nct01225/collatz/MNv3/

# Upload files
# (already done from local machine via scp)

# Build
./build_v15_checkpoint.sh

# Submit overnight job
sbatch run_overnight.slurm

# Check status
squeue -u nct01225

# Monitor progress (while running)
tail -f checkpoint_overnight_job*.log

# When complete, check results
cat v15_overnight_*.out
cat completion_*.txt
```

## Configuration

In `run_overnight.slurm`:
- **TOTAL_RANGE=400000000000** (400 billion numbers)
- At 137M/sec = **~48 minutes runtime**
- Adjust based on how long you want it to run

## Checkpoint Log Format

```
[CHECKPOINT] Completed: 1000000000 / 400000000000 (0.25%) | Throughput: 137.2M nums/sec | Elapsed: 7.3s
[CHECKPOINT] Completed: 2000000000 / 400000000000 (0.50%) | Throughput: 137.1M nums/sec | Elapsed: 14.6s
[CHECKPOINT] Completed: 3000000000 / 400000000000 (0.75%) | Throughput: 137.0M nums/sec | Elapsed: 21.9s
...
```

## Files Created

- `checkpoint_overnight_job<JOBID>.log` - Progress log (updated every ~7 seconds)
- `v15_overnight_<JOBID>.out` - Final results  
- `v15_overnight_<JOBID>.err` - stderr output (config, precompute, etc.)
- `completion_<JOBID>.txt` - Success marker (only if completed)
- Log files: `overflow_seeds_*.txt`, `fuse_seeds_*.txt`, `cold_*.txt`

## For Longer Runs

To run for full 8 hours:
```bash
# Edit run_overnight.slurm, change:
TOTAL_RANGE=4000000000000   # 4 trillion numbers (~8 hours)
```

## Recovery

If job crashes:
1. Check checkpoint log to see last completed point
2. Adjust START_OFFSET in run_overnight.slurm to resume
3. Resubmit

## Performance Expectations

- **Throughput**: ~137M nums/sec (proven baseline)
- **Checkpoint overhead**: Negligible (~0.1% performance impact)
- **400B numbers**: ~48 minutes
- **4T numbers**: ~8 hours (max for 8-hour time limit)
