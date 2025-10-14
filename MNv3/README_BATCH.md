# Batch Overnight Run Strategy

## Problem
Your account only has `acc_debug` QOS which limits jobs to 2 hours.

## Solution
Submit **5 array jobs** that run sequentially, each processing 980 billion numbers (~2 hours).

## Coverage
- **Job 0**: Range [0, 980B)
- **Job 1**: Range [980B, 1960B)  
- **Job 2**: Range [1960B, 2940B)
- **Job 3**: Range [2940B, 3920B)
- **Job 4**: Range [3920B, 4900B)

**Total: 4.9 trillion numbers in ~10 hours**

## How to Use

### Submit All Jobs
```bash
cd /gpfs/projects/nct_352/nct01225/collatz/overnight13
sbatch run_batch_overnight.slurm
```

This submits 5 jobs at once. SLURM will schedule them (they may run in parallel if nodes are available, or sequentially if not).

### Monitor Progress
```bash
./monitor_batch.sh
```

Shows:
- Which jobs are running/completed
- Latest checkpoint from each job
- Total progress

### Check Individual Job
```bash
# See live progress from job 0
tail -f checkpoint_batch0_job*.log

# Check output from job 2
cat v15_batch_*_2.out
```

## Array Job Behavior

**Key difference from sequential submission:**
- All 5 jobs are submitted simultaneously
- SLURM scheduler decides when each runs
- If nodes available: Multiple jobs run in PARALLEL (faster!)
- If nodes busy: Jobs queue and run when resources free
- Each job is independent - no dependencies

## Expected Timeline

**If running in parallel** (5 nodes available):
- All 5 jobs finish in ~2 hours
- Total throughput: 685M nums/sec (5 × 137M)

**If running sequentially** (1 node available):
- Jobs run one after another
- Total time: ~10 hours
- Throughput: 137M nums/sec per job

## Files Generated

Per job:
- `v15_batch_<JOBID>_<0-4>.out` - Results
- `v15_batch_<JOBID>_<0-4>.err` - stderr/config
- `checkpoint_batch<0-4>_job<JOBID>.log` - Progress log
- `completion_batch<0-4>_<JOBID>.txt` - Completion marker

## Advantages

✅ Works with 2-hour QOS limit
✅ Can run in parallel if nodes available  
✅ Each job independent (one fails, others continue)
✅ Easy to track progress
✅ Can submit and sleep!

## To Stop All Jobs

```bash
scancel <JOBID>
```

Where JOBID is the main array job ID (shown when you submit).
