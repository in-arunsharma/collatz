# 🚀 MARENOSTRUM FILESYSTEM GUIDE

## 📁 Where to Put Your Code & Data

### ❌ DON'T USE: `/gpfs/home/nct01/nct01225/`
- **Why:** Not designed for running jobs (slow, I/O bottleneck)
- **Use for:** Source code storage only
- **Quota:** Limited per user

### ✅ USE THIS: `/gpfs/scratch/<GROUP>/nct01225/`
- **Why:** Designed for job execution (fast, optimized for I/O)
- **Use for:** Running jobs, temporary output files
- **Quota:** Group-based (larger than home)
- **Example:** `/gpfs/scratch/bsc01/nct01225/` (replace `bsc01` with your group)

### 🚀 FASTEST: `$TMPDIR` (Local SSD on compute node)
- **Why:** Local NVMe SSD = ultra-fast I/O
- **Location:** `/scratch/tmp/$SLURM_JOB_ID`
- **Auto-set:** `$TMPDIR` environment variable
- **Use for:** Temporary computation files during job execution
- **WARNING:** Deleted automatically when job finishes!

---

## 🔧 Updated Deployment Steps

### Step 1: Transfer Code
```bash
# From your local machine
scp collatz_mn5.tar.gz nct01225@glogin1.bsc.es:~/
```

### Step 2: SSH and Setup in SCRATCH (not home!)
```bash
ssh nct01225@glogin1.bsc.es

# Check your group (replace bsc01 with your actual group)
groups

# Create your scratch working directory
# IMPORTANT: Replace 'bscXX' with your actual group!
mkdir -p /gpfs/scratch/bscXX/nct01225/collatz
cd /gpfs/scratch/bscXX/nct01225/collatz

# Extract code here
tar xzf ~/collatz_mn5.tar.gz
cd marenostrum/

# Check quota
bsc_quota
```

### Step 3: Build in Scratch Directory
```bash
# Make scripts executable
chmod +x *.sh

# Build
./build_phase1_gpp.sh
```

### Step 4: Submit Job (already optimized for scratch)
```bash
sbatch slurm_phase1_gpp.slurm
```

---

## 📊 Filesystem Performance Comparison

| Filesystem | Speed | Use Case | Auto-deleted? |
|-----------|-------|----------|---------------|
| `/gpfs/home` | Slow | Code storage | No |
| `/gpfs/scratch` | Fast | Job execution | No |
| `$TMPDIR` (local SSD) | **Fastest** | Per-node temp files | **Yes** |

---

## 🎯 Optimal Setup for Your Collatz Job

### Where Each Part Goes:

1. **Source code:** `/gpfs/scratch/bscXX/nct01225/collatz/marenostrum/`
   - This is where you extract and build
   - SLURM runs from here

2. **Output JSON files:** `/gpfs/scratch/bscXX/nct01225/collatz/marenostrum/`
   - Results written here by default
   - Persistent (not deleted after job)

3. **Temporary computation:** `$TMPDIR` (optional optimization)
   - For memo tables, intermediate data
   - Ultra-fast local SSD
   - Auto-cleaned after job

---

## 🔧 Finding Your Group

```bash
# After SSH to MareNostrum
groups

# Example output:
# nct01225 bsc01 other_groups
#          ^^^^^ This is your group!

# Then use:
cd /gpfs/scratch/bsc01/nct01225/
```

---

## ⚠️ Important Notes

1. **Job submission location matters:**
   - SLURM's `$SLURM_SUBMIT_DIR` = where you ran `sbatch`
   - Make sure you're in `/gpfs/scratch/...` not `/gpfs/home`

2. **Local SSD cleanup:**
   - `$TMPDIR` is **automatically deleted** when job finishes
   - Don't put final results there!
   - Good for: memo tables, intermediate computation

3. **Quota checking:**
   ```bash
   bsc_quota  # Check your disk usage
   ```

4. **Backups:**
   - `/gpfs/scratch` is **NOT backed up**
   - Copy important results back to your local machine!

---

## 🚀 Complete Deployment Workflow

```bash
# 1. Transfer from local machine
scp collatz_mn5.tar.gz nct01225@glogin1.bsc.es:~/

# 2. SSH
ssh nct01225@glogin1.bsc.es

# 3. Check your group
groups  # Note the group name (e.g., bsc01)

# 4. Create scratch workspace (REPLACE bscXX!)
mkdir -p /gpfs/scratch/bscXX/nct01225/collatz
cd /gpfs/scratch/bscXX/nct01225/collatz

# 5. Extract
tar xzf ~/collatz_mn5.tar.gz
cd marenostrum/

# 6. Build
chmod +x *.sh
./build_phase1_gpp.sh

# 7. Submit job (from scratch directory!)
sbatch slurm_phase1_gpp.slurm

# 8. Monitor
squeue -u $USER
tail -f collatz_phase1_*.out

# 9. After job completes, download results
# (from your local machine)
scp nct01225@glogin1.bsc.es:/gpfs/scratch/bscXX/nct01225/collatz/marenostrum/*.json ./results/
```

---

## 📈 Performance Optimization (Future)

For maximum performance, you can use `$TMPDIR` for hot data:

```bash
# In your worker code (future optimization):
# 1. Copy memo table to $TMPDIR at job start
# 2. Process using local SSD (ultra-fast)
# 3. Write results back to /gpfs/scratch
```

**Current setup:** Already optimized to run from `/gpfs/scratch` ✅

---

## 🎯 Quick Reference

**Your working directory:**
```
/gpfs/scratch/<YOUR_GROUP>/nct01225/collatz/marenostrum/
```

**Check quota:**
```bash
bsc_quota
```

**Check group:**
```bash
groups
```

**Filesystem hierarchy:**
```
/gpfs/scratch/<GROUP>/<USERNAME>/
    └── collatz/
        └── marenostrum/
            ├── collatz_mpi_gpp (executable)
            ├── *.json (output files)
            └── *.out, *.err (SLURM logs)
```

---

**REMEMBER: Run from `/gpfs/scratch/<GROUP>/`, NOT from `/gpfs/home`!** 🚀
