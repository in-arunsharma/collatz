# MareNostrum 5 Filesystem Best Practices for MNv2

## 📂 Correct Directory Structure

### ✅ What We're Using Now (Corrected)

```
/gpfs/projects/nct_352/nct01225/collatz/MNv2/
├── collatz_mpi_gpp           # Executable (built here)
├── collatz_mpi_gpp.cpp       # Source code
├── config.hpp                # Configuration
├── build_gpp.sh              # Build script
├── slurm_gpp.slurm           # Job script
└── *.md                      # Documentation

/gpfs/scratch/nct_352/nct01225/collatz_output/
├── mnv2_gpp_12345.out        # Job stdout
├── mnv2_gpp_12345.err        # Job stderr
└── output/                   # Overflow/fuse logs (optional)

$TMPDIR (auto-set by SLURM, e.g., /scratch/tmp/12345/)
└── <temp files during job execution>
```

---

## 🎯 MareNostrum 5 Filesystem Types

### 1. `/gpfs/projects/<GROUP>/<USER>/` - Project Space
**Purpose:** Long-term project data, source code, executables  
**Backed up:** ✅ Yes (every 3-4 days)  
**Quota:** Group-based (shared with team)  
**Access:** From all nodes (login, compute)  
**Speed:** Good for source code and small files

**✅ Use for:**
- Source code (`.cpp`, `.hpp`, `.cu`)
- Compiled executables
- Build scripts (`.sh`)
- Job scripts (`.slurm`)
- Documentation (`.md`)

**❌ Don't use for:**
- Large temporary files
- Job output logs (use `/gpfs/scratch`)
- High-frequency I/O during job execution (use `$TMPDIR`)

### 2. `/gpfs/scratch/<GROUP>/<USER>/` - Temporary Scratch
**Purpose:** Temporary job outputs, large intermediate files  
**Backed up:** ❌ No backups  
**Quota:** Group-based  
**Access:** From all nodes  
**Speed:** Good for large files shared across nodes

**✅ Use for:**
- Job stdout/stderr logs (`.out`, `.err`)
- Overflow/fuse/cycle logs
- Large result files
- Checkpoint files
- Data that needs to persist across multiple jobs

**❌ Don't use for:**
- Source code (no backup!)
- Executables (no backup!)
- Anything you want to keep long-term

### 3. `$TMPDIR` - Local SSD (Per-Node)
**Purpose:** Fast temporary storage during job execution  
**Backed up:** ❌ No (auto-deleted after job)  
**Quota:** Per-node limit (varies by partition)  
**Access:** Only from that specific compute node  
**Speed:** ⚡ Very fast (NVMe SSD)

**✅ Use for:**
- Intermediate computation files
- Temporary checkpoints
- High-frequency I/O operations
- Files that don't need to survive job completion

**❌ Don't use for:**
- Anything you want to keep after job ends
- Files that need to be accessed from other nodes
- Large files (limited space per node)

### 4. `/gpfs/home/<USER>/` - Home Directory
**Purpose:** Personal configuration, small scripts  
**Backed up:** ✅ Yes (daily)  
**Quota:** Per-user (default: small)  
**Access:** From all nodes  
**Speed:** Slow (avoid for jobs)

**✅ Use for:**
- Personal `.bashrc`, config files
- Small utility scripts
- Personal notes

**❌ Don't use for:**
- Running jobs (too slow!)
- Source code (use `/gpfs/projects`)
- Job outputs (use `/gpfs/scratch`)

---

## 🚀 MNv2 Deployment Commands (Corrected)

### Step-by-Step

```bash
# 1. Local: Package code
cd ~/Desktop/MN25
tar czf mnv2.tar.gz MNv2/

# 2. Local: Transfer to MareNostrum
scp mnv2.tar.gz nct01225@glogin1.bsc.es:~/

# 3. MareNostrum: Extract in PROJECTS directory
ssh nct01225@glogin1.bsc.es
cd /gpfs/projects/nct_352/nct01225/
mkdir -p collatz
cd collatz
tar xzf ~/mnv2.tar.gz
cd MNv2/

# 4. MareNostrum: Create output directory in SCRATCH
mkdir -p /gpfs/scratch/nct_352/nct01225/collatz_output

# 5. MareNostrum: Build
./build_gpp.sh

# 6. MareNostrum: Submit
sbatch slurm_gpp.slurm

# 7. MareNostrum: Monitor
tail -f /gpfs/scratch/nct_352/nct01225/collatz_output/mnv2_gpp_*.out
```

---

## 📊 Why This Structure?

### Projects Directory (Source Code)
✅ **Backed up** - Your source code is safe  
✅ **Persistent** - Won't be deleted  
✅ **Shared** - Team members can access  
✅ **Version controlled** - Good for git repositories

### Scratch Directory (Job Outputs)
✅ **No quota pressure** on projects  
✅ **Large files OK** - Better quota limits  
✅ **Temporary** - Auto-cleanup policies  
✅ **Visible across nodes** - Can check from login node

### $TMPDIR (Job-Local Temp)
✅ **Fastest I/O** - Local NVMe SSD  
✅ **No network overhead** - Direct disk access  
✅ **Auto-cleanup** - Deleted after job (no manual cleanup needed)  
✅ **Per-job isolation** - No conflicts with other jobs

---

## ⚠️ Common Mistakes (FIXED in MNv2)

| ❌ Wrong | ✅ Correct | Why |
|---------|-----------|-----|
| Source in `/gpfs/scratch` | Source in `/gpfs/projects` | No backup in scratch! |
| Build in `/gpfs/home` | Build in `/gpfs/projects` | Home is too slow |
| Logs in `/gpfs/projects` | Logs in `/gpfs/scratch` | Don't waste project quota |
| Use `/tmp` for temp files | Use `$TMPDIR` | `/tmp` is NFS (slow) |
| Run from home directory | Run from projects | Home has strict quotas |

---

## 🔍 Check Your Quotas

```bash
# On MareNostrum
bsc_quota

# Example output:
# Filesystem     Size  Used  Avail  Use%  Path
# projects       100G   25G    75G   25%   /gpfs/projects/nct_352
# scratch        500G  120G   380G   24%   /gpfs/scratch/nct_352
# home            50G   10G    40G   20%   /gpfs/home/nct01225
```

---

## 📝 Updated SLURM Script Paths

The `slurm_gpp.slurm` now correctly uses:

```bash
#SBATCH --output=/gpfs/scratch/nct_352/nct01225/collatz_output/mnv2_gpp_%j.out
#SBATCH --error=/gpfs/scratch/nct_352/nct01225/collatz_output/mnv2_gpp_%j.err

# Create output directory (idempotent)
mkdir -p /gpfs/scratch/nct_352/nct01225/collatz_output

# Run executable from projects directory
./collatz_mpi_gpp 0 1000000000 mnv2_gpp_${SLURM_JOB_ID}
```

---

## ✅ Verification Checklist

Before deploying, verify:

- [x] Source code in `/gpfs/projects/nct_352/nct01225/collatz/MNv2/`
- [x] Executable built in projects directory
- [x] Output directory created: `/gpfs/scratch/nct_352/nct01225/collatz_output/`
- [x] SLURM script points to correct output directory
- [x] No files in `/gpfs/home` (except tarball)
- [x] No hardcoded `/tmp` paths (use `$TMPDIR` if needed)

---

## 🎯 Performance Impact

Using correct filesystems improves performance:

| Filesystem | Latency | Throughput | Best for |
|------------|---------|------------|----------|
| Local SSD (`$TMPDIR`) | **~100 μs** | **3+ GB/s** | Temp files during job |
| GPFS (`/gpfs/*`) | ~1 ms | ~10 GB/s | Shared data, logs |
| NFS (`/tmp`, `/`) | ~10 ms | ~100 MB/s | OS only, avoid! |

**Bottom line:** Use local SSD for hot I/O, GPFS for shared data!

---

## 📚 References

- [MareNostrum 5 Filesystems](https://www.bsc.es/user-support/mn5.php#filesystems)
- [GPFS Best Practices](https://www.ibm.com/docs/en/gpfs)
- Your account: `nct01225@glogin1.bsc.es`
- Your group: `nct_352`

**All documentation has been updated!** ✅
