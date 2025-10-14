# MNv2 Quick Reference

## 🚀 Deploy in 5 Minutes

```bash
# 1. Package
cd ~/Desktop/MN25 && tar czf mnv2.tar.gz MNv2/

# 2. Transfer
scp mnv2.tar.gz nct01225@glogin1.bsc.es:~/

# 3. Setup on MareNostrum
ssh nct01225@glogin1.bsc.es
cd /gpfs/projects/nct_352/nct01225/collatz  # Projects dir (backed up!)
tar xzf ~/mnv2.tar.gz && cd MNv2/

# Create output directory for job logs (in scratch)
mkdir -p /gpfs/scratch/nct_352/nct01225/collatz_output

# 4. Build
./build_gpp.sh

# 5. Submit
sbatch slurm_gpp.slurm
```

---

## ⚙️ Configuration

### Edit `config.hpp`
```cpp
constexpr int GPP_NODES = 2;     // Your node allocation
constexpr int ACC_NODES = 0;     // Phase 2
```

### Edit `slurm_gpp.slurm`
```bash
#SBATCH --nodes=2          # Must match GPP_NODES
#SBATCH --ntasks=3         # nodes + 1
#SBATCH --time=00:15:00    # Walltime

START_OFFSET=0
COUNT=1000000000           # 1B seeds (~2 seconds)
```

---

## 📊 Performance Expectations

| Nodes | Cores | Throughput | Time (1B seeds) |
|-------|-------|------------|-----------------|
| 2     | 224   | 500M/s     | ~2 sec          |
| 5     | 560   | 1.2B/s     | ~0.8 sec        |
| 10    | 1,120 | 2.5B/s     | ~0.4 sec        |

---

## 🐛 Troubleshooting

### Build fails
```bash
module load intel/2023.2 impi/2021.10
./build_gpp.sh
```

### Job fails
```bash
squeue -u $USER       # Check status
sacct -j <JOBID>      # Check errors
sinfo -p gpp          # Check node availability
```

### Low performance
- ✅ Check thread count (should be 112 for GPP)
- ✅ Verify `OMP_PROC_BIND=close` in slurm script
- ✅ Use large workload (1B+ seeds minimum)

---

## 📝 Monitor Job

```bash
# Watch output
tail -f mnv2_gpp_<JOBID>.out

# Check queue
squeue -u $USER

# View results
sacct -j <JOBID> --format=JobID,State,ExitCode,Elapsed,MaxRSS
```

---

## ✅ Success Checklist

- [ ] Transferred to MareNostrum
- [ ] Built successfully
- [ ] Job submitted
- [ ] Results show ~500M nums/sec (2 nodes)
- [ ] Ready to scale!

---

## 📚 Documentation

- **README.md** - Architecture overview
- **DEPLOY.md** - Detailed deployment guide
- **SUMMARY.md** - Technical implementation details
- **QUICKREF.md** - This file

---

## 🎯 Quick Commands

```bash
# Local test
./test_local.sh

# Build
./build_gpp.sh

# Submit
sbatch slurm_gpp.slurm

# Monitor
squeue -u $USER
tail -f mnv2_gpp_*.out

# Cancel job
scancel <JOBID>

# Check allocation
sshare -U
```

---

## 💡 Tips

1. **Start small:** 2 nodes, 1B seeds, 15 min walltime
2. **Verify correctness:** Check results match expected
3. **Scale gradually:** 2 → 5 → 10 nodes
4. **Monitor efficiency:** Should be >80%
5. **Filesystem usage:**
   - Source code → `/gpfs/projects/nct_352/$USER/` (backed up)
   - Job outputs → `/gpfs/scratch/nct_352/$USER/` (temporary)
   - Temp files → `$TMPDIR` (local SSD, fast I/O)

---

## 🆘 Help

- MareNostrum 5 docs: https://www.bsc.es/marenostrum/marenostrum5
- Check DEPLOY.md for common issues
- Test locally with `./test_local.sh` first

**Good luck! 🚀**
