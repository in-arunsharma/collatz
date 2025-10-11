# Perf Tools Guide - Instructions and Cycles Analysis

## Overview

`perf` has different tools for different purposes:

---

## 1. perf stat - PERFORMANCE COUNTERS (for your table)

**Use this for:** Instructions, cycles, IPC, branch misses, cache misses

### Basic command:
```bash
perf stat ./V1.1 0 1000000
```

### Detailed counters:
```bash
perf stat -d ./V1.1 0 1000000
```

### What you get:
```
Performance counter stats for './V1.1 0 1000000':

   655.90 msec task-clock
2,711,287,554 cpu_core/cycles/                # 4.134 GHz
8,443,523,416 cpu_core/instructions/          # IPC: 3.11
1,571,593,913 cpu_core/branches/
    1,029,218 cpu_core/branch-misses/         # 0.065%
      
   13,104,593 L1-dcache-loads
       94,226 L1-dcache-load-misses
       54,327 LLC-loads
       11,370 LLC-load-misses
```

**This is what you put in your PERFORMANCE.md table!**

---

## 2. perf record + perf report - HOTSPOT ANALYSIS

**Use this for:** Finding which functions are slow

### Step 1: Record
```bash
perf record ./V1.1 0 100000
# Creates perf.data file
```

### Step 2: View report
```bash
perf report
# Interactive TUI

# OR text output:
perf report --stdio
```

### What you get:
```
# Overhead  Function
  94.74%    compute_collatz()
   2.34%    main()
   ...
```

**Interpretation:**
- 94.74% of CPU time in `compute_collatz` = good (that's our hot loop)
- 2.34% in main = good (minimal overhead)
- If you see high % in print functions = bad (I/O bottleneck)

---

## 3. perf annotate - ASSEMBLY ANALYSIS

**Use this for:** Seeing which CPU instructions are slow

### Command:
```bash
perf annotate compute_collatz
```

### What you get:
Assembly code with % time per instruction:
```
  2.45 │   shr    %cl,%rdx          # Right shift (from >>)
  0.89 │   add    $0x1,%r12         # Increment steps
 15.23 │   imul   $0x3,%rdx,%rax    # 3*n (expensive!)
  8.45 │   add    $0x1,%rax         # +1
  3.12 │   tzcnt  %rax,%rcx         # Count trailing zeros (CTZ)
```

**This shows which operations take most cycles!**

---

## 4. perf list - AVAILABLE EVENTS

```bash
perf list
```

Shows all hardware events you can measure:
- cycles, instructions, branches, cache-misses, etc.

---

## Quick Reference for Your Use Case

### For PERFORMANCE.md table:
```bash
perf stat -d ./V1.1 0 1000000 2>&1 | tee perf_v1.1.txt
```

Extract:
- instructions
- cycles  
- IPC (instructions ÷ cycles)
- branch-misses
- L1D/LLC misses

### For finding bottlenecks:
```bash
perf record ./V1.1 0 1000000
perf report --stdio | head -20
```

Look for:
- High % in unexpected functions = problem
- High % in compute_collatz = expected (hot loop)

### For instruction-level optimization:
```bash
perf record ./V1.1 0 1000000
perf annotate compute_collatz
```

Look for:
- High % on multiply/divide = consider optimization
- High % on branches = predict better or eliminate

---

## Advanced: Specific Counter Events

### Custom events:
```bash
perf stat -e cycles,instructions,branches,branch-misses,cache-misses,L1-dcache-load-misses ./V1.1 0 1000000
```

### CPU-specific events:
```bash
# List available:
perf list | grep -i cache

# Measure specific:
perf stat -e L1-dcache-loads,L1-dcache-load-misses,LLC-loads,LLC-load-misses ./V1.1 0 1000000
```

---

## Summary Table

| Tool | Purpose | Output | When to Use |
|------|---------|--------|-------------|
| `perf stat` | **Counters** | Instructions, cycles, IPC, misses | **For performance table** |
| `perf report` | **Hotspots** | % time per function | Find slow functions |
| `perf annotate` | **Assembly** | Cycles per instruction | Micro-optimization |
| `perf list` | **Events** | Available counters | Discover what to measure |

---

## Your Workflow

1. **Benchmark new version:**
   ```bash
   perf stat -d ./V1.X 0 1000000 2>&1 | tee perf_v1.X.txt
   ```

2. **Update PERFORMANCE.md** with numbers from perf stat

3. **If performance is bad, find hotspot:**
   ```bash
   perf record ./V1.X 0 1000000
   perf report --stdio
   ```

4. **If specific function is slow, analyze assembly:**
   ```bash
   perf annotate function_name
   ```

5. **Optimize and repeat!**
