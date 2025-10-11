# Collatz Conjecture - MareNostrum 5 Hackathon

**Goal:** Push tested limits from 2^71 using optimization techniques and MareNostrum 5 HPC

---

## Quick Start

```bash
cd secuncialTech
g++ -O3 -march=native -std=c++17 -o V1.0 V1.0.cpp
./V1.0 0 1000000  # Test 1M numbers from 2^71
```

---

## Current Status

**V1.0 Baseline:** 1,219,512 numbers/sec (820ms, 1M test)  
**V1.1 CTZ+Bitwise:** 1,736,111 numbers/sec (576ms, 1M test) - **1.42x faster**

See [`docs/PERFORMANCE.md`](docs/PERFORMANCE.md) for complete performance tracking table.  
See [`docs/V1.1_ANALYSIS.md`](docs/V1.1_ANALYSIS.md) for detailed V1.1 analysis.

---

## Project Structure

```
MN25/
├── secuncialTech/          # Source code versions
│   ├── V1.0.cpp           # Current: baseline
│   ├── build.sh           # Build script
│   └── bench.sh           # Benchmark script
├── docs/                   # Documentation
│   ├── PERFORMANCE.md     # Main performance tracking table
│   ├── STATUS.md          # Current status and next steps
│   └── README.md          # Detailed optimization roadmap
└── MareNostrumINFO/       # MareNostrum 5 specs
```

---

## Documentation

- **[docs/PERFORMANCE.md](docs/PERFORMANCE.md)** - Performance tracking table (THE MAIN TABLE)
- **[docs/STATUS.md](docs/STATUS.md)** - Current status and next optimization steps
- **[docs/README.md](docs/README.md)** - Detailed optimization roadmap and techniques

---

**Focus:** PURE PERFORMANCE
