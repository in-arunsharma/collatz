# MareNostrum 5 Technical Specifications
## EuroHPC Supercomputer - BSC-CNS Barcelona

### System Overview
- **Peak Performance:** 314 PFlops total
- **Architecture:** Bull Sequana XH3000 + Lenovo ThinkSystem hybrid
- **Location:** Barcelona Supercomputing Center (BSC-CNS)
- **Classification:** Pre-exascale EuroHPC supercomputer

### System Partitions

#### 1. General Purpose Partition (GPP)
**Node Count:** 6,408 standard nodes + 72 HBM nodes

**Standard Node Configuration:**
- **CPUs:** 2× Intel Sapphire Rapids 8480+ @ 2.0GHz
- **Cores:** 56 cores per CPU × 2 = 112 cores per node
- **Memory:** 256GB DDR5 (216 nodes have 1024GB)
- **Storage:** 960GB NVMe local storage (/scratch)
- **Network:** 1× NDR200 shared by 2 nodes (100Gb/s per node)

**HBM Nodes (72 nodes):**
- **CPUs:** Intel Sapphire Rapids 03H-LC @ 1.9GHz
- **Cores:** 112 cores per node
- **Memory:** 128GB HBM (High Bandwidth Memory)
- **Bandwidth:** 2TB/s memory bandwidth per node

**GPP Performance:**
- **Peak:** 45.9 PFlops
- **Total Cores:** ~717,000 cores
- **Total Memory:** ~1.6 PB DDR5 + ~9TB HBM

#### 2. Accelerated Partition (ACC) - **TARGET FOR COLLATZ PROJECT**
**Node Count:** 1,120 nodes

**Node Configuration:**
- **CPUs:** 2× Intel Sapphire Rapids 8460Y+ @ 2.3GHz  
- **CPU Cores:** 40 cores per CPU × 2 = 80 cores per node
- **System Memory:** 512GB DDR5 per node
- **GPUs:** 4× NVIDIA Hopper GPUs per node
- **GPU Memory:** 64GB HBM2 per GPU × 4 = 256GB GPU memory per node
- **Local Storage:** 480GB NVMe (/scratch)
- **Network:** 4× NDR200 (800Gb/s total bandwidth per node)

**ACC Performance:**
- **Peak:** 260 PFlops
- **Total GPUs:** 4,480 NVIDIA Hopper GPUs
- **Total GPU Memory:** ~1.1 PB HBM2
- **Total System Memory:** ~560 TB DDR5

#### 3. General Purpose - Next Generation Partition
- **Status:** Based on NVIDIA GRACE CPU
- **Details:** Limited information available

#### 4. Accelerated - Next Generation Partition  
- **Status:** Not fully defined
- **Timeline:** More information expected in coming months

### Network Architecture

#### GPP Network Topology:
- **Type:** Fat-tree architecture
- **Island Size:** 2,160 nodes per island (full fat-tree, no contention)
- **Inter-island:** 2/3 contention ratio
- **Bandwidth:** NDR200 (shared between 2 nodes)

#### ACC Network Topology:
- **Type:** Fat-tree architecture  
- **Island Size:** 160 nodes per island (full fat-tree, no contention)
- **Inter-island:** 1/2 contention ratio
- **Bandwidth:** 4× NDR200 per node (800Gb/s total)

### Storage System
- **Capacity:** 248PB net capacity
- **Technology:** SSD/Flash + Hard disks hybrid
- **Read Performance:** 1.6TB/s aggregate
- **Write Performance:** 1.2TB/s aggregate
- **Filesystem:** IBM Storage Scale (parallel filesystem)
- **Archive:** 402PB tape-based long-term storage

### Software Environment

#### Compilers and Development:
- **OS:** Red Hat Enterprise Server
- **Intel Tools:**
  - Intel OneAPI toolkit
  - C/C++/Fortran compilers
  - Intel MKL (Math Kernel Library)
  - Intel MPI
  - Intel Trace Analyzer and Collector
- **NVIDIA Tools:**
  - NVIDIA HPC SDK
  - NVIDIA CUDA Toolkit
- **Other:**
  - OpenMPI
  - DDT parallel debugger

#### HPC and Performance Tools:
- **Scheduler:** Slurm batch scheduler
- **Performance:** BSC performance tools
- **Energy:** EAR (Energy management framework for HPC)

### Hardware Specifications Deep Dive

#### Intel Sapphire Rapids 8480+ (GPP)
- **Architecture:** x86-64, 4th Gen Intel Xeon Scalable
- **Cores:** 56 cores per socket
- **Base Frequency:** 2.0GHz
- **Memory:** DDR5 support
- **Features:** AVX-512, Intel AMX (Advanced Matrix Extensions)

#### Intel Sapphire Rapids 8460Y+ (ACC)
- **Architecture:** x86-64, 4th Gen Intel Xeon Scalable  
- **Cores:** 40 cores per socket
- **Base Frequency:** 2.3GHz
- **Memory:** DDR5 support
- **Features:** Optimized for GPU workloads

#### NVIDIA Hopper GPUs (ACC Partition)
**Assumed H100 specifications:**
- **Architecture:** Hopper (4th Gen Tensor Cores)
- **Memory:** 64GB HBM2 per GPU
- **Memory Bandwidth:** ~3TB/s per GPU
- **Compute Units:** ~16,896 CUDA cores per GPU
- **Tensor Performance:** Enhanced for AI/ML workloads
- **Integer Performance:** Excellent for computational workloads
- **Interconnect:** NVLink 4.0 for GPU-to-GPU communication

### Collatz Project Resource Mapping

#### Recommended Resource Allocation:
**Target Partition:** Accelerated (ACC)
- **Rationale:** NVIDIA Hopper GPUs ideal for integer-heavy Collatz computations

**Node Allocation Strategy:**
- **Phase 1:** 4-10 nodes (16-40 GPUs)
- **Phase 2:** 10-20 nodes (40-80 GPUs)  
- **Phase 3:** 20-50 nodes (80-200 GPUs)
- **Maximum:** 50-100 nodes (200-400 GPUs)

#### Performance Projections:
**Single Hopper GPU:**
- **Theoretical:** ~1M Collatz evaluations/second
- **128-bit arithmetic overhead:** ~100K-500K evaluations/second
- **Memory bound operations:** Dependent on sequence lengths

**Multi-GPU Scaling:**
- **40 GPUs:** 4-20M evaluations/second
- **200 GPUs:** 20-100M evaluations/second
- **400 GPUs:** 40-200M evaluations/second

#### Memory Requirements:
**Per GPU (64GB HBM2):**
- **Input ranges:** ~1GB (millions of starting numbers)
- **Intermediate calculations:** Variable (depends on sequence lengths)
- **Results storage:** ~1GB (step counts and anomalies)
- **Working space:** ~62GB available

**Per Node (512GB DDR5):**
- **Coordination data:** Node-level work distribution
- **Result aggregation:** Cross-GPU result collection
- **Backup storage:** Checkpoint data for fault tolerance

### Network Performance for Collatz:
**Intra-node (4 GPUs):**
- **NVLink:** Ultra-high bandwidth for GPU-to-GPU
- **PCIe Gen5:** CPU-GPU communication

**Inter-node:**
- **NDR200:** 800Gb/s per node - excellent for work distribution
- **Fat-tree:** Low latency for coordination messages
- **Embarrassingly parallel:** Minimal inter-node communication needed

### Development Environment Setup:
```bash
# Expected module system
module load cuda/12.x
module load intel/oneapi
module load mpi/intel

# Compilation example
nvcc -O3 -arch=sm_90 collatz.cu -o collatz_gpu
mpicc -O3 -fopenmp collatz_mpi.c -o collatz_mpi
```

### Resource Request Justification:
**Why Collatz Deserves Significant ACC Resources:**
1. **Perfect fit:** Integer arithmetic ideal for Hopper GPUs
2. **Embarrassingly parallel:** Excellent scaling characteristics
3. **Time-bound:** Hackathon requires maximum parallelism
4. **Educational:** Demonstrates supercomputer capabilities
5. **Scientific:** Advancing mathematical knowledge
6. **Technical showcase:** Highlights MareNostrum 5 power

### Comparison with Other Supercomputers:
- **Frontier (ORNL):** ~1.1 exaflops, AMD-based
- **Fugaku (RIKEN):** ~442 petaflops, ARM-based  
- **Summit (ORNL):** ~200 petaflops, IBM Power + NVIDIA V100
- **MareNostrum 5:** 314 petaflops, Intel + NVIDIA Hopper

**MareNostrum 5 Advantages for Collatz:**
- **Latest Hopper GPUs:** Best integer performance
- **High memory bandwidth:** Critical for 128-bit arithmetic
- **Modern network:** Excellent for coordination
- **European location:** Accessible for European students