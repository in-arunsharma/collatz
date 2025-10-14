// MNv2 Configuration - MareNostrum 5 Hackathon
// 
// HACKATHON TEAM: Adjust these values based on your node allocation!
// Shared pool across all teams - coordinate resource usage
//
// Quick Start:
// 1. Set GPP_NODES and ACC_NODES based on what you get allocated
// 2. Run ./build_gpp.sh for Phase 1 (GPP-only testing)
// 3. Run ./build_acc.sh for Phase 2 (GPU acceleration)
// 4. Use sbatch slurm_gpp.slurm or slurm_acc.slurm

#pragma once

#include <cstdint>

namespace MNv2 {

// ============================================================================
// CLUSTER CONFIGURATION - MODIFY THESE FOR YOUR ALLOCATION
// ============================================================================

// Phase 1: GPP Nodes (CPU-only, simpler to validate first)
constexpr int GPP_NODES = 1;              // Single node optimization first
constexpr int GPP_CORES_PER_NODE = 112;   // MareNostrum 5 GPP: 2× 56-core Sapphire Rapids
constexpr int GPP_THREADS_PER_CORE = 1;  // No hyperthreading for computational work

// Phase 2: ACC Nodes (GPU+CPU hybrid, more complex)
constexpr int ACC_NODES = 0;              // Set to 1-3 once GPP is working
constexpr int ACC_CORES_PER_NODE = 80;    // MareNostrum 5 ACC: 2× 40-core Sapphire Rapids
constexpr int ACC_GPUS_PER_NODE = 4;      // 4× NVIDIA H100 per ACC node
constexpr int ACC_THREADS_PER_CORE = 1;  // No hyperthreading

// Total resources (auto-calculated)
constexpr int TOTAL_GPP_CORES = GPP_NODES * GPP_CORES_PER_NODE;
constexpr int TOTAL_ACC_CORES = ACC_NODES * ACC_CORES_PER_NODE;
constexpr int TOTAL_ACC_GPUS = ACC_NODES * ACC_GPUS_PER_NODE;
constexpr int TOTAL_NODES = GPP_NODES + ACC_NODES;

// ============================================================================
// WORKLOAD CONFIGURATION
// ============================================================================

// Collatz search range
constexpr int BASE_POWER = 71;                    // Search from 2^71
constexpr uint64_t DEFAULT_RANGE_SIZE = 1000000000ULL;  // 1 billion seeds (adjust per run)

// Memo table (thread-safe read-only after precompute)
// CACHE OPTIMIZATION: Sweet spot between L2 and L3 cache
constexpr uint32_t MEMO_TABLE_BITS = 19;         // 2^19 = 512K entries = 2MB (L2 boundary)
constexpr uint64_t MEMO_TABLE_SIZE = (1ULL << MEMO_TABLE_BITS);

// Safety fuses
constexpr uint64_t HOT_PATH_FUSE = 100000;       // Hot path iteration limit
constexpr uint64_t COLD_QUEUE_FUSE = 1000000;    // Cold queue extended limit (10×)

// ============================================================================
// MPI CONFIGURATION
// ============================================================================

// MPI process layout:
// - Rank 0: Master (coordinator, no computation)
// - Rank 1..N: Workers (GPP or ACC nodes)
constexpr int MASTER_RANK = 0;

// How to partition work across workers
enum class WorkDistribution {
    EQUAL_RANGE,     // Divide range equally by count (simple, load-balanced)
    BLOCK_CYCLIC,    // Interleave blocks (better for irregular workloads)
    DYNAMIC          // Master dispatches chunks on-demand (future: load balancing)
};
constexpr WorkDistribution WORK_MODE = WorkDistribution::EQUAL_RANGE;

// ============================================================================
// PERFORMANCE TUNING
// ============================================================================

// OpenMP settings (exported in SLURM script)
// export OMP_NUM_THREADS=$CORES_PER_NODE
// export OMP_PROC_BIND=close
// export OMP_PLACES=cores

// Progress reporting
constexpr bool ENABLE_PROGRESS = false;          // Set true for debugging, false for benchmarks
constexpr uint64_t PROGRESS_INTERVAL = 16384;   // Report every 2^14 seeds (if enabled)

// Output control
constexpr bool VERBOSE_WORKERS = true;           // Workers print progress
constexpr bool VERBOSE_MASTER = true;            // Master prints aggregation

// ============================================================================
// HARDWARE SPECS (MareNostrum 5) - DO NOT MODIFY
// ============================================================================

// GPP Node Specs
struct GPPNodeSpec {
    static constexpr const char* name = "General Purpose Partition";
    static constexpr int cpus = 2;                          // 2× Intel Sapphire Rapids 8480+
    static constexpr int cores_per_cpu = 56;
    static constexpr int total_cores = 112;
    static constexpr double cpu_freq_ghz = 2.0;            // Base frequency
    static constexpr int memory_gb = 256;                  // Some nodes have 1TB
    static constexpr int l3_cache_mb = 105;                // 105MB L3 per CPU
};

// ACC Node Specs  
struct ACCNodeSpec {
    static constexpr const char* name = "Accelerated Partition";
    static constexpr int cpus = 2;                          // 2× Intel Sapphire Rapids 8460Y+
    static constexpr int cores_per_cpu = 40;
    static constexpr int total_cores = 80;
    static constexpr double cpu_freq_ghz = 2.3;
    static constexpr int memory_gb = 512;
    static constexpr int gpus = 4;                          // 4× NVIDIA H100
    static constexpr int gpu_memory_gb = 64;               // 64GB HBM2 per GPU
    static constexpr double gpu_tflops_fp64 = 34;          // H100 FP64 peak
};

// ============================================================================
// VERSION INFO
// ============================================================================

constexpr const char* VERSION = "MNv2-Phase1";
constexpr const char* BUILD_DATE = __DATE__;
constexpr const char* BUILD_TIME = __TIME__;

// Git commit (injected at compile time via -DGIT_COMMIT="hash")
#ifndef GIT_COMMIT
#define GIT_COMMIT "unknown"
#endif

// ============================================================================
// HELPER FUNCTIONS
// ============================================================================

inline void print_configuration() {
    printf("\n");
    printf("╔═══════════════════════════════════════════════════════════════╗\n");
    printf("║          MareNostrum 5 Collatz - Configuration               ║\n");
    printf("╚═══════════════════════════════════════════════════════════════╝\n");
    printf("\n");
    printf("Version:           %s\n", VERSION);
    printf("Build:             %s %s\n", BUILD_DATE, BUILD_TIME);
    printf("Git commit:        %s\n", GIT_COMMIT);
    printf("\n");
    printf("--- Resource Allocation ---\n");
    printf("GPP Nodes:         %d × %d cores = %d total cores\n", 
           GPP_NODES, GPP_CORES_PER_NODE, TOTAL_GPP_CORES);
    if (ACC_NODES > 0) {
        printf("ACC Nodes:         %d × %d cores = %d total cores\n",
               ACC_NODES, ACC_CORES_PER_NODE, TOTAL_ACC_CORES);
        printf("ACC GPUs:          %d × %d GPUs = %d total GPUs\n",
               ACC_NODES, ACC_GPUS_PER_NODE, TOTAL_ACC_GPUS);
    }
    printf("Total Nodes:       %d\n", TOTAL_NODES);
    printf("\n");
    printf("--- Workload ---\n");
    printf("Base:              2^%d\n", BASE_POWER);
    printf("Memo table:        2^%u = %llu entries (%.1f MB)\n",
           MEMO_TABLE_BITS, MEMO_TABLE_SIZE, 
           (MEMO_TABLE_SIZE * 4.0) / (1024 * 1024));
    printf("Hot path fuse:     %llu iterations\n", HOT_PATH_FUSE);
    printf("Cold queue fuse:   %llu iterations\n", COLD_QUEUE_FUSE);
    printf("\n");
    printf("--- MPI Configuration ---\n");
    printf("Master rank:       %d (coordinator only)\n", MASTER_RANK);
    printf("Worker ranks:      1..%d\n", TOTAL_NODES);
    printf("Work distribution: ");
    switch (WORK_MODE) {
        case WorkDistribution::EQUAL_RANGE: printf("Equal range\n"); break;
        case WorkDistribution::BLOCK_CYCLIC: printf("Block-cyclic\n"); break;
        case WorkDistribution::DYNAMIC: printf("Dynamic\n"); break;
    }
    printf("\n");
}

inline int get_openmp_threads_for_node_type(bool is_gpp) {
    return is_gpp ? GPP_CORES_PER_NODE : ACC_CORES_PER_NODE;
}

inline const char* get_node_type_name(bool is_gpp) {
    return is_gpp ? "GPP" : "ACC";
}

} // namespace MNv2
