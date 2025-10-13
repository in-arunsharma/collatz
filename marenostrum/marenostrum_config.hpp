// marenostrum_config.hpp - Configuration Constants for MareNostrum 5 Deployment
// Centralized configuration for hackathon resource allocation and parameters

#ifndef MARENOSTRUM_CONFIG_HPP
#define MARENOSTRUM_CONFIG_HPP

#include <cstdint>

namespace mn5 {

// ============================================================================
// HACKATHON RESOURCE ALLOCATION
// ============================================================================

// Maximum resources allocated to your team
// Adjust these based on actual allocation from organizers
struct ResourceLimits {
    // GPP Nodes (General Purpose - CPU only)
    static constexpr int MAX_GPP_NODES = 10;
    static constexpr int REQUESTED_GPP_NODES = 5;  // Start conservative
    
    // ACC Nodes (Accelerated - GPU + CPU)
    static constexpr int MAX_ACC_NODES = 5;
    static constexpr int REQUESTED_ACC_NODES = 2;  // Start with 2, scale up
    
    // Total resources
    static constexpr int MAX_TOTAL_NODES = MAX_GPP_NODES + MAX_ACC_NODES;
    static constexpr int REQUESTED_TOTAL_NODES = REQUESTED_GPP_NODES + REQUESTED_ACC_NODES;
    
    // Time limits (in hours)
    static constexpr int MAX_WALLTIME_HOURS = 6;
    static constexpr int INITIAL_WALLTIME_HOURS = 2;  // Start with 2h, extend if needed
};

// ============================================================================
// HARDWARE SPECIFICATIONS (MareNostrum 5 official specs)
// ============================================================================

struct GPPNodeSpec {
    static constexpr int CPUS_PER_NODE = 2;
    static constexpr int CORES_PER_CPU = 56;
    static constexpr int TOTAL_CORES = CPUS_PER_NODE * CORES_PER_CPU;  // 112
    static constexpr double CPU_FREQ_GHZ = 2.0;
    static constexpr int MEMORY_GB = 256;  // Standard nodes (some have 1024GB)
    static constexpr const char* CPU_MODEL = "Intel Sapphire Rapids 8480+";
};

struct ACCNodeSpec {
    static constexpr int CPUS_PER_NODE = 2;
    static constexpr int CORES_PER_CPU = 40;
    static constexpr int TOTAL_CORES = CPUS_PER_NODE * CORES_PER_CPU;  // 80
    static constexpr double CPU_FREQ_GHZ = 2.3;
    static constexpr int MEMORY_GB = 512;
    
    static constexpr int GPUS_PER_NODE = 4;
    static constexpr int GPU_MEMORY_GB = 64;
    static constexpr const char* GPU_MODEL = "NVIDIA Hopper H100";
    static constexpr const char* CPU_MODEL = "Intel Sapphire Rapids 8460Y+";
};

// ============================================================================
// COLLATZ COMPUTATION PARAMETERS
// ============================================================================

struct CollatzParameters {
    // Starting point: 2^71
    static constexpr uint64_t BASE_POWER = 71;
    
    // Memo table configuration
    static constexpr int MEMO_BITS = 20;  // 2^20 = 1M entries
    static constexpr uint64_t MEMO_SIZE = 1ULL << MEMO_BITS;
    
    // Safety fuses
    static constexpr uint64_t HOT_FUSE = 100000;      // Fast path limit
    static constexpr uint64_t COLD_FUSE = 1000000;    // Extended limit for cold queue
    
    // Batch processing (GPU)
    static constexpr uint64_t GPU_BATCH_SIZE = 1000000;  // 1M seeds per GPU batch
    
    // Mod-6 filter: only test n ≡ 1,5 (mod 6)
    static constexpr int MOD_FILTER = 6;
    static constexpr int MOD_CANDIDATES[2] = {1, 5};
    
    // Progress reporting interval
    static constexpr uint64_t PROGRESS_INTERVAL = 100000;
};

// ============================================================================
// PERFORMANCE ESTIMATES (for planning)
// ============================================================================

struct PerformanceEstimates {
    // CPU throughput (based on local testing)
    static constexpr double CPU_CORE_THROUGHPUT = 2.3e6;  // 2.3M nums/sec per core
    
    // GPU throughput (conservative estimate)
    static constexpr double GPU_THROUGHPUT = 50e6;  // 50M nums/sec per GPU
    
    // Expected cluster performance
    static constexpr double GPP_NODE_THROUGHPUT = GPPNodeSpec::TOTAL_CORES * CPU_CORE_THROUGHPUT;
    static constexpr double ACC_NODE_THROUGHPUT = 
        (ACCNodeSpec::TOTAL_CORES * CPU_CORE_THROUGHPUT * 0.2) +  // 20% CPU for cold queue
        (ACCNodeSpec::GPUS_PER_NODE * GPU_THROUGHPUT);              // GPU hot path
    
    // Total expected throughput
    static constexpr double TOTAL_THROUGHPUT = 
        (ResourceLimits::REQUESTED_GPP_NODES * GPP_NODE_THROUGHPUT) +
        (ResourceLimits::REQUESTED_ACC_NODES * ACC_NODE_THROUGHPUT);
    
    // Expected: ~1.3B (GPP) + 0.4B (ACC) = ~1.7 billion nums/sec at start
};

// ============================================================================
// FILE NAMING CONVENTIONS
// ============================================================================

struct FileNaming {
    // Output directory on MareNostrum
    static constexpr const char* OUTPUT_DIR = "/gpfs/scratch/bsc/output/collatz";
    
    // Log file prefixes
    static constexpr const char* OVERFLOW_LOG_PREFIX = "overflow_seeds";
    static constexpr const char* FUSE_LOG_PREFIX = "fuse_seeds";
    static constexpr const char* CYCLE_LOG_PREFIX = "cycle_seeds";
    static constexpr const char* BIG_OVERFLOW_LOG_PREFIX = "256bit_overflow";
    
    // Metadata
    static constexpr const char* METADATA_PREFIX = "run_metadata";
    static constexpr const char* METADATA_EXTENSION = ".json";
    
    // Per-rank output pattern: <prefix>_rank<rank>_<tag>.txt
};

// ============================================================================
// MPI CONFIGURATION
// ============================================================================

struct MPIConfig {
    // Rank 0 is always the master/coordinator
    static constexpr int MASTER_RANK = 0;
    
    // MPI tags for different message types
    static constexpr int TAG_WORK_ASSIGNMENT = 100;
    static constexpr int TAG_RESULTS = 101;
    static constexpr int TAG_PROGRESS = 102;
    static constexpr int TAG_SHUTDOWN = 103;
};

// ============================================================================
// SLURM PARTITION NAMES (MareNostrum 5 specific)
// ============================================================================

struct SlurmPartitions {
    static constexpr const char* GPP = "gpp";  // General Purpose Partition
    static constexpr const char* ACC = "acc";  // Accelerated Partition
};

// ============================================================================
// HELPER FUNCTIONS
// ============================================================================

// Calculate total expected cores
inline constexpr int get_total_cores() {
    return (ResourceLimits::REQUESTED_GPP_NODES * GPPNodeSpec::TOTAL_CORES) +
           (ResourceLimits::REQUESTED_ACC_NODES * ACCNodeSpec::TOTAL_CORES);
}

// Calculate total GPUs
inline constexpr int get_total_gpus() {
    return ResourceLimits::REQUESTED_ACC_NODES * ACCNodeSpec::GPUS_PER_NODE;
}

// Calculate expected throughput
inline constexpr double get_expected_throughput() {
    return PerformanceEstimates::TOTAL_THROUGHPUT;
}

// Print configuration summary
inline void print_config_summary() {
    fprintf(stderr, "\n=== MareNostrum 5 Configuration ===\n");
    fprintf(stderr, "Requested Resources:\n");
    fprintf(stderr, "  GPP Nodes: %d (max %d)\n", 
            ResourceLimits::REQUESTED_GPP_NODES, ResourceLimits::MAX_GPP_NODES);
    fprintf(stderr, "  ACC Nodes: %d (max %d)\n", 
            ResourceLimits::REQUESTED_ACC_NODES, ResourceLimits::MAX_ACC_NODES);
    fprintf(stderr, "  Total Nodes: %d\n", ResourceLimits::REQUESTED_TOTAL_NODES);
    fprintf(stderr, "  Total CPU Cores: %d\n", get_total_cores());
    fprintf(stderr, "  Total GPUs: %d\n", get_total_gpus());
    fprintf(stderr, "  Walltime: %d hours\n", ResourceLimits::INITIAL_WALLTIME_HOURS);
    fprintf(stderr, "\n");
    fprintf(stderr, "Expected Performance:\n");
    fprintf(stderr, "  Total Throughput: %.2f billion nums/sec\n", 
            get_expected_throughput() / 1e9);
    fprintf(stderr, "  Time for 1 trillion: %.2f minutes\n",
            1e12 / get_expected_throughput() / 60.0);
    fprintf(stderr, "\n");
    fprintf(stderr, "Collatz Parameters:\n");
    fprintf(stderr, "  Base: 2^%lu\n", CollatzParameters::BASE_POWER);
    fprintf(stderr, "  Memo size: 2^%d = %lu entries\n", 
            CollatzParameters::MEMO_BITS, CollatzParameters::MEMO_SIZE);
    fprintf(stderr, "  Hot fuse: %lu steps\n", CollatzParameters::HOT_FUSE);
    fprintf(stderr, "  Cold fuse: %lu steps\n", CollatzParameters::COLD_FUSE);
    fprintf(stderr, "\n");
}

} // namespace mn5

#endif // MARENOSTRUM_CONFIG_HPP
