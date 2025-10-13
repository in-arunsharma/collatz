// worker_gpp_node.cpp - GPP Node Worker Implementation
// Pure CPU OpenMP parallelization for General Purpose Partition nodes
// This is your proven V1.4b-openmp code adapted for MPI

#include <iostream>
#include <cstdio>
#include <cstdint>
#include <omp.h>
#include <chrono>
#include <string>
#include <vector>

// Include node configuration
#include "node_config.hpp"

// V1.4b-openmp core integration
#include "collatz_core_v14b.cpp"
#include "marenostrum_config.hpp"

// Import V1.4b-openmp core functionality
// NOTE: You'll need to extract the core functions from V1.4b-openmp.cpp
// and put them in a header or link against them

using namespace mn5;

// Forward declarations (these would come from your V1.4b-openmp core)
// For now, we'll create simplified versions

typedef unsigned __int128 uint128_t;

// Simplified Collatz computation result
struct CollatzResult {
    uint64_t steps;
    uint128_t peak;
    bool overflow;
    bool fuse_hit;
};

// Thread-local statistics (from V1.4b-openmp)
struct alignas(64) ThreadStats {
    uint64_t tested = 0;
    uint64_t total_steps = 0;
    uint64_t max_steps_seen = 0;
    uint128_t hardest_n = 0;
    uint128_t max_peak = 0;
    uint64_t overflow_count = 0;
    uint64_t fuse_count = 0;
};

// Deterministic seed generation (mod-6 filter)
inline uint128_t seed_from_index(uint128_t base_start, uint64_t idx) {
    // From V1.4b-openmp: map index to seed with mod-6 filter
    uint64_t pair = idx / 2;
    uint64_t odd = idx % 2;
    return base_start + (pair * 6) + (odd ? 4 : 0);
}

// Count total candidates after mod-6 filter
inline uint64_t count_candidates(uint128_t start, uint128_t end) {
    // Simplified: approximately 1/3 of range after mod-6 filter
    uint64_t range = (uint64_t)(end - start);
    return (range + 2) / 3;
}

// Mock Collatz computation (replace with actual V1.4b-openmp function)
inline CollatzResult compute_collatz_readonly(
    uint128_t n,
    const uint32_t* memo,
    uint64_t memo_size,
    uint64_t max_steps
) {
    CollatzResult result;
    result.steps = 0;
    result.peak = n;
    result.overflow = false;
    result.fuse_hit = false;
    
    // TODO: Import actual Collatz computation from V1.4b-openmp.cpp
    // For now, this is a placeholder
    
    return result;
}

// GPP Node Worker - processes assigned seed range with OpenMP
class GPPNodeWorker {
private:
    NodeConfig config_;
    uint128_t start_;
    uint128_t end_;
    std::string run_tag_;
    
    // Global statistics
    uint64_t global_tested_ = 0;
    uint64_t global_max_steps_ = 0;
    uint128_t global_max_steps_seed_ = 0;
    uint128_t global_max_peak_ = 0;
    uint64_t global_overflow_ = 0;
    uint64_t global_fuse_ = 0;
    
    // Memo table (shared read-only)
    std::vector<uint32_t> memo_;
    
public:
    GPPNodeWorker(const NodeConfig& config, uint128_t start, uint128_t end, 
                  const std::string& run_tag)
        : config_(config), start_(start), end_(end), run_tag_(run_tag) {
        
        // Initialize memo table
        memo_.resize(CollatzParameters::MEMO_SIZE, 0);
        // TODO: Fill memo table (from V1.4b-openmp precomputation)
    }
    
    // Process the assigned range
    void process() {
        fprintf(stderr, "[Rank %d] GPP Worker starting...\n", config_.mpi_rank);
        fprintf(stderr, "[Rank %d] Range size: %lu seeds\n", config_.mpi_rank,
                (uint64_t)(end_ - start_));
        fprintf(stderr, "[Rank %d] OpenMP threads: %d\n", config_.mpi_rank,
                config_.openmp_threads);
        
        auto t_start = std::chrono::high_resolution_clock::now();
        
        // Calculate total candidates after mod-6 filter
        uint64_t total_candidates = count_candidates(start_, end_);
        
        // OpenMP parallel processing (V1.4b-openmp style)
        #pragma omp parallel num_threads(config_.openmp_threads)
        {
            ThreadStats local;
            
            #pragma omp for schedule(static, 10000) nowait
            for (uint64_t idx = 0; idx < total_candidates; ++idx) {
                uint128_t n = seed_from_index(start_, idx);
                
                // Compute Collatz sequence
                CollatzResult res = compute_collatz_readonly(
                    n, memo_.data(), memo_.size(), 
                    CollatzParameters::HOT_FUSE
                );
                
                local.tested++;
                local.total_steps += res.steps;
                
                if (res.steps > local.max_steps_seen) {
                    local.max_steps_seen = res.steps;
                    local.hardest_n = n;
                }
                
                if (res.peak > local.max_peak) {
                    local.max_peak = res.peak;
                }
                
                if (res.overflow) local.overflow_count++;
                if (res.fuse_hit) local.fuse_count++;
                
                // Progress reporting
                if (local.tested % CollatzParameters::PROGRESS_INTERVAL == 0) {
                    #pragma omp critical
                    {
                        fprintf(stderr, "[Rank %d, Thread %d] Progress: %lu tested\n",
                                config_.mpi_rank, omp_get_thread_num(), local.tested);
                    }
                }
            }
            
            // Merge thread-local statistics
            #pragma omp critical
            {
                global_tested_ += local.tested;
                
                if (local.max_steps_seen > global_max_steps_) {
                    global_max_steps_ = local.max_steps_seen;
                    global_max_steps_seed_ = local.hardest_n;
                }
                
                if (local.max_peak > global_max_peak_) {
                    global_max_peak_ = local.max_peak;
                }
                
                global_overflow_ += local.overflow_count;
                global_fuse_ += local.fuse_count;
            }
        }
        
        auto t_end = std::chrono::high_resolution_clock::now();
        auto duration_ms = std::chrono::duration_cast<std::chrono::milliseconds>(
            t_end - t_start).count();
        
        // Print results
        print_results(duration_ms);
    }
    
    // Get results for MPI communication
    void get_results(uint64_t& tested, uint64_t& time_ms, uint64_t& max_steps,
                     uint128_t& max_steps_seed, uint128_t& max_peak,
                     uint64_t& overflow, uint64_t& fuse) const {
        tested = global_tested_;
        max_steps = global_max_steps_;
        max_steps_seed = global_max_steps_seed_;
        max_peak = global_max_peak_;
        overflow = global_overflow_;
        fuse = global_fuse_;
        // time_ms needs to be passed separately
    }
    
private:
    void print_results(uint64_t time_ms) {
        double throughput = global_tested_ / (time_ms / 1000.0);
        
        fprintf(stderr, "\n[Rank %d] === GPP Node Results ===\n", config_.mpi_rank);
        fprintf(stderr, "  Tested:       %lu numbers\n", global_tested_);
        fprintf(stderr, "  Time:         %lu ms\n", time_ms);
        fprintf(stderr, "  Throughput:   %.2f M nums/sec\n", throughput / 1e6);
        fprintf(stderr, "  Max steps:    %lu\n", global_max_steps_);
        fprintf(stderr, "  Overflow:     %lu\n", global_overflow_);
        fprintf(stderr, "  Fuse hits:    %lu\n", global_fuse_);
        fprintf(stderr, "\n");
    }
};

// Export for MPI main
extern "C" {
    void* create_gpp_worker(int mpi_rank, int mpi_size, 
                            uint64_t start_lo, uint64_t start_hi,
                            uint64_t end_lo, uint64_t end_hi,
                            const char* run_tag) {
        NodeConfig config = detect_node_configuration(mpi_rank, mpi_size);
        
        uint128_t start = ((uint128_t)start_hi << 64) | start_lo;
        uint128_t end = ((uint128_t)end_hi << 64) | end_lo;
        
        return new GPPNodeWorker(config, start, end, run_tag);
    }
    
    void gpp_worker_process(void* worker_ptr) {
        GPPNodeWorker* worker = static_cast<GPPNodeWorker*>(worker_ptr);
        worker->process();
    }
    
    void destroy_gpp_worker(void* worker_ptr) {
        delete static_cast<GPPNodeWorker*>(worker_ptr);
    }
}

#endif // WORKER_GPP_NODE_CPP
