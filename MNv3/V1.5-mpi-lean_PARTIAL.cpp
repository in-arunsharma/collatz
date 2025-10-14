// V1.5-mpi-lean: Minimal MPI wrapper around proven V1.5-openmp kernel
//
// CHANGES FROM V1.5-openmp:
// 1. Added MPI_Init/Finalize
// 2. Work distribution: Simple offset splitting per rank
// 3. Seeds computed on-the-fly (NO big vector) - friend's key optimization
// 4. Static scheduling (uniform work) - friend's recommendation
// 5. Optional --cold off flag for pure throughput
//
// UNCHANGED FROM V1.5:
// - All computation kernels (compute_collatz_readonly, Brent, etc.)
// - Thread-local statistics
// - Logging infrastructure
// - Memo table precomputation
//
// TARGET: 2 nodes × 137M/sec = 274M/sec

#include <mpi.h>
#include <iostream>
#include <vector>
#include <cstdint>
#include <cstring>
#include <cstdio>
#include <ctime>
#include <time.h>
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <mutex>
#include <deque>
#include <omp.h>

typedef __uint128_t uint128_t;

// ----- All V1.5 structures and functions (UNCHANGED) -----
// (Keeping your exact V1.5 code - only showing the additions)

// [Insert all your V1.5 code here - uint256_t, print_u128, etc.]
// For brevity, I'll just show the new/changed parts:

// ... [ALL YOUR V1.5 CODE FROM LINE 13 TO LINE 643] ...

// NEW: Helper to count candidates without building vector
static inline uint64_t count_candidates_u128(uint128_t start, uint128_t end) {
    if (end <= start) return 0;
    // Mod-6 filtering: 2 out of every 6 numbers (n odd and n%3!=0)
    return (uint64_t)(((end - start) * 2) / 3);
}

// REPLACED: MPI main()
int main(int argc, char** argv) {
    // MPI initialization
    MPI_Init(&argc, &argv);
    int rank = 0, size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    if (argc < 3) {
        if (rank == 0) {
            fprintf(stderr, "Usage: %s <start_offset> <count> [--threads N] [--small-limit BITS] [--tag NAME] [--cold off]\n", argv[0]);
        }
        MPI_Finalize();
        return 1;
    }
    
    // Parse arguments
    uint64_t global_start_offset = std::stoull(argv[1]);
    uint64_t global_count = std::stoull(argv[2]);
    
    const char* run_tag = "v15mpi";
    int num_threads = 0;
    bool cold_enabled = true;
    
    for (int i = 3; i < argc; i++) {
        if (!strcmp(argv[i], "--threads") && i+1 < argc) {
            num_threads = std::stoi(argv[++i]);
        } else if (!strcmp(argv[i], "--small-limit") && i+1 < argc) {
            SMALL_LIMIT_BITS = std::stoul(argv[++i]);
        } else if (!strcmp(argv[i], "--tag") && i+1 < argc) {
            run_tag = argv[++i];
        } else if (!strcmp(argv[i], "--cold") && i+1 < argc) {
            cold_enabled = strcmp(argv[++i], "off") != 0;
        }
    }
    
    // Work distribution per rank
    uint64_t block = global_count / size;
    uint64_t my_offset = global_start_offset + (uint64_t)rank * block;
    uint64_t my_count = (rank == size-1) ? (global_count - (uint64_t)rank * block) : block;
    
    // OpenMP setup
    if (num_threads > 0) omp_set_num_threads(num_threads);
    num_threads = omp_get_max_threads();
    
    if (rank == 0) {
        fprintf(stderr, "[CONFIG] ranks=%d  threads/rank=%d  total-threads=%d\n", 
                size, num_threads, size*num_threads);
        fprintf(stderr, "[HINT] Use: OMP_PLACES=cores  OMP_PROC_BIND=spread  and  --ntasks-per-node=1 --cpus-per-task=<all>\n");
    }
    
    // Logs (rank-tagged)
    char rank_tag[128];
    snprintf(rank_tag, sizeof(rank_tag), "%s_r%d", run_tag, rank);
    open_logs_once(rank_tag);
    
    // Memo table per rank
    small_limit = (1ULL << SMALL_LIMIT_BITS);
    if (rank == 0) {
        fprintf(stderr, "[CONFIG] Memo table: 2^%u = %llu entries (%.2f MB)\n",
                SMALL_LIMIT_BITS, (unsigned long long)small_limit,
                (small_limit*4.0)/(1024*1024));
    }
    memo.assign(small_limit, UNKNOWN);
    precompute_small_table();
    if (rank == 0) validate_memo_table();
    
    // Prefault
    volatile uint64_t pf = 0;
    for (uint64_t i = 0; i < small_limit; i += 1024) pf += memo[i];
    
    // Build my numeric range
    uint128_t base = (uint128_t)1 << 71;
    uint128_t start = base + my_offset;
    uint128_t end = start + my_count;
    
    // Align start to first candidate
    uint64_t delta0 = 4;
    align_start_and_delta(start, delta0);
    
    // Count candidates (NO VECTOR)
    const uint64_t total_candidates = count_candidates_u128(start, end);
    const uint32_t* __restrict memo_ptr = memo.data();
    
    // Stats
    ThreadStats stats{};
    struct timespec t0, t1;
    clock_gettime(CLOCK_MONOTONIC, &t0);
    
    // ===== KEY OPTIMIZATION: On-the-fly seed generation + static scheduling =====
    #pragma omp parallel
    {
        ThreadStats local{};
        std::deque<ColdQueueEntry> q1, q2;
        
        // STATIC scheduling (friend's advice)
        #pragma omp for schedule(static, 8192)
        for (uint64_t idx = 0; idx < total_candidates; ++idx) {
            // Compute seed from index (NO vector lookup)
            // n = start + floor(idx/2)*6 + (idx&1 ? 4 : 0)
            uint128_t n = start + ((uint128_t)(idx >> 1) * 6u) + ((idx & 1u) ? 4u : 0u);
            
            // V1.5 hot path (UNCHANGED)
            CollatzResult res = compute_collatz_readonly(n, memo_ptr, small_limit, SAFETY_FUSE);
            
            local.tested++;
            local.total_steps += res.steps;
            
            if (res.steps > local.max_steps_seen) {
                local.max_steps_seen = res.steps;
                local.hardest_n = n;
            }
            
            if (res.peak > local.max_peak) {
                local.max_peak = res.peak;
            }
            
            // Cold queue handling (optional)
            if (cold_enabled) {
                if (res.overflow) {
                    local.overflow_count++;
                    local.cold_q2_triggers++;
                    q2.emplace_back(n, res.steps, res.peak);
                } else if (res.steps >= SAFETY_FUSE) {
                    local.fuse_count++;
                    local.cold_q1_triggers++;
                    q1.emplace_back(n, res.steps, res.peak);
                }
                
                if (local.tested % COLD_BATCH_SIZE == 0 ||
                    q1.size() >= MAX_COLD_QUEUE_SIZE ||
                    q2.size() >= MAX_COLD_QUEUE_SIZE) {
                    process_cold_queue_fuse(q1, local, memo_ptr);
                    process_cold_queue_overflow(q2, local);
                }
            }
        }
        
        if (cold_enabled) {
            process_cold_queue_fuse(q1, local, memo_ptr);
            process_cold_queue_overflow(q2, local);
        }
        
        // Aggregate thread stats
        #pragma omp critical
        {
            stats.tested += local.tested;
            stats.total_steps += local.total_steps;
            
            if (local.max_steps_seen > stats.max_steps_seen) {
                stats.max_steps_seen = local.max_steps_seen;
                stats.hardest_n = local.hardest_n;
            }
            
            if (local.max_peak > stats.max_peak) {
                stats.max_peak = local.max_peak;
            }
            
            stats.overflow_count += local.overflow_count;
            stats.fuse_count += local.fuse_count;
            stats.cold_q1_triggers += local.cold_q1_triggers;
            stats.cold_q1_cycles_found += local.cold_q1_cycles_found;
            stats.cold_q1_verified_ok += local.cold_q1_verified_ok;
            stats.cold_q2_triggers += local.cold_q2_triggers;
            stats.cold_q2_cycles_found += local.cold_q2_cycles_found;
            stats.cold_q2_256bit_overflow += local.cold_q2_256bit_overflow;
            stats.cold_q2_verified_ok += local.cold_q2_verified_ok;
            stats.cold_processing_time_ms += local.cold_processing_time_ms;
        }
    }
    
    clock_gettime(CLOCK_MONOTONIC, &t1);
    double my_time = (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec) * 1e-9;
    
    // MPI reductions
    uint64_t sum_tested = 0, sum_steps = 0, max_steps_all = 0;
    double max_time = 0.0;
    
    MPI_Reduce(&stats.tested, &sum_tested, 1, MPI_UINT64_T, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&stats.total_steps, &sum_steps, 1, MPI_UINT64_T, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&stats.max_steps_seen, &max_steps_all, 1, MPI_UINT64_T, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&my_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    
    // Output (rank 0 only)
    if (rank == 0) {
        double throughput = sum_tested / max_time;
        double avg_steps = (double)sum_steps / (double)sum_tested;
        
        std::cout << "\n=== V1.5-mpi (lean) Results ===\n";
        std::cout << "MPI ranks:    " << size << "\n";
        std::cout << "Threads/rank: " << num_threads << "\n";
        std::cout << "Tested:       " << sum_tested << "\n";
        std::cout << "Wall time:    " << (uint64_t)(max_time*1000) << " ms\n";
        std::cout << "Throughput:   " << (uint64_t)throughput << " nums/sec\n";
        std::cout << "Per-rank:     " << (uint64_t)(throughput/size) << " nums/sec\n";
        std::cout << "Avg steps:    " << avg_steps << "\n";
        std::cout << "Max steps:    " << max_steps_all << "\n";
        std::cout << "Cold queues:  " << (cold_enabled ? "ON" : "OFF (hot run)") << "\n";
        std::cout << "\n=== SUCCESS ===\n";
    }
    
    // Close logs
    if (g_overflow_log.is_open()) g_overflow_log.close();
    if (g_fuse_log.is_open()) g_fuse_log.close();
    if (g_cycle_log.is_open()) g_cycle_log.close();
    if (g_256bit_overflow_log.is_open()) g_256bit_overflow_log.close();
    
    MPI_Finalize();
    return 0;
}
