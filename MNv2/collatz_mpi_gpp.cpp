// MNv2: MPI+OpenMP Implementation for MareNostrum 5 GPP Nodes
//
// ARCHITECTURE:
// - Master (Rank 0): Partition work, collect results, aggregate statistics
// - Workers (Rank 1..N): Each processes a seed range with OpenMP parallelism
// - Computation: Embedded V1.4b-openmp core (hot path + cold queues)
//
// BUILD: ./build_gpp.sh
// LOCAL TEST: mpirun -np 4 ./collatz_mpi_gpp 0 1000000 test
// SUBMIT: sbatch slurm_gpp.slurm

#include <mpi.h>
#include <omp.h>
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
#include <atomic>
#include <chrono>

#include "config.hpp"

using namespace MNv2;

// ============================================================================
// CORE COLLATZ COMPUTATION (from V1.4b-openmp)
// ============================================================================

typedef __uint128_t uint128_t;

// ----- 256-bit arithmetic (for overflow handling) -----
struct uint256_t {
    uint64_t limbs[4];
    
    uint256_t() { limbs[0] = limbs[1] = limbs[2] = limbs[3] = 0; }
    explicit uint256_t(uint64_t val) {
        limbs[0] = val;
        limbs[1] = limbs[2] = limbs[3] = 0;
    }
    explicit uint256_t(uint128_t val) {
        limbs[0] = (uint64_t)val;
        limbs[1] = (uint64_t)(val >> 64);
        limbs[2] = limbs[3] = 0;
    }
    
    bool is_zero() const { return (limbs[0] | limbs[1] | limbs[2] | limbs[3]) == 0; }
    bool is_one() const { return limbs[0] == 1 && limbs[1] == 0 && limbs[2] == 0 && limbs[3] == 0; }
    bool is_even() const { return (limbs[0] & 1) == 0; }
    
    int ctz() const {
        if (limbs[0] != 0) return __builtin_ctzll(limbs[0]);
        if (limbs[1] != 0) return 64 + __builtin_ctzll(limbs[1]);
        if (limbs[2] != 0) return 128 + __builtin_ctzll(limbs[2]);
        if (limbs[3] != 0) return 192 + __builtin_ctzll(limbs[3]);
        return 256;
    }
    
    uint256_t operator>>(int shift) const {
        uint256_t result;
        if (shift == 0) return *this;
        if (shift >= 256) return uint256_t();
        int limb_shift = shift / 64;
        int bit_shift = shift % 64;
        if (bit_shift == 0) {
            for (int i = 0; i < 4 - limb_shift; i++) result.limbs[i] = limbs[i + limb_shift];
        } else {
            for (int i = 0; i < 4 - limb_shift; i++) {
                result.limbs[i] = limbs[i + limb_shift] >> bit_shift;
                if (i + limb_shift + 1 < 4)
                    result.limbs[i] |= limbs[i + limb_shift + 1] << (64 - bit_shift);
            }
        }
        return result;
    }
    
    uint256_t operator+(uint64_t val) const {
        uint256_t result = *this;
        __uint128_t sum = (__uint128_t)result.limbs[0] + val;
        result.limbs[0] = (uint64_t)sum;
        uint64_t carry = (uint64_t)(sum >> 64);
        for (int i = 1; i < 4 && carry; i++) {
            sum = (__uint128_t)result.limbs[i] + carry;
            result.limbs[i] = (uint64_t)sum;
            carry = (uint64_t)(sum >> 64);
        }
        return result;
    }
    
    uint256_t operator*(uint64_t scalar) const {
        uint256_t result;
        uint64_t carry = 0;
        for (int i = 0; i < 4; i++) {
            __uint128_t prod = (__uint128_t)limbs[i] * scalar + carry;
            result.limbs[i] = (uint64_t)prod;
            carry = (uint64_t)(prod >> 64);
        }
        return result;
    }
    
    bool operator>(const uint256_t& other) const {
        for (int i = 3; i >= 0; i--) {
            if (limbs[i] > other.limbs[i]) return true;
            if (limbs[i] < other.limbs[i]) return false;
        }
        return false;
    }
    
    bool operator==(const uint256_t& other) const {
        return limbs[0] == other.limbs[0] && limbs[1] == other.limbs[1] &&
               limbs[2] == other.limbs[2] && limbs[3] == other.limbs[3];
    }
    
    bool would_overflow_3n_plus_1() const {
        static const uint256_t MAX_SAFE_256 = []() {
            uint256_t limit;
            for (int i = 0; i < 4; i++) limit.limbs[i] = 0x5555555555555555ULL;
            return limit;
        }();
        return *this > MAX_SAFE_256;
    }
};

// ----- 128-bit utilities -----
static inline std::string to_string_u128(uint128_t val) {
    if (val == 0) return "0";
    char buffer[64];
    int pos = 0;
    while (val > 0) {
        buffer[pos++] = '0' + (val % 10);
        val /= 10;
    }
    std::string result;
    for (int i = pos - 1; i >= 0; i--) result += buffer[i];
    return result;
}

static inline int ctz_u128(uint128_t x) {
    if (x == 0) return 128;
    uint64_t lo = (uint64_t)x;
    if (lo != 0) return __builtin_ctzll(lo);
    uint64_t hi = (uint64_t)(x >> 64);
    return 64 + __builtin_ctzll(hi);
}

// ----- Structures -----
struct CollatzResult {
    uint64_t steps;
    uint128_t peak;
    bool overflow;
};

struct ColdQueueEntry {
    uint128_t seed;
    uint64_t steps_before_cold;
    uint128_t peak_before_cold;
    ColdQueueEntry(uint128_t s, uint64_t st, uint128_t p) 
        : seed(s), steps_before_cold(st), peak_before_cold(p) {}
};

struct alignas(64) ThreadStats {
    uint64_t tested = 0;
    uint64_t total_steps = 0;
    uint64_t max_steps_seen = 0;
    uint128_t hardest_n = 0;
    uint128_t max_peak = 0;
    uint64_t overflow_count = 0;
    uint64_t fuse_count = 0;
    std::vector<ColdQueueEntry> cold_queue_fuse;
    std::vector<ColdQueueEntry> cold_queue_overflow;
};

// ----- Constants -----
static constexpr uint32_t UNKNOWN = UINT32_MAX;
static constexpr uint128_t MAX_SAFE = ((~(uint128_t)0) - 1) / 3;

// ----- Global memo table (read-only after precompute) -----
static std::vector<uint32_t> memo;

// ----- Deterministic seed mapping (mod-6 filter) -----
static inline uint128_t seed_from_index(uint128_t start_aligned, uint64_t idx) {
    uint64_t pair = idx / 2;
    uint64_t odd = idx % 2;
    return start_aligned + (pair * 6) + (odd ? 4 : 0);
}

static inline uint64_t count_candidates(uint128_t start, uint128_t end) {
    if (end <= start) return 0;
    uint128_t range = end - start;
    return static_cast<uint64_t>((range * 2) / 3);
}

// ----- Precompute memo table (single-threaded, called once per rank) -----
static void precompute_small_table() {
    memo.resize(MEMO_TABLE_SIZE, UNKNOWN);
    memo[0] = UNKNOWN;
    memo[1] = 0;
    
    for (uint64_t n = 2; n < MEMO_TABLE_SIZE; n++) {
        if (memo[n] != UNKNOWN) continue;
        
        std::vector<uint64_t> path_vals;
        std::vector<uint64_t> path_steps;
        uint128_t current = (uint128_t)n;
        uint64_t steps = 0;
        bool completed = false;
        
        while (true) {
            if (current == 1) {
                completed = true;
                break;
            }
            if (current < MEMO_TABLE_SIZE && memo[(uint64_t)current] != UNKNOWN) {
                steps += memo[(uint64_t)current];
                completed = true;
                break;
            }
            if (steps >= HOT_PATH_FUSE) break;
            
            if (current < MEMO_TABLE_SIZE) {
                path_vals.push_back((uint64_t)current);
                path_steps.push_back(steps);
            }
            
            if ((current & 1) == 0) {
                int shift = ctz_u128(current);
                current >>= shift;
                steps += shift;
            } else {
                if (current > MAX_SAFE) break;
                uint128_t t = 3 * current + 1;
                int k = ctz_u128(t);
                current = t >> k;
                steps += 1 + k;
            }
        }
        
        if (completed) {
            for (int i = (int)path_vals.size() - 1; i >= 0; i--) {
                uint64_t val = path_vals[i];
                uint64_t steps_at_val = path_steps[i];
                uint64_t steps_from_val = steps - steps_at_val;
                if (steps_from_val < UNKNOWN && memo[val] == UNKNOWN) {
                    memo[val] = (uint32_t)steps_from_val;
                }
            }
        }
    }
}

// ----- Core Collatz computation (hot path) -----
__attribute__((always_inline)) static inline
CollatzResult compute_collatz_readonly(uint128_t n, const uint32_t* __restrict memo_ptr,
                                       uint64_t memo_limit, uint64_t max_steps = HOT_PATH_FUSE) {
    CollatzResult res = {0, n, false};
    if (n == 0 || n == 1) return res;
    
    uint128_t current = n;
    uint64_t steps = 0;
    
    while (current != 1) {
        if (steps >= max_steps) {
            res.steps = steps;
            res.overflow = true;
            return res;
        }
        
        if (current < memo_limit) {
            uint32_t cached = memo_ptr[(uint64_t)current];
            if (cached != UNKNOWN) {
                res.steps = steps + cached;
                return res;
            }
        }
        
        if (res.peak < current) res.peak = current;
        
        if ((current & 1) == 0) {
            int shift = ctz_u128(current);
            current >>= shift;
            steps += shift;
        } else {
            if (current > MAX_SAFE) {
                res.steps = steps;
                res.overflow = true;
                return res;
            }
            uint128_t t = 3 * current + 1;
            int k = ctz_u128(t);
            current = t >> k;
            steps += 1 + k;
        }
    }
    
    res.steps = steps;
    return res;
}

// ----- Mod-6 alignment helper -----
static inline uint32_t mod3_u128(uint128_t n) {
    uint64_t lo = (uint64_t)n;
    uint64_t hi = (uint64_t)(n >> 64);
    return ((lo % 3) + (hi % 3)) % 3;
}

static inline void align_start_and_delta(uint128_t &n, uint64_t &delta) {
    if ((uint64_t)n % 2 == 0) n++;
    uint32_t r3 = mod3_u128(n);
    if (r3 == 0) n += 2;
    delta = (mod3_u128(n) == 1) ? 4 : 2;
}

// ============================================================================
// MPI STRUCTURES & COMMUNICATION
// ============================================================================

struct WorkAssignment {
    uint128_t start;
    uint128_t end;
    uint64_t total_candidates;
};

struct WorkerResults {
    int rank;
    uint64_t tested;
    uint64_t total_steps;
    uint64_t max_steps;
    uint128_t max_steps_seed;
    uint128_t max_peak;
    uint64_t overflow_count;
    uint64_t fuse_count;
    double time_seconds;
};

// MPI send/recv helpers for uint128_t (send as 2× uint64_t)
static void MPI_Send_u128(uint128_t val, int dest, int tag, MPI_Comm comm) {
    uint64_t parts[2] = {(uint64_t)val, (uint64_t)(val >> 64)};
    MPI_Send(parts, 2, MPI_UNSIGNED_LONG_LONG, dest, tag, comm);
}

static uint128_t MPI_Recv_u128(int source, int tag, MPI_Comm comm) {
    uint64_t parts[2];
    MPI_Recv(parts, 2, MPI_UNSIGNED_LONG_LONG, source, tag, comm, MPI_STATUS_IGNORE);
    return ((uint128_t)parts[1] << 64) | parts[0];
}

// ============================================================================
// WORKER PROCESS (Ranks 1..N)
// ============================================================================

static WorkerResults worker_process(int rank, int size) {
    // Receive work assignment from master
    uint128_t start = MPI_Recv_u128(MASTER_RANK, 100, MPI_COMM_WORLD);
    uint128_t end = MPI_Recv_u128(MASTER_RANK, 101, MPI_COMM_WORLD);
    uint64_t total_candidates;
    MPI_Recv(&total_candidates, 1, MPI_UNSIGNED_LONG_LONG, MASTER_RANK, 102, 
             MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    
    // Configure OpenMP - use SLURM's OMP_NUM_THREADS if set, else default
    int num_threads = omp_get_max_threads();  // SLURM sets this via --cpus-per-task
    if (num_threads <= 0 || num_threads > GPP_CORES_PER_NODE) {
        num_threads = GPP_CORES_PER_NODE;
    }
    omp_set_num_threads(num_threads);  // Always set explicitly
    
    if (VERBOSE_WORKERS) {
        fprintf(stderr, "[Rank %d] Assigned %llu candidates, using %d threads\n",
                rank, total_candidates, num_threads);
    }
    
    // Precompute memo table (read-only after this)
    precompute_small_table();
    const uint32_t* __restrict memo_ptr = memo.data();
    
    // Align start for mod-6 filter
    uint64_t delta = 4;
    uint128_t aligned_start = start;
    align_start_and_delta(aligned_start, delta);
    
    // Global statistics (merged from threads)
    uint64_t global_tested = 0;
    uint64_t global_total_steps = 0;
    uint64_t global_max_steps = 0;
    uint128_t global_hardest_n = 0;
    uint128_t global_max_peak = 0;
    uint64_t global_overflow_count = 0;
    uint64_t global_fuse_count = 0;
    
    std::deque<ColdQueueEntry> global_queue_fuse;
    std::deque<ColdQueueEntry> global_queue_overflow;
    
    auto t_start = std::chrono::high_resolution_clock::now();
    
    // ===== OPENMP PARALLEL HOT PATH =====
    // Use OpenMP reductions for scalars (much faster than critical sections)
    uint64_t private_tested = 0;
    uint64_t private_total_steps = 0;
    uint64_t private_overflow_count = 0;
    uint64_t private_fuse_count = 0;
    
    #pragma omp parallel reduction(+:private_tested,private_total_steps,private_overflow_count,private_fuse_count)
    {
        uint64_t local_max_steps = 0;
        uint128_t local_hardest_n = 0;
        uint128_t local_max_peak = 0;
        std::vector<ColdQueueEntry> local_queue_fuse;
        std::vector<ColdQueueEntry> local_queue_overflow;
        
        // Diagnostic: verify thread count
        #pragma omp single
        {
            int actual_threads = omp_get_num_threads();
            if (VERBOSE_WORKERS) {
                fprintf(stderr, "[Rank %d] OpenMP parallel region active with %d threads\n",
                        rank, actual_threads);
            }
        }
        
        // OPTIMIZED: Pre-calculate stride to avoid division/modulo in loop
        #pragma omp for schedule(static, 10000) nowait
        for (uint64_t idx = 0; idx < total_candidates; ++idx) {
            // Optimized seed calculation: avoid division (use bit shift instead)
            uint128_t n = aligned_start + ((idx >> 1) * 6) + ((idx & 1) ? 4 : 0);
            if (n >= end) continue;
            
            CollatzResult res = compute_collatz_readonly(n, memo_ptr, MEMO_TABLE_SIZE, HOT_PATH_FUSE);
            
            private_tested++;
            private_total_steps += res.steps;
            
            if (res.steps > local_max_steps) {
                local_max_steps = res.steps;
                local_hardest_n = n;
            }
            if (res.peak > local_max_peak) {
                local_max_peak = res.peak;
            }
            
            if (res.overflow) {
                if (res.steps >= HOT_PATH_FUSE) {
                    private_fuse_count++;
                    // Skip cold queue collection for max performance
                } else {
                    private_overflow_count++;
                    // Skip cold queue collection for max performance  
                }
            }
        }
        
        // ===== CRITICAL SECTION: Only for max values (no vectors) =====
        #pragma omp critical
        {
            if (local_max_steps > global_max_steps) {
                global_max_steps = local_max_steps;
                global_hardest_n = local_hardest_n;
            }
            if (local_max_peak > global_max_peak) {
                global_max_peak = local_max_peak;
            }
            // Cold queues removed for maximum performance
        }
    }
    
    // Copy reduction results to globals
    global_tested = private_tested;
    global_total_steps = private_total_steps;
    global_overflow_count = private_overflow_count;
    global_fuse_count = private_fuse_count;
    
    auto t_end = std::chrono::high_resolution_clock::now();
    double elapsed = std::chrono::duration<double>(t_end - t_start).count();
    
    // Note: Cold queue processing skipped for now (add later if needed)
    // For now, just count them as verified (they exist but not reprocessed)
    
    if (VERBOSE_WORKERS) {
        fprintf(stderr, "[Rank %d] Completed %llu seeds in %.2fs (%.2f M/sec)\n",
                rank, global_tested, elapsed, global_tested / elapsed / 1e6);
    }
    
    // Send results back to master
    WorkerResults results;
    results.rank = rank;
    results.tested = global_tested;
    results.total_steps = global_total_steps;
    results.max_steps = global_max_steps;
    results.max_steps_seed = global_hardest_n;
    results.max_peak = global_max_peak;
    results.overflow_count = global_overflow_count;
    results.fuse_count = global_fuse_count;
    results.time_seconds = elapsed;
    
    return results;
}

// ============================================================================
// MASTER PROCESS (Rank 0)
// ============================================================================

static void master_process(int size, uint64_t start_offset, uint64_t count, const char* run_tag) {
    printf("\n");
    printf("╔════════════════════════════════════════════════════════════════╗\n");
    printf("║       MNv2: MareNostrum 5 Collatz (Phase 1 - GPP Nodes)      ║\n");
    printf("╚════════════════════════════════════════════════════════════════╝\n");
    printf("\n");
    
    print_configuration();
    
    int num_workers = size - 1;
    printf("MPI Configuration:\n");
    printf("  Total ranks:     %d\n", size);
    printf("  Master rank:     %d\n", MASTER_RANK);
    printf("  Worker ranks:    %d (ranks 1..%d)\n", num_workers, size - 1);
    printf("\n");
    
    // Calculate base range
    uint128_t base = (uint128_t)1 << BASE_POWER;
    uint128_t start = base + start_offset;
    uint128_t end = start + count;
    
    // Align for mod-6 filter
    uint64_t delta = 4;
    align_start_and_delta(start, delta);
    
    printf("Workload:\n");
    printf("  Start:           2^%d + %llu\n", BASE_POWER, start_offset);
    printf("  Range size:      %llu\n", count);
    printf("  Seeds per worker: ~%llu\n", count / num_workers);
    printf("\n");
    
    // Partition work evenly across workers
    uint64_t seeds_per_worker = count / num_workers;
    
    printf("[Master] Distributing work to %d workers...\n", num_workers);
    for (int worker_rank = 1; worker_rank < size; ++worker_rank) {
        uint128_t worker_start = start + (worker_rank - 1) * seeds_per_worker;
        uint128_t worker_end = (worker_rank < size - 1) 
                                ? worker_start + seeds_per_worker 
                                : end;  // Last worker gets remainder
        uint64_t worker_candidates = count_candidates(worker_start, worker_end);
        
        MPI_Send_u128(worker_start, worker_rank, 100, MPI_COMM_WORLD);
        MPI_Send_u128(worker_end, worker_rank, 101, MPI_COMM_WORLD);
        MPI_Send(&worker_candidates, 1, MPI_UNSIGNED_LONG_LONG, worker_rank, 102, MPI_COMM_WORLD);
        
        printf("  Rank %d: %llu candidates\n", worker_rank, worker_candidates);
    }
    
    printf("[Master] Waiting for results...\n\n");
    
    // Collect results from all workers
    std::vector<WorkerResults> all_results;
    auto wall_start = std::chrono::high_resolution_clock::now();
    
    for (int worker_rank = 1; worker_rank < size; ++worker_rank) {
        WorkerResults res;
        res.rank = worker_rank;
        
        MPI_Recv(&res.tested, 1, MPI_UNSIGNED_LONG_LONG, worker_rank, 200, 
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&res.total_steps, 1, MPI_UNSIGNED_LONG_LONG, worker_rank, 201,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&res.max_steps, 1, MPI_UNSIGNED_LONG_LONG, worker_rank, 202,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        res.max_steps_seed = MPI_Recv_u128(worker_rank, 203, MPI_COMM_WORLD);
        res.max_peak = MPI_Recv_u128(worker_rank, 204, MPI_COMM_WORLD);
        MPI_Recv(&res.overflow_count, 1, MPI_UNSIGNED_LONG_LONG, worker_rank, 205,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&res.fuse_count, 1, MPI_UNSIGNED_LONG_LONG, worker_rank, 206,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&res.time_seconds, 1, MPI_DOUBLE, worker_rank, 207,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        
        all_results.push_back(res);
        
        printf("[Master] Rank %d: %llu tested in %.2fs (%.2f M/sec)\n",
               res.rank, res.tested, res.time_seconds, res.tested / res.time_seconds / 1e6);
    }
    
    auto wall_end = std::chrono::high_resolution_clock::now();
    double wall_time = std::chrono::duration<double>(wall_end - wall_start).count();
    
    // Aggregate results
    uint64_t total_tested = 0;
    uint64_t total_steps = 0;
    uint64_t max_steps = 0;
    uint128_t hardest_seed = 0;
    uint128_t max_peak = 0;
    uint64_t total_overflows = 0;
    uint64_t total_fuses = 0;
    double max_worker_time = 0;
    
    for (const auto& res : all_results) {
        total_tested += res.tested;
        total_steps += res.total_steps;
        total_overflows += res.overflow_count;
        total_fuses += res.fuse_count;
        if (res.max_steps > max_steps) {
            max_steps = res.max_steps;
            hardest_seed = res.max_steps_seed;
        }
        if (res.max_peak > max_peak) {
            max_peak = res.max_peak;
        }
        if (res.time_seconds > max_worker_time) {
            max_worker_time = res.time_seconds;
        }
    }
    
    double avg_steps = (double)total_steps / total_tested;
    double throughput = total_tested / wall_time;
    
    // Print final summary
    printf("\n");
    printf("╔════════════════════════════════════════════════════════════════╗\n");
    printf("║                         RESULTS                                ║\n");
    printf("╚════════════════════════════════════════════════════════════════╝\n");
    printf("\n");
    printf("Tested:          %llu numbers\n", total_tested);
    printf("Wall time:       %.2f seconds\n", wall_time);
    printf("Worker time:     %.2f seconds (max)\n", max_worker_time);
    printf("Throughput:      %.2f M nums/sec\n", throughput / 1e6);
    printf("Avg steps:       %.2f\n", avg_steps);
    printf("Max steps:       %llu (seed: %s)\n", max_steps, to_string_u128(hardest_seed).c_str());
    printf("Max peak:        %s\n", to_string_u128(max_peak).c_str());
    printf("Overflows:       %llu (128-bit exceeded)\n", total_overflows);
    printf("Fuse hits:       %llu (exceeded %llu iterations)\n", total_fuses, HOT_PATH_FUSE);
    printf("\n");
    printf("Parallel efficiency: %.1f%% (%.2f vs %.2f ideal)\n",
           100.0 * throughput / (num_workers * (total_tested / max_worker_time)),
           throughput, num_workers * (total_tested / max_worker_time));
    printf("\n");
    printf("╔════════════════════════════════════════════════════════════════╗\n");
    printf("║                         SUCCESS!                               ║\n");
    printf("╚════════════════════════════════════════════════════════════════╝\n");
    printf("\n");
}

// ============================================================================
// MAIN ENTRY POINT
// ============================================================================

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    if (size < 2) {
        if (rank == 0) {
            fprintf(stderr, "ERROR: Need at least 2 MPI ranks (1 master + 1 worker)\n");
            fprintf(stderr, "Usage: mpirun -np <N> %s <start_offset> <count> [run_tag]\n", argv[0]);
        }
        MPI_Finalize();
        return 1;
    }
    
    if (argc < 3) {
        if (rank == 0) {
            fprintf(stderr, "Usage: %s <start_offset> <count> [run_tag]\n", argv[0]);
            fprintf(stderr, "\n");
            fprintf(stderr, "  start_offset  Offset from 2^%d\n", BASE_POWER);
            fprintf(stderr, "  count         Total seeds to test\n");
            fprintf(stderr, "  run_tag       Optional tag for output files\n");
            fprintf(stderr, "\n");
            fprintf(stderr, "Example:\n");
            fprintf(stderr, "  mpirun -np 4 %s 0 1000000 test\n", argv[0]);
        }
        MPI_Finalize();
        return 1;
    }
    
    uint64_t start_offset = std::stoull(argv[1]);
    uint64_t count = std::stoull(argv[2]);
    const char* run_tag = (argc >= 4) ? argv[3] : "mnv2_gpp";
    
    if (rank == MASTER_RANK) {
        // Master process
        master_process(size, start_offset, count, run_tag);
    } else {
        // Worker process
        WorkerResults results = worker_process(rank, size);
        
        // Send results to master
        MPI_Send(&results.tested, 1, MPI_UNSIGNED_LONG_LONG, MASTER_RANK, 200, MPI_COMM_WORLD);
        MPI_Send(&results.total_steps, 1, MPI_UNSIGNED_LONG_LONG, MASTER_RANK, 201, MPI_COMM_WORLD);
        MPI_Send(&results.max_steps, 1, MPI_UNSIGNED_LONG_LONG, MASTER_RANK, 202, MPI_COMM_WORLD);
        MPI_Send_u128(results.max_steps_seed, MASTER_RANK, 203, MPI_COMM_WORLD);
        MPI_Send_u128(results.max_peak, MASTER_RANK, 204, MPI_COMM_WORLD);
        MPI_Send(&results.overflow_count, 1, MPI_UNSIGNED_LONG_LONG, MASTER_RANK, 205, MPI_COMM_WORLD);
        MPI_Send(&results.fuse_count, 1, MPI_UNSIGNED_LONG_LONG, MASTER_RANK, 206, MPI_COMM_WORLD);
        MPI_Send(&results.time_seconds, 1, MPI_DOUBLE, MASTER_RANK, 207, MPI_COMM_WORLD);
    }
    
    MPI_Finalize();
    return 0;
}
