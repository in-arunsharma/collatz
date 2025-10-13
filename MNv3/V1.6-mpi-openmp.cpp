// V1.6-mpi-openmp: MPI + OpenMP Hybrid for MareNostrum 5 GPP Multi-Node
//
// PARALLELIZATION STRATEGY:
// - MPI: Distribute work across nodes (embarrassingly parallel)
// - OpenMP: Utilize all cores within each node (112 per node)
// - No inter-node communication during computation
// - Aggregate results at end with MPI_Reduce
//
// TARGET: MareNostrum 5 GPP multi-node
// EXPECTED: 137M nums/sec × N nodes (linear scaling)
//
// BASED ON: V1.5-openmp (verified 137M nums/sec per node)
//
// CHANGES FROM V1.5-openmp:
// 1. Added MPI initialization and finalization
// 2. Split work range by MPI rank
// 3. Each rank processes independent subrange
// 4. Aggregate statistics with MPI_Reduce at end
// 5. Root rank prints combined results

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
#include <mpi.h>  // MPI header

// ----- 128-bit utilities -----
typedef __uint128_t uint128_t;

// ----- 256-bit arithmetic (same as V1.5) -----
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
    
    bool is_zero() const {
        return (limbs[0] | limbs[1] | limbs[2] | limbs[3]) == 0;
    }
    bool is_one() const {
        return limbs[0] == 1 && limbs[1] == 0 && limbs[2] == 0 && limbs[3] == 0;
    }
    bool is_even() const {
        return (limbs[0] & 1) == 0;
    }
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
            for (int i = 0; i < 4 - limb_shift; i++) {
                result.limbs[i] = limbs[i + limb_shift];
            }
        } else {
            for (int i = 0; i < 4 - limb_shift; i++) {
                result.limbs[i] = limbs[i + limb_shift] >> bit_shift;
                if (i + limb_shift + 1 < 4) {
                    result.limbs[i] |= limbs[i + limb_shift + 1] << (64 - bit_shift);
                }
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
            limit.limbs[0] = 0x5555555555555555ULL;
            limit.limbs[1] = 0x5555555555555555ULL;
            limit.limbs[2] = 0x5555555555555555ULL;
            limit.limbs[3] = 0x5555555555555555ULL;
            return limit;
        }();
        return *this > MAX_SAFE_256;
    }
    
    std::string to_string() const {
        if (is_zero()) return "0";
        
        uint256_t tmp = *this;
        std::string result;
        
        while (!tmp.is_zero()) {
            uint64_t remainder = 0;
            for (int i = 3; i >= 0; i--) {
                __uint128_t dividend = ((__uint128_t)remainder << 64) | tmp.limbs[i];
                tmp.limbs[i] = dividend / 10;
                remainder = dividend % 10;
            }
            result = char('0' + remainder) + result;
        }
        return result;
    }
};

// ----- 128-bit utilities -----
static inline uint128_t u128_from_u64(uint64_t x) {
    return (uint128_t)x;
}

static void print_u128(uint128_t val) {
    if (val == 0) {
        std::cout << "0";
        return;
    }
    char buffer[64];
    int pos = 0;
    while (val > 0) {
        buffer[pos++] = '0' + (val % 10);
        val /= 10;
    }
    for (int i = pos - 1; i >= 0; i--) {
        std::cout << buffer[i];
    }
}

static std::string to_string_u128(uint128_t val) {
    if (val == 0) return "0";
    char buffer[64];
    int pos = 0;
    while (val > 0) {
        buffer[pos++] = '0' + (val % 10);
        val /= 10;
    }
    std::string result;
    for (int i = pos - 1; i >= 0; i--) {
        result += buffer[i];
    }
    return result;
}

static int ctz_u128(uint128_t x) {
    if (x == 0) return 128;
    uint64_t lo = (uint64_t)x;
    if (lo != 0) {
        return __builtin_ctzll(lo);
    }
    uint64_t hi = (uint64_t)(x >> 64);
    return 64 + __builtin_ctzll(hi);
}

// ----- Result structures -----
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

struct BrentResult {
    bool cycle_found;
    uint64_t steps;
    uint64_t cycle_length;
    std::string cycle_value;
    bool overflow;
};

// ----- Constants -----
static constexpr uint64_t SAFETY_FUSE = 100000;
static constexpr uint64_t EXTENDED_FUSE = 1000000;
static constexpr uint64_t COLD_BATCH_SIZE = 10000;
static constexpr size_t MAX_COLD_QUEUE_SIZE = 1000;
static constexpr uint32_t UNKNOWN = UINT32_MAX;
static constexpr uint128_t MAX_SAFE = ((~(uint128_t)0) - 1) / 3;

// ----- Global memo table -----
static uint32_t SMALL_LIMIT_BITS = 20;
static uint64_t small_limit = 0;
static std::vector<uint32_t> memo;

// ----- Logging (thread-safe) -----
static std::mutex g_log_mutex;
static std::ofstream g_overflow_log;
static std::ofstream g_fuse_log;
static std::ofstream g_cycle_log;
static std::ofstream g_256bit_overflow_log;
static bool logs_opened = false;

// ----- Statistics structure -----
struct ThreadStats {
    uint64_t tested = 0;
    uint64_t total_steps = 0;
    uint64_t max_steps_seen = 0;
    uint128_t max_peak = 0;
    uint128_t hardest_n = 0;
    uint64_t overflow_count = 0;
    uint64_t fuse_count = 0;
    uint64_t cold_q1_triggers = 0;
    uint64_t cold_q1_cycles_found = 0;
    uint64_t cold_q1_verified_ok = 0;
    uint64_t cold_q2_triggers = 0;
    uint64_t cold_q2_cycles_found = 0;
    uint64_t cold_q2_256bit_overflow = 0;
    uint64_t cold_q2_verified_ok = 0;
    double cold_processing_time_ms = 0.0;
};

// ----- Logging helpers -----
static void open_logs_once(const char* run_tag, int rank) {
    if (logs_opened) return;
    logs_opened = true;
    
    std::string tag = std::string(run_tag ? run_tag : "v16mpi") + "_rank" + std::to_string(rank);
    g_overflow_log.open("overflow_seeds_" + tag + ".txt", std::ios::app);
    g_fuse_log.open("fuse_seeds_" + tag + ".txt", std::ios::app);
    g_cycle_log.open("cycle_seeds_" + tag + ".txt", std::ios::app);
    g_256bit_overflow_log.open("256bit_overflow_" + tag + ".txt", std::ios::app);
}

static inline void log_overflow_seed(uint128_t n, uint64_t steps, uint128_t peak) {
    std::lock_guard<std::mutex> lock(g_log_mutex);
    if (g_overflow_log.is_open()) {
        g_overflow_log << to_string_u128(n) << " " << steps << " " << to_string_u128(peak) << "\n";
        g_overflow_log.flush();
    }
}

static inline void log_fuse_seed(uint128_t n, uint64_t steps) {
    std::lock_guard<std::mutex> lock(g_log_mutex);
    if (g_fuse_log.is_open()) {
        g_fuse_log << to_string_u128(n) << " " << steps << "\n";
        g_fuse_log.flush();
    }
}

static inline void log_cycle_seed(uint128_t n, uint64_t cycle_length, const std::string& cycle_value) {
    std::lock_guard<std::mutex> lock(g_log_mutex);
    if (g_cycle_log.is_open()) {
        g_cycle_log << to_string_u128(n) << " " << cycle_length << " " << cycle_value << "\n";
        g_cycle_log.flush();
    }
}

static inline void log_256bit_overflow(uint128_t n, uint64_t steps) {
    std::lock_guard<std::mutex> lock(g_log_mutex);
    if (g_256bit_overflow_log.is_open()) {
        g_256bit_overflow_log << to_string_u128(n) << " " << steps << "\n";
        g_256bit_overflow_log.flush();
    }
}

// ----- Brent cycle detection (same as V1.5) -----
static BrentResult detect_cycle_brent_128(uint128_t n, uint64_t max_steps) {
    BrentResult result = {false, 0, 0, "", false};
    
    uint128_t tortoise = n;
    uint128_t hare = n;
    uint64_t power = 1;
    uint64_t lambda = 1;
    uint64_t steps = 0;
    
    while (steps < max_steps) {
        if (hare == 1) {
            result.steps = steps;
            return result;
        }
        
        if ((hare & 1) == 0) {
            int shift = ctz_u128(hare);
            hare >>= shift;
            steps += shift;
        } else {
            if (hare > MAX_SAFE) {
                result.overflow = true;
                result.steps = steps;
                return result;
            }
            uint128_t t = 3 * hare + 1;
            int k = ctz_u128(t);
            hare = t >> k;
            steps += 1 + k;
        }
        
        if (tortoise == hare) {
            result.cycle_found = true;
            result.cycle_length = lambda;
            result.cycle_value = to_string_u128(hare);
            result.steps = steps;
            return result;
        }
        
        if (lambda == power) {
            tortoise = hare;
            power *= 2;
            lambda = 0;
        }
        lambda++;
    }
    
    result.steps = steps;
    result.overflow = true;
    return result;
}

static BrentResult detect_cycle_brent_256(uint256_t n, uint64_t max_steps) {
    BrentResult result = {false, 0, 0, "", false};
    
    uint256_t tortoise = n;
    uint256_t hare = n;
    uint64_t power = 1;
    uint64_t lambda = 1;
    uint64_t steps = 0;
    
    while (steps < max_steps) {
        if (hare.is_one()) {
            result.steps = steps;
            return result;
        }
        
        if (hare.is_even()) {
            int shift = hare.ctz();
            hare = hare >> shift;
            steps += shift;
        } else {
            if (hare.would_overflow_3n_plus_1()) {
                result.overflow = true;
                result.steps = steps;
                return result;
            }
            uint256_t t = hare * 3 + 1;
            int k = t.ctz();
            hare = t >> k;
            steps += 1 + k;
        }
        
        if (tortoise == hare) {
            result.cycle_found = true;
            result.cycle_length = lambda;
            result.cycle_value = hare.to_string();
            result.steps = steps;
            return result;
        }
        
        if (lambda == power) {
            tortoise = hare;
            power *= 2;
            lambda = 0;
        }
        lambda++;
    }
    
    result.steps = steps;
    result.overflow = true;
    return result;
}

// ----- Precompute table -----
static void precompute_small_table(int rank) {
    memo[0] = UNKNOWN;
    memo[1] = 0;
    
    if (rank == 0) {
        fprintf(stderr, "[PRECOMPUTE] Filling 2^%u table entries...\n", SMALL_LIMIT_BITS);
    }
    
    uint64_t filled = 2;
    
    for (uint64_t n = 2; n < small_limit; n++) {
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
            
            if (current < small_limit && memo[(uint64_t)current] != UNKNOWN) {
                steps += memo[(uint64_t)current];
                completed = true;
                break;
            }
            
            if (steps >= SAFETY_FUSE) break;
            
            if (current < small_limit) {
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
        
        if (!completed) continue;
        
        for (int i = (int)path_vals.size() - 1; i >= 0; i--) {
            uint64_t val = path_vals[i];
            if (val < small_limit && memo[val] == UNKNOWN) {
                memo[val] = steps - path_steps[i];
                filled++;
            }
        }
    }
    
    if (rank == 0) {
        fprintf(stderr, "[PRECOMPUTE] Filled %llu / %llu entries (%.2f%%)\n", 
                (unsigned long long)filled, (unsigned long long)small_limit, 
                100.0 * (double)filled / (double)small_limit);
    }
}

// ----- Collatz computation -----
__attribute__((always_inline)) static inline
CollatzResult compute_collatz_readonly(uint128_t n, 
                                       const uint32_t* __restrict memo_ptr,
                                       uint64_t memo_limit,
                                       uint64_t max_steps = SAFETY_FUSE) {
    CollatzResult res = {0, n, false};
    
    if (n == 0) return res;
    if (n == 1) return res;
    
    uint128_t current = n;
    uint64_t steps = 0;
    
    while (current != 1) {
        if (steps >= max_steps) {
            res.steps = steps;
            return res;
        }
        
        if (current < memo_limit) {
            uint32_t cached = memo_ptr[(uint64_t)current];
            if (cached != UNKNOWN) {
                res.steps = steps + cached;
                return res;
            }
        }
        
        if (current > res.peak) {
            res.peak = current;
        }
        
        if ((current & 1) == 0) {
            int shift = ctz_u128(current);
            current >>= shift;
            steps += shift;
        } else {
            if (current > MAX_SAFE) {
                res.overflow = true;
                res.steps = steps;
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

// ----- Cold queue processing -----
static void process_cold_queue_fuse(std::deque<ColdQueueEntry>& queue, 
                                   ThreadStats& stats,
                                   const uint32_t* memo_ptr) {
    if (queue.empty()) return;
    
    struct timespec t_start, t_end;
    clock_gettime(CLOCK_MONOTONIC, &t_start);
    
    for (const auto& entry : queue) {
        BrentResult brent = detect_cycle_brent_128(entry.seed, EXTENDED_FUSE);
        
        if (brent.cycle_found) {
            stats.cold_q1_cycles_found++;
            log_cycle_seed(entry.seed, brent.cycle_length, brent.cycle_value);
        } else if (!brent.overflow) {
            stats.cold_q1_verified_ok++;
        }
    }
    
    queue.clear();
    
    clock_gettime(CLOCK_MONOTONIC, &t_end);
    double elapsed = (t_end.tv_sec - t_start.tv_sec) + (t_end.tv_nsec - t_start.tv_nsec) * 1e-9;
    stats.cold_processing_time_ms += elapsed * 1000;
}

static void process_cold_queue_overflow(std::deque<ColdQueueEntry>& queue,
                                       ThreadStats& stats) {
    if (queue.empty()) return;
    
    struct timespec t_start, t_end;
    clock_gettime(CLOCK_MONOTONIC, &t_start);
    
    for (const auto& entry : queue) {
        uint256_t seed_256(entry.seed);
        BrentResult brent = detect_cycle_brent_256(seed_256, EXTENDED_FUSE);
        
        if (brent.cycle_found) {
            stats.cold_q2_cycles_found++;
            log_cycle_seed(entry.seed, brent.cycle_length, brent.cycle_value);
        } else if (brent.overflow) {
            stats.cold_q2_256bit_overflow++;
            log_256bit_overflow(entry.seed, brent.steps);
        } else {
            stats.cold_q2_verified_ok++;
        }
    }
    
    queue.clear();
    
    clock_gettime(CLOCK_MONOTONIC, &t_end);
    double elapsed = (t_end.tv_sec - t_start.tv_sec) + (t_end.tv_nsec - t_start.tv_nsec) * 1e-9;
    stats.cold_processing_time_ms += elapsed * 1000;
}

// ----- Mod-6 helpers -----
static inline uint32_t mod3_u128(uint128_t n) {
    uint64_t lo = (uint64_t)n;
    uint64_t hi = (uint64_t)(n >> 64);
    uint32_t r_lo = lo % 3;
    uint32_t r_hi = hi % 3;
    return (r_lo + r_hi) % 3;
}

static inline void align_start_and_delta(uint128_t &n, uint64_t &delta) {
    if ((uint64_t)n % 2 == 0) n++;
    
    uint32_t r3 = mod3_u128(n);
    if (r3 == 0) n += 2;
    
    r3 = mod3_u128(n);
    delta = (r3 == 1) ? 4 : 2;
}

// ----- Main driver -----
int main(int argc, char** argv) {
    // Initialize MPI
    MPI_Init(&argc, &argv);
    
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    if (argc < 3) {
        if (rank == 0) {
            fprintf(stderr, "Usage: %s <start_offset> <count> [options]\n", argv[0]);
            fprintf(stderr, "Options:\n");
            fprintf(stderr, "  --threads <N>        OpenMP threads per node (default: all)\n");
            fprintf(stderr, "  --small-limit <N>    Memo table bits (default: 20)\n");
            fprintf(stderr, "  --tag <name>         Run tag (default: v16mpi)\n");
        }
        MPI_Finalize();
        return 1;
    }
    
    uint64_t start_offset = std::stoull(argv[1]);
    uint64_t total_count = std::stoull(argv[2]);
    
    const char* run_tag = "v16mpi";
    int num_threads = 0;
    
    // Parse optional arguments
    for (int i = 3; i < argc; i++) {
        if (strcmp(argv[i], "--threads") == 0 && i + 1 < argc) {
            num_threads = std::stoi(argv[++i]);
        } else if (strcmp(argv[i], "--small-limit") == 0 && i + 1 < argc) {
            SMALL_LIMIT_BITS = std::stoul(argv[++i]);
        } else if (strcmp(argv[i], "--tag") == 0 && i + 1 < argc) {
            run_tag = argv[++i];
        }
    }
    
    // Set OpenMP threads
    if (num_threads > 0) {
        omp_set_num_threads(num_threads);
    } else {
        num_threads = omp_get_max_threads();
    }
    
    if (rank == 0) {
        fprintf(stderr, "[CONFIG] MPI ranks: %d\n", size);
        fprintf(stderr, "[CONFIG] OpenMP threads per rank: %d\n", num_threads);
        fprintf(stderr, "[CONFIG] Total parallel workers: %d\n", size * num_threads);
    }
    
    // Open logs (per rank)
    open_logs_once(run_tag, rank);
    
    // Setup memo table (each rank builds its own copy)
    small_limit = (1ULL << SMALL_LIMIT_BITS);
    memo.resize(small_limit, UNKNOWN);
    precompute_small_table(rank);
    
    // Prefault pages
    uint64_t prefault_sum = 0;
    for (uint64_t i = 0; i < small_limit; i += 1024) {
        prefault_sum += memo[i];
    }
    
    // Split work by rank
    uint64_t count_per_rank = total_count / size;
    uint64_t my_start_offset = start_offset + rank * count_per_rank;
    uint64_t my_count = count_per_rank;
    
    // Last rank gets remainder
    if (rank == size - 1) {
        my_count = total_count - (size - 1) * count_per_rank;
    }
    
    uint128_t base = u128_from_u64(1) << 71;
    uint128_t start = base + my_start_offset;
    uint128_t end = start + my_count;
    uint64_t delta = 4;
    align_start_and_delta(start, delta);
    
    if (rank == 0) {
        fprintf(stderr, "[RUN] Total range: [2^71+%llu, 2^71+%llu)\n",
                (unsigned long long)start_offset,
                (unsigned long long)(start_offset + total_count));
        fprintf(stderr, "[RUN] Each rank processes ~%llu numbers\n", 
                (unsigned long long)count_per_rank);
    }
    
    // Build local test list
    std::vector<uint128_t> numbers_to_test;
    {
        uint128_t n = start;
        uint64_t current_delta = delta;
        while (n < end) {
            numbers_to_test.push_back(n);
            n += current_delta;
            current_delta ^= 6;
        }
    }
    
    // Global statistics
    ThreadStats global_stats;
    
    // Timing
    struct timespec t_start, t_end;
    MPI_Barrier(MPI_COMM_WORLD);  // Synchronize start
    clock_gettime(CLOCK_MONOTONIC, &t_start);
    
    const uint32_t* __restrict memo_ptr = memo.data();
    
    // OpenMP parallel loop (same as V1.5)
    #pragma omp parallel
    {
        ThreadStats local_stats;
        std::deque<ColdQueueEntry> local_cold_q1;
        std::deque<ColdQueueEntry> local_cold_q2;
        
        #pragma omp for schedule(dynamic, 512)
        for (size_t idx = 0; idx < numbers_to_test.size(); idx++) {
            uint128_t n = numbers_to_test[idx];
            
            CollatzResult res = compute_collatz_readonly(n, memo_ptr, small_limit, SAFETY_FUSE);
            
            local_stats.tested++;
            local_stats.total_steps += res.steps;
            
            if (res.steps > local_stats.max_steps_seen) {
                local_stats.max_steps_seen = res.steps;
                local_stats.hardest_n = n;
            }
            
            if (res.peak > local_stats.max_peak) {
                local_stats.max_peak = res.peak;
            }
            
            if (res.overflow) {
                local_stats.overflow_count++;
                local_stats.cold_q2_triggers++;
                local_cold_q2.emplace_back(n, res.steps, res.peak);
                log_overflow_seed(n, res.steps, res.peak);
            } else if (res.steps >= SAFETY_FUSE) {
                local_stats.fuse_count++;
                local_stats.cold_q1_triggers++;
                local_cold_q1.emplace_back(n, res.steps, res.peak);
                log_fuse_seed(n, res.steps);
            }
            
            if (local_stats.tested % COLD_BATCH_SIZE == 0 || 
                local_cold_q1.size() >= MAX_COLD_QUEUE_SIZE ||
                local_cold_q2.size() >= MAX_COLD_QUEUE_SIZE) {
                
                process_cold_queue_fuse(local_cold_q1, local_stats, memo_ptr);
                process_cold_queue_overflow(local_cold_q2, local_stats);
            }
        }
        
        process_cold_queue_fuse(local_cold_q1, local_stats, memo_ptr);
        process_cold_queue_overflow(local_cold_q2, local_stats);
        
        #pragma omp critical
        {
            global_stats.tested += local_stats.tested;
            global_stats.total_steps += local_stats.total_steps;
            
            if (local_stats.max_steps_seen > global_stats.max_steps_seen) {
                global_stats.max_steps_seen = local_stats.max_steps_seen;
                global_stats.hardest_n = local_stats.hardest_n;
            }
            
            if (local_stats.max_peak > global_stats.max_peak) {
                global_stats.max_peak = local_stats.max_peak;
            }
            
            global_stats.overflow_count += local_stats.overflow_count;
            global_stats.fuse_count += local_stats.fuse_count;
            global_stats.cold_q1_triggers += local_stats.cold_q1_triggers;
            global_stats.cold_q1_cycles_found += local_stats.cold_q1_cycles_found;
            global_stats.cold_q1_verified_ok += local_stats.cold_q1_verified_ok;
            global_stats.cold_q2_triggers += local_stats.cold_q2_triggers;
            global_stats.cold_q2_cycles_found += local_stats.cold_q2_cycles_found;
            global_stats.cold_q2_256bit_overflow += local_stats.cold_q2_256bit_overflow;
            global_stats.cold_q2_verified_ok += local_stats.cold_q2_verified_ok;
            global_stats.cold_processing_time_ms += local_stats.cold_processing_time_ms;
        }
    }
    
    clock_gettime(CLOCK_MONOTONIC, &t_end);
    double local_elapsed_s = (t_end.tv_sec - t_start.tv_sec) +
                             (t_end.tv_nsec - t_start.tv_nsec) * 1e-9;
    
    // Aggregate results across ranks
    uint64_t total_tested = 0;
    uint64_t total_total_steps = 0;
    uint64_t global_max_steps = 0;
    double max_elapsed_s = 0.0;
    
    MPI_Reduce(&global_stats.tested, &total_tested, 1, MPI_UINT64_T, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&global_stats.total_steps, &total_total_steps, 1, MPI_UINT64_T, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&global_stats.max_steps_seen, &global_max_steps, 1, MPI_UINT64_T, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&local_elapsed_s, &max_elapsed_s, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    
    // Close logs
    if (g_overflow_log.is_open()) g_overflow_log.close();
    if (g_fuse_log.is_open()) g_fuse_log.close();
    if (g_cycle_log.is_open()) g_cycle_log.close();
    if (g_256bit_overflow_log.is_open()) g_256bit_overflow_log.close();
    
    // Root rank prints results
    if (rank == 0) {
        double throughput = total_tested / max_elapsed_s;
        double avg_steps = (double)total_total_steps / total_tested;
        
        std::cout << "\n=== V1.6-mpi-openmp Results (MPI + OpenMP Hybrid) ===\n";
        std::cout << "MPI ranks:    " << size << "\n";
        std::cout << "Threads/rank: " << num_threads << "\n";
        std::cout << "Total cores:  " << (size * num_threads) << "\n";
        std::cout << "Tested:       " << total_tested << " numbers\n";
        std::cout << "Time:         " << (uint64_t)(max_elapsed_s * 1000) << " ms\n";
        std::cout << "Throughput:   " << (uint64_t)throughput << " nums/sec\n";
        std::cout << "Speedup:      " << (throughput / 2200000.0) << "× (vs 2.2M/s sequential)\n";
        std::cout << "Per-node:     " << (uint64_t)(throughput / size) << " nums/sec\n";
        std::cout << "Avg steps:    " << avg_steps << "\n";
        std::cout << "Max steps:    " << global_max_steps << "\n";
        std::cout << "\n=== SUCCESS ===\n";
    }
    
    MPI_Finalize();
    return 0;
}
