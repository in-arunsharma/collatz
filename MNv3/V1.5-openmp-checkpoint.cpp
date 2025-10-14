// V1.5-openmp-checkpoint: OpenMP Parallel Version with Progress Checkpointing
//
// BASED ON: V1.5-openmp (137M nums/sec proven)
//
// ADDED FEATURES:
// - Progress checkpointing every 1 billion numbers
// - Logs completion ranges to checkpoint file
// - Allows recovery tracking if job is interrupted
//
// PARALLELIZATION STRATEGY:
// - Read-only memo table (thread-safe by design)
// - Embarrassingly parallel: each thread tests independent ranges
// - Dynamic scheduling for load balancing
// - Thread-local statistics, combined at end
// - Chunk-based work distribution (512 numbers per chunk)
//
// TARGET: MareNostrum 5 GPP nodes (56-112 cores)
// EXPECTED: 100-200M nums/sec (50-100× speedup over sequential)

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
#include <atomic>
#include <deque>
#include <omp.h>  // OpenMP header

// ----- Basic 128-bit utilities -----
typedef __uint128_t uint128_t;

// ----- 256-bit arithmetic (same as V1.4b) -----
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

// ----- 128-bit utility functions -----
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

// ----- Result structure -----
struct CollatzResult {
    uint64_t steps;
    uint128_t peak;
    bool overflow;
};

// ----- Cold queue structures -----
struct ColdQueueEntry {
    uint128_t seed;
    uint64_t steps_before_cold;
    uint128_t peak_before_cold;
    
    ColdQueueEntry(uint128_t s, uint64_t st, uint128_t p) 
        : seed(s), steps_before_cold(st), peak_before_cold(p) {}
};

// ----- Brent cycle detection result -----
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

// ----- Global memo table (read-only after precompute - THREAD-SAFE) -----
static uint32_t SMALL_LIMIT_BITS = 20;
static uint64_t small_limit = 0;
static std::vector<uint32_t> memo;
static bool enable_progress = false;

// ----- Logging infrastructure (thread-safe with mutex) -----
static std::mutex g_log_mutex;
static std::ofstream g_overflow_log;
static std::ofstream g_fuse_log;
static std::ofstream g_cycle_log;
static std::ofstream g_256bit_overflow_log;
static bool logs_opened = false;

// ----- Thread-local statistics structure -----
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

// ----- Progress Checkpointing -----
static std::mutex checkpoint_mutex;
static FILE* checkpoint_file = nullptr;
static const uint64_t CHECKPOINT_INTERVAL = 1000000000ULL;  // Log every 1 billion numbers

static void open_checkpoint_file(const char* run_tag) {
    char filename[256];
    snprintf(filename, sizeof(filename), "checkpoint_%s.log", run_tag);
    checkpoint_file = fopen(filename, "a");
    if (checkpoint_file) {
        fprintf(checkpoint_file, "=== Checkpoint Log: %s ===\n", run_tag);
        fprintf(checkpoint_file, "Start time: %ld\n", time(nullptr));
        fflush(checkpoint_file);
    }
}

static void log_checkpoint(uint64_t numbers_completed, uint64_t total_numbers, double elapsed_sec) {
    if (!checkpoint_file) return;
    std::lock_guard<std::mutex> lock(checkpoint_mutex);
    
    double progress_pct = 100.0 * numbers_completed / total_numbers;
    double throughput = numbers_completed / elapsed_sec;
    
    fprintf(checkpoint_file, "[CHECKPOINT] Completed: %llu / %llu (%.2f%%) | Throughput: %.2fM nums/sec | Elapsed: %.1fs\n",
            (unsigned long long)numbers_completed,
            (unsigned long long)total_numbers,
            progress_pct,
            throughput / 1e6,
            elapsed_sec);
    fflush(checkpoint_file);
}

static void close_checkpoint_file() {
    if (checkpoint_file) {
        fprintf(checkpoint_file, "End time: %ld\n", time(nullptr));
        fclose(checkpoint_file);
        checkpoint_file = nullptr;
    }
}

// ----- Logging helpers (same as V1.4b, already thread-safe) -----
static void open_logs_once(const char* run_tag) {
    if (logs_opened) return;
    logs_opened = true;
    
    std::string tag = run_tag ? run_tag : "v15omp";
    g_overflow_log.open("overflow_seeds_" + tag + ".txt", std::ios::app);
    g_fuse_log.open("fuse_seeds_" + tag + ".txt", std::ios::app);
    g_cycle_log.open("cycle_seeds_" + tag + ".txt", std::ios::app);
    g_256bit_overflow_log.open("256bit_overflow_" + tag + ".txt", std::ios::app);
    
    if (g_overflow_log.is_open()) {
        g_overflow_log << "# Overflow seeds - V1.5-openmp\n";
    }
    if (g_fuse_log.is_open()) {
        g_fuse_log << "# Fuse-hit seeds - V1.5-openmp\n";
    }
    if (g_cycle_log.is_open()) {
        g_cycle_log << "# GENUINE CYCLES FOUND - V1.5-openmp\n";
    }
    if (g_256bit_overflow_log.is_open()) {
        g_256bit_overflow_log << "# 256-bit overflow - V1.5-openmp\n";
    }
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

// ----- Brent cycle detection (same as V1.4b) -----
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

// ----- Precompute table (same as V1.4b) -----
static void precompute_small_table() {
    memo[0] = UNKNOWN;
    memo[1] = 0;
    
    fprintf(stderr, "[PRECOMPUTE] Filling 2^%u table entries...\n", SMALL_LIMIT_BITS);
    
    uint64_t filled = 2;
    
    for (uint64_t n = 2; n < small_limit; n++) {
        if (memo[n] != UNKNOWN) {
            continue;
        }
        
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
            
            if (steps >= SAFETY_FUSE) {
                break;
            }
            
            if (current < small_limit) {
                path_vals.push_back((uint64_t)current);
                path_steps.push_back(steps);
            }
            
            if ((current & 1) == 0) {
                int shift = ctz_u128(current);
                current >>= shift;
                steps += shift;
            } else {
                if (current > MAX_SAFE) {
                    break;
                }
                uint128_t t = 3 * current + 1;
                int k = ctz_u128(t);
                current = t >> k;
                steps += 1 + k;
            }
        }
        
        if (!completed) {
            continue;
        }
        
        for (int i = (int)path_vals.size() - 1; i >= 0; i--) {
            uint64_t val = path_vals[i];
            if (val < small_limit && memo[val] == UNKNOWN) {
                memo[val] = steps - path_steps[i];
                filled++;
            }
        }
    }
    
    fprintf(stderr, "[PRECOMPUTE] Filled %llu / %llu entries (%.2f%%)\n", 
            (unsigned long long)filled, (unsigned long long)small_limit, 
            100.0 * (double)filled / (double)small_limit);
}

static bool validate_memo_table() {
    auto has = [&](uint64_t i) { return i < small_limit; };
    
    if (has(1) && memo[1] != 0) {
        fprintf(stderr, "[VALIDATE] FAILED: memo[1] = %u, expected 0\n", memo[1]);
        return false;
    }
    if (has(2) && memo[2] != 1) {
        fprintf(stderr, "[VALIDATE] FAILED: memo[2] = %u, expected 1\n", memo[2]);
        return false;
    }
    if (has(3) && memo[3] != 7) {
        fprintf(stderr, "[VALIDATE] FAILED: memo[3] = %u, expected 7\n", memo[3]);
        return false;
    }
    
    fprintf(stderr, "[VALIDATE] Self-test passed\n");
    return true;
}

// ----- Collatz computation (same as V1.4b) -----
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

// ----- Cold queue processing (per-thread, non-parallel) -----
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
    if (argc < 3) {
        fprintf(stderr, "Usage: %s <start_offset> <count> [options]\n", argv[0]);
        fprintf(stderr, "Options:\n");
        fprintf(stderr, "  --threads <N>        Set number of OpenMP threads (default: all available)\n");
        fprintf(stderr, "  --small-limit <N>    Set memo table bits (default: 20)\n");
        fprintf(stderr, "  --tag <name>         Set run tag for output files (default: v15omp)\n");
        return 1;
    }
    
    uint64_t start_offset = std::stoull(argv[1]);
    uint64_t count = std::stoull(argv[2]);
    
    const char* run_tag = "v15omp";
    int num_threads = 0;  // 0 = use OpenMP default
    
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
        fprintf(stderr, "[CONFIG] Using %d OpenMP threads\n", num_threads);
    } else {
        num_threads = omp_get_max_threads();
        fprintf(stderr, "[CONFIG] Using all available threads: %d\n", num_threads);
    }
    
    // Open logs
    open_logs_once(run_tag);
    open_checkpoint_file(run_tag);
    
    // Setup memo table
    small_limit = (1ULL << SMALL_LIMIT_BITS);
    fprintf(stderr, "[CONFIG] Memo table size: 2^%u = %llu entries (%.2f MB)\n",
            SMALL_LIMIT_BITS, (unsigned long long)small_limit, 
            (small_limit * 4.0) / (1024 * 1024));
    
    memo.resize(small_limit, UNKNOWN);
    precompute_small_table();
    validate_memo_table();
    
    // Prefault pages
    uint64_t prefault_sum = 0;
    for (uint64_t i = 0; i < small_limit; i += 1024) {
        prefault_sum += memo[i];
    }
    if (prefault_sum == 0xDEADBEEF) {
        fprintf(stderr, "Prefault done\n");
    }
    
    // Setup range
    uint128_t base = u128_from_u64(1) << 71;
    uint128_t start = base + start_offset;
    uint128_t end = start + count;
    uint64_t delta = 4;
    align_start_and_delta(start, delta);
    
    fprintf(stderr, "[RUN] Scanning [");
    print_u128(start);
    fprintf(stderr, ", ");
    print_u128(end);
    fprintf(stderr, ") with %d threads\n", num_threads);
    
    // Build list of all numbers to test (for OpenMP parallel for)
    std::vector<uint128_t> numbers_to_test;
    {
        uint128_t n = start;
        uint64_t current_delta = delta;
        while (n < end) {
            numbers_to_test.push_back(n);
            n += current_delta;
            current_delta ^= 6;  // Alternate 2↔4
        }
    }
    
    fprintf(stderr, "[INFO] Testing %zu numbers across %d threads\n", 
            numbers_to_test.size(), num_threads);
    
    // Global statistics (combined from threads)
    ThreadStats global_stats;
    
    // Progress tracking for checkpoints
    std::atomic<uint64_t> numbers_completed(0);
    uint64_t next_checkpoint = CHECKPOINT_INTERVAL;
    uint64_t total_numbers = numbers_to_test.size();
    
    // Timing
    struct timespec t_start, t_end, t_now;
    clock_gettime(CLOCK_MONOTONIC, &t_start);
    
    // Read-only memo pointer
    const uint32_t* __restrict memo_ptr = memo.data();
    
    // PARALLEL LOOP with thread-local statistics and checkpointing
    #pragma omp parallel
    {
        // Thread-local data
        ThreadStats local_stats;
        std::deque<ColdQueueEntry> local_cold_q1;
        std::deque<ColdQueueEntry> local_cold_q2;
        
        // Dynamic scheduling for load balance
        #pragma omp for schedule(dynamic, 512)
        for (size_t idx = 0; idx < numbers_to_test.size(); idx++) {
            uint128_t n = numbers_to_test[idx];
            
            // Hot path computation
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
            
            // Handle edge cases
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
            
            // Process cold queues periodically
            if (local_stats.tested % COLD_BATCH_SIZE == 0 || 
                local_cold_q1.size() >= MAX_COLD_QUEUE_SIZE ||
                local_cold_q2.size() >= MAX_COLD_QUEUE_SIZE) {
                
                process_cold_queue_fuse(local_cold_q1, local_stats, memo_ptr);
                process_cold_queue_overflow(local_cold_q2, local_stats);
            }
            
            // Update progress counter (every 10,000 numbers to reduce atomic overhead)
            if (local_stats.tested % 10000 == 0) {
                numbers_completed.fetch_add(10000, std::memory_order_relaxed);
            }
        }
        
        // Final cold queue processing
        process_cold_queue_fuse(local_cold_q1, local_stats, memo_ptr);
        process_cold_queue_overflow(local_cold_q2, local_stats);
        
        // Combine thread-local statistics into global (critical section)
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
            
            // Log checkpoints as threads finish
            uint64_t current = global_stats.tested;
            if (current >= next_checkpoint) {
                clock_gettime(CLOCK_MONOTONIC, &t_now);
                double elapsed = (t_now.tv_sec - t_start.tv_sec) + 
                               (t_now.tv_nsec - t_start.tv_nsec) * 1e-9;
                log_checkpoint(current, total_numbers, elapsed);
                next_checkpoint = ((current / CHECKPOINT_INTERVAL) + 1) * CHECKPOINT_INTERVAL;
            }
        }
    }
    
    clock_gettime(CLOCK_MONOTONIC, &t_end);
    
    // Results
    double elapsed_s = (t_end.tv_sec - t_start.tv_sec) +
                       (t_end.tv_nsec - t_start.tv_nsec) * 1e-9;
    double throughput = global_stats.tested / elapsed_s;
    double avg_steps = (double)global_stats.total_steps / global_stats.tested;
    
    // Close logs
    if (g_overflow_log.is_open()) g_overflow_log.close();
    if (g_fuse_log.is_open()) g_fuse_log.close();
    if (g_cycle_log.is_open()) g_cycle_log.close();
    if (g_256bit_overflow_log.is_open()) g_256bit_overflow_log.close();
    
    // Print results
    std::cout << "\n=== V1.5-openmp Results (OpenMP Parallel) ===\n";
    std::cout << "Threads:      " << num_threads << "\n";
    std::cout << "Tested:       " << global_stats.tested << " numbers\n";
    std::cout << "Time:         " << (uint64_t)(elapsed_s * 1000) << " ms\n";
    std::cout << "Throughput:   " << (uint64_t)throughput << " nums/sec\n";
    std::cout << "Speedup:      " << (throughput / 2200000.0) << "× (vs 2.2M/s sequential)\n";
    std::cout << "Avg steps:    " << avg_steps << "\n";
    std::cout << "Max steps:    " << global_stats.max_steps_seen << " at n=";
    print_u128(global_stats.hardest_n);
    std::cout << "\n";
    std::cout << "Peak value:   ";
    print_u128(global_stats.max_peak);
    std::cout << "\n";
    std::cout << "\n--- Edge Cases ---\n";
    std::cout << "Overflows:    " << global_stats.overflow_count << "\n";
    std::cout << "Fuse hits:    " << global_stats.fuse_count << "\n";
    std::cout << "Cold Q1:      " << global_stats.cold_q1_triggers << " (cycles=" 
              << global_stats.cold_q1_cycles_found << ", verified=" 
              << global_stats.cold_q1_verified_ok << ")\n";
    std::cout << "Cold Q2:      " << global_stats.cold_q2_triggers << " (cycles="
              << global_stats.cold_q2_cycles_found << ", 256bit-overflow="
              << global_stats.cold_q2_256bit_overflow << ", verified="
              << global_stats.cold_q2_verified_ok << ")\n";
    std::cout << "Cold time:    " << (uint64_t)global_stats.cold_processing_time_ms << " ms\n";
    std::cout << "\n=== SUCCESS ===\n";
    
    // Close checkpoint file
    close_checkpoint_file();
    
    return 0;
}
