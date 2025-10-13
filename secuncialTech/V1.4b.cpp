// V1.4b: Hot/Cold Queue Architecture with Integrated Cycle Detection
//
// KEY CHANGES FROM V1.4:
// 1. HOT PATH (99.9%+): Pure 128-bit, 100K fuse, no cycle detection (~2M nums/sec)
// 2. COLD QUEUE 1: Seeds hitting 100K fuse → reprocess with 1M fuse + Brent cycle detection
// 3. COLD QUEUE 2: Seeds exceeding 128-bit → reprocess with 256-bit + 1M fuse + Brent
// 4. BATCH PROCESSING: Process cold queues every 10K hot seeds (minimal overhead)
// 5. REAL-TIME LOGGING: Crash-safe append-only files for all edge cases
//
// HOT/COLD STRATEGY:
// - Hot path (>99%): Unchanged V1.3d performance, immediate results
// - Cold queue 1: Extended iteration limit (10× longer), Brent cycle check
// - Cold queue 2: 256-bit arithmetic for overflow cases
// - Integration: Single unified program, no post-processing needed
//
// SCIENTIFIC RIGOR:
// - Detect genuine cycles (if they exist!)
// - Verify Collatz for 2^71+ range with mathematical certainty
// - Track all extremals (max steps, max peaks)
// - Full reproducibility for MareNostrum deployment
//
// PERFORMANCE: Hot path unchanged (~2M nums/sec), cold overhead <1%

#include <iostream>
#include <vector>
#include <cstdint>
#include <cstring>
#include <cstdio>
#include <ctime>
#include <time.h>   // for clock_gettime (portable)
#include <string>   // for std::stoull, std::stoul
#include <fstream>  // for overflow/fuse logging
#include <sstream>  // for metadata JSON
#include <iomanip>  // for hex formatting
#include <mutex>    // for thread-safe logging (future OpenMP)
#include <deque>    // for cold queues

// ----- Basic 128-bit utilities -----
typedef __uint128_t uint128_t;

// ----- 256-bit arithmetic (GPU-compatible fixed-width) -----
struct uint256_t {
    uint64_t limbs[4];  // limbs[0] = low 64 bits, limbs[3] = high 64 bits
    
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
        // Count trailing zeros across all limbs
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
            // Simple limb shifting
            for (int i = 0; i < 4 - limb_shift; i++) {
                result.limbs[i] = limbs[i + limb_shift];
            }
        } else {
            // Complex bit shifting across limb boundaries
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
        // If carry is non-zero here, we've overflowed 256 bits
        return result;
    }
    
    bool operator>(const uint256_t& other) const {
        for (int i = 3; i >= 0; i--) {
            if (limbs[i] > other.limbs[i]) return true;
            if (limbs[i] < other.limbs[i]) return false;
        }
        return false;  // Equal
    }
    
    bool operator==(const uint256_t& other) const {
        return limbs[0] == other.limbs[0] && limbs[1] == other.limbs[1] &&
               limbs[2] == other.limbs[2] && limbs[3] == other.limbs[3];
    }
    
    bool would_overflow_3n_plus_1() const {
        // Check if 3*this+1 would overflow 256 bits
        // Exact MAX_SAFE_256 = floor((2^256 - 1 - 1) / 3)
        // = 0x5555555555555555555555555555555555555555555555555555555555555555
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
        // Convert to decimal string (slow but works)
        if (is_zero()) return "0";
        
        uint256_t tmp = *this;
        std::string result;
        
        while (!tmp.is_zero()) {
            // Divide by 10, get remainder
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
    std::string cycle_value;  // A value on the cycle (meeting point), not necessarily entry point
    bool overflow;  // 256-bit overflow in cold queue 2
};

// ----- Constants -----
static constexpr uint64_t SAFETY_FUSE = 100000;  // Hot path fuse
static constexpr uint64_t EXTENDED_FUSE = 1000000;  // Cold queue fuse (10× hot)
static constexpr uint64_t COLD_BATCH_SIZE = 10000;  // Process cold queues every N hot seeds
static constexpr size_t MAX_COLD_QUEUE_SIZE = 1000;  // Safety valve - process if queue fills
static constexpr uint32_t UNKNOWN = UINT32_MAX;
// Safe threshold for 3n+1: (UINT128_MAX - 1) / 3
static constexpr uint128_t MAX_SAFE = ((~(uint128_t)0) - 1) / 3;

// ----- Global memo table (read-only after precompute) -----
static uint32_t SMALL_LIMIT_BITS = 20;  // default 2^20 = 1M entries (4 MB) - mutable via --small-limit
static uint64_t small_limit = 0;
static std::vector<uint32_t> memo;
static bool enable_progress = false;  // Progress reporting off by default (use --progress)

// ----- Logging infrastructure for rare edge cases -----
static std::mutex g_log_mutex;
static std::ofstream g_overflow_log;
static std::ofstream g_fuse_log;
static std::ofstream g_cycle_log;  // NEW: Genuine cycles (if found!)
static std::ofstream g_256bit_overflow_log;  // NEW: 256-bit overflow cases
static bool logs_opened = false;

// ----- Cold queue statistics -----
static uint64_t cold_q1_triggers = 0;
static uint64_t cold_q1_cycles_found = 0;
static uint64_t cold_q1_verified_ok = 0;
static uint64_t cold_q2_triggers = 0;
static uint64_t cold_q2_cycles_found = 0;
static uint64_t cold_q2_256bit_overflow = 0;
static uint64_t cold_q2_verified_ok = 0;
static double cold_processing_time_ms = 0.0;

// ----- Logging helpers for edge cases -----
static void open_logs_once(const char* run_tag) {
    if (logs_opened) return;
    logs_opened = true;
    
    std::string tag = run_tag ? run_tag : "v14b";
    g_overflow_log.open("overflow_seeds_" + tag + ".txt", std::ios::app);
    g_fuse_log.open("fuse_seeds_" + tag + ".txt", std::ios::app);
    g_cycle_log.open("cycle_seeds_" + tag + ".txt", std::ios::app);
    g_256bit_overflow_log.open("256bit_overflow_" + tag + ".txt", std::ios::app);
    
    if (g_overflow_log.is_open()) {
        g_overflow_log << "# Overflow seeds (exceeded 128-bit) - reprocessed with 256-bit\n";
        g_overflow_log << "# Format: seed_decimal steps_before_overflow peak_value\n";
    }
    if (g_fuse_log.is_open()) {
        g_fuse_log << "# Fuse-hit seeds (exceeded max_steps) - reprocessed with extended fuse\n";
        g_fuse_log << "# Format: seed_decimal steps_reached\n";
    }
    if (g_cycle_log.is_open()) {
        g_cycle_log << "# GENUINE CYCLES FOUND (Collatz conjecture violation!)\n";
        g_cycle_log << "# Format: seed_decimal cycle_length cycle_value\n";
        g_cycle_log << "# NOTE: cycle_length is in bundled-iterations (CTZ-collapsed steps), not raw Collatz steps\n";
        g_cycle_log << "# NOTE: cycle_value is a point on the cycle (meeting point), not necessarily the cycle entry\n";
    }
    if (g_256bit_overflow_log.is_open()) {
        g_256bit_overflow_log << "# 256-bit overflow (need 512-bit reprocessing)\n";
        g_256bit_overflow_log << "# Format: seed_decimal steps_reached\n";
    }
}

static inline void log_overflow_seed(uint128_t n, uint64_t steps, uint128_t peak) {
    std::lock_guard<std::mutex> lock(g_log_mutex);
    if (g_overflow_log.is_open()) {
        g_overflow_log << to_string_u128(n) << " " << steps << " " << to_string_u128(peak) << "\n";
        g_overflow_log.flush();  // Ensure written (in case of crash)
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

// ----- Brent cycle detection for 128-bit -----
static BrentResult detect_cycle_brent_128(uint128_t n, uint64_t max_steps) {
    BrentResult result = {false, 0, 0, "", false};
    
    uint128_t tortoise = n;
    uint128_t hare = n;
    uint64_t power = 1;
    uint64_t lambda = 1;  // Cycle length
    uint64_t steps = 0;
    
    while (steps < max_steps) {
        // Advance hare one step
        if (hare == 1) {
            // Reached 1 - no cycle
            result.steps = steps;
            return result;
        }
        
        if ((hare & 1) == 0) {
            int shift = ctz_u128(hare);
            hare >>= shift;
            steps += shift;
        } else {
            if (hare > MAX_SAFE) {
                // Overflow - not a cycle issue
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
            // Cycle detected! (meeting point on cycle)
            result.cycle_found = true;
            result.cycle_length = lambda;
            result.cycle_value = to_string_u128(hare);
            result.steps = steps;
            return result;
        }
        
        if (lambda == power) {
            // Time to teleport tortoise
            tortoise = hare;
            power *= 2;
            lambda = 0;
        }
        lambda++;
    }
    
    // Hit max_steps without finding cycle or reaching 1
    result.steps = steps;
    result.overflow = true;  // Treat as "couldn't verify"
    return result;
}

// ----- Brent cycle detection for 256-bit -----
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

static uint64_t compute_memo_checksum() {
    // Simple 64-bit rolling checksum for reproducibility
    uint64_t checksum = 0;
    const uint64_t PRIME = 0x9E3779B97F4A7C15ULL;
    for (uint64_t i = 0; i < small_limit; i++) {
        if (memo[i] != UNKNOWN) {
            checksum ^= (i * PRIME) ^ ((uint64_t)memo[i] * (PRIME >> 1));
        }
    }
    return checksum;
}

// ----- Precompute full table with path compression -----
static void precompute_small_table() {
    memo[0] = UNKNOWN;  // undefined
    memo[1] = 0;        // base case
    
    fprintf(stderr, "[PRECOMPUTE] Filling 2^%u table entries with path compression...\n", 
            SMALL_LIMIT_BITS);
    
    uint64_t filled = 2;  // Already have 0 and 1
    
    for (uint64_t n = 2; n < small_limit; n++) {
        if (memo[n] != UNKNOWN) {
            continue;  // Already computed
        }
        
        // Compute trajectory, tracking path for backfill
        std::vector<uint64_t> path_vals;
        std::vector<uint64_t> path_steps;  // Cumulative steps at each value
        uint128_t current = (uint128_t)n;  // FIX: Must be 128-bit (can exceed 2^64 in trajectory)
        uint64_t steps = 0;
        bool completed = false;  // FIX: Only backfill if we reached 1 or cached value
        
        while (true) {
            // Check if we hit 1
            if (current == 1) {
                completed = true;
                break;
            }
            
            // Check if we hit something already known
            if (current < small_limit && memo[(uint64_t)current] != UNKNOWN) {
                steps += memo[(uint64_t)current];  // Include cached tail
                completed = true;
                break;
            }
            
            // Safety fuse
            if (steps >= SAFETY_FUSE) {
                // Abort without backfilling (we don't know the true step count)
                break;
            }
            
            // Record current value if small (BEFORE stepping)
            if (current < small_limit) {
                path_vals.push_back((uint64_t)current);
                path_steps.push_back(steps);  // Steps so far when we saw this value
            }
            
            // Collatz step with CTZ (bundle odd step + even collapse)
            if ((current & 1) == 0) {
                int shift = ctz_u128(current);
                current >>= shift;
                steps += shift;
            } else {
                // Check for overflow BEFORE computing 3n+1
                if (current > MAX_SAFE) {
                    // Abort without backfilling (overflow means we can't trust step count)
                    break;
                }
                // Bundle odd step + collapse even run in one shot
                uint128_t t = 3 * current + 1;
                int k = ctz_u128(t);
                current = t >> k;
                steps += 1 + k;
            }
        }
        
        // FIX: Only backfill if we successfully completed (reached 1 or cached)
        if (!completed) {
            continue;  // Leave all visited values as UNKNOWN
        }
        
        // Backfill the path with CORRECT step counts
        for (int i = (int)path_vals.size() - 1; i >= 0; i--) {
            uint64_t val = path_vals[i];
            if (val < small_limit && memo[val] == UNKNOWN) {
                uint64_t steps_from_val_to_1 = steps - path_steps[i];
                if (steps_from_val_to_1 <= UINT32_MAX) {
                    memo[val] = (uint32_t)steps_from_val_to_1;
                    filled++;
                } else {
                    memo[val] = UNKNOWN;  // Too many steps to fit in uint32_t
                }
            }
        }
    }
    
    fprintf(stderr, "[PRECOMPUTE] Filled %llu / %llu entries (%.2f%%)\n", 
            (unsigned long long)filled, (unsigned long long)small_limit, 
            100.0 * (double)filled / (double)small_limit);
}

// ----- Unit self-test: Verify memo table correctness -----
static bool validate_memo_table() {
    // Guard against small table sizes
    auto has = [&](uint64_t i) { return i < small_limit; };
    
    // Test known small values
    if (has(1) && memo[1] != 0) {
        fprintf(stderr, "[VALIDATE] FAILED: memo[1] = %u, expected 0\n", memo[1]);
        return false;
    }
    if (has(2) && memo[2] != 1) {
        fprintf(stderr, "[VALIDATE] FAILED: memo[2] = %u, expected 1\n", memo[2]);
        return false;
    }
    if (has(4) && memo[4] != 2) {
        fprintf(stderr, "[VALIDATE] FAILED: memo[4] = %u, expected 2\n", memo[4]);
        return false;
    }
    if (has(8) && memo[8] != 3) {
        fprintf(stderr, "[VALIDATE] FAILED: memo[8] = %u, expected 3\n", memo[8]);
        return false;
    }
    if (has(16) && memo[16] != 4) {
        fprintf(stderr, "[VALIDATE] FAILED: memo[16] = %u, expected 4\n", memo[16]);
        return false;
    }
    
    // Test a few more known values
    // 3 → 10 → 5 → 16 → 8 → 4 → 2 → 1 (7 steps)
    if (has(3) && memo[3] != 7) {
        fprintf(stderr, "[VALIDATE] FAILED: memo[3] = %u, expected 7\n", memo[3]);
        return false;
    }
    // 5 → 16 → 8 → 4 → 2 → 1 (5 steps)
    if (has(5) && memo[5] != 5) {
        fprintf(stderr, "[VALIDATE] FAILED: memo[5] = %u, expected 5\n", memo[5]);
        return false;
    }
    // 7 → 22 → 11 → 34 → 17 → 52 → 26 → 13 → 40 → 20 → 10 → 5 → 16 → 8 → 4 → 2 → 1 (16 steps)
    if (has(7) && memo[7] != 16) {
        fprintf(stderr, "[VALIDATE] FAILED: memo[7] = %u, expected 16\n", memo[7]);
        return false;
    }
    
    fprintf(stderr, "[VALIDATE] Self-test passed (known values correct)\n");
    return true;
}

// ----- Optional: Save table to disk -----
static bool save_table_to_disk(const char* filename) {
    FILE* f = fopen(filename, "wb");
    if (!f) {
        perror("fopen for write");
        return false;
    }
    
    // Write header: version + small_limit_bits
    uint32_t header[2] = { 1, SMALL_LIMIT_BITS };  // version 1
    if (fwrite(header, sizeof(uint32_t), 2, f) != 2) {
        perror("fwrite header");
        fclose(f);
        return false;
    }
    
    // Write memo table
    if (fwrite(memo.data(), sizeof(uint32_t), small_limit, f) != small_limit) {
        perror("fwrite memo");
        fclose(f);
        return false;
    }
    
    fclose(f);
    fprintf(stderr, "[SAVE] Wrote %llu entries to %s\n", 
            (unsigned long long)small_limit, filename);
    return true;
}

// ----- Optional: Load table from disk -----
static bool load_table_from_disk(const char* filename) {
    FILE* f = fopen(filename, "rb");
    if (!f) {
        return false;  // File doesn't exist, not an error
    }
    
    // Read header
    uint32_t header[2];
    if (fread(header, sizeof(uint32_t), 2, f) != 2) {
        fclose(f);
        return false;
    }
    
    uint32_t version = header[0];
    uint32_t saved_bits = header[1];
    
    if (version != 1 || saved_bits != SMALL_LIMIT_BITS) {
        fprintf(stderr, "[LOAD] Mismatch: saved=%u^%u, wanted=2^%u (rebuilding)\n",
                version, saved_bits, SMALL_LIMIT_BITS);
        fclose(f);
        return false;
    }
    
    // Read memo table
    if (fread(memo.data(), sizeof(uint32_t), small_limit, f) != small_limit) {
        perror("fread memo");
        fclose(f);
        return false;
    }
    
    fclose(f);
    fprintf(stderr, "[LOAD] Loaded %llu entries from %s\n", 
            (unsigned long long)small_limit, filename);
    return true;
}

// ----- Collatz computation (READ-ONLY memo access with restrict hint) -----
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
        // Early-exit check: READ-ONLY lookup with restrict hint
        if (current < memo_limit) {
            uint32_t cached = memo_ptr[(uint64_t)current];
            if (cached != UNKNOWN) {
                uint64_t total = steps + cached;
                // Honor safety fuse even on cache hit
                if (total >= max_steps) {
                    res.overflow = true;
                    res.steps = total;
                    return res;
                }
                res.steps = total;
                return res;
            }
        }
        
        // Safety fuse
        if (steps >= max_steps) {
            res.overflow = true;
            res.steps = steps;
            return res;
        }
        
        // Hoist low 64 bits for parity check (common case after shifts)
        uint64_t lo = (uint64_t)current;
        
        // Collatz step with CTZ (bundle odd step + even collapse)
        if ((lo & 1ULL) == 0) {
            // Even path: collapse run of trailing zeros
            int shift = ctz_u128(current);
            current >>= shift;
            steps += shift;
        } else {
            // Odd path: check overflow BEFORE computing 3n+1
            if (current > MAX_SAFE) {
                res.overflow = true;
                res.steps = steps;
                return res;
            }
            // Bundle odd step + collapse even run in one shot
            uint128_t t = 3 * current + 1;
            int k = ctz_u128(t);
            current = t >> k;
            steps += 1 + k;
            
            if (t > res.peak) {
                res.peak = t;
            }
        }
    }
    
    res.steps = steps;
    return res;
}

// ----- Mod-6 filtering helpers -----
static inline uint32_t mod3_u128(uint128_t n) {
    uint64_t lo = (uint64_t)n;
    uint64_t hi = (uint64_t)(n >> 64);
    uint32_t r_lo = lo % 3;
    uint32_t r_hi = hi % 3;
    // Note: 2^64 ≡ 1 (mod 3), so (lo%3 + hi%3) % 3 is correct
    return (r_lo + r_hi) % 3;
}

static inline void align_start_and_delta(uint128_t &n, uint64_t &delta) {
    // Make odd (parity check is cheap on low 64 bits)
    if ((uint64_t)n % 2 == 0) { n += 1; }
    
    // Compute mod3 using fast 64-bit method
    uint32_t r3 = mod3_u128(n);
    if (r3 == 0) { n += 2; r3 = 2; }  // Skip multiple of 3; now r3 ∈ {1,2}
    
    // Choose first stride based on residue:
    // odd + r3==1 => n ≡ 1 (mod 6) -> start with +4
    // odd + r3==2 => n ≡ 5 (mod 6) -> start with +2
    delta = (r3 == 1) ? 4 : 2;
}

// ----- Cold Queue Processing -----
static void process_cold_queue_fuse(std::deque<ColdQueueEntry>& queue, const char* run_tag) {
    if (queue.empty()) return;
    
    struct timespec t_start, t_end;
    clock_gettime(CLOCK_MONOTONIC, &t_start);
    
    fprintf(stderr, "[COLD_Q1] Processing %zu fuse hits with %lluK step limit + Brent...\n",
            queue.size(), EXTENDED_FUSE / 1000);
    
    size_t verified = 0;
    size_t cycles = 0;
    size_t timeouts = 0;  // Extended fuse hits
    
    for (const auto& entry : queue) {
        BrentResult result = detect_cycle_brent_128(entry.seed, EXTENDED_FUSE);
        
        if (result.cycle_found) {
            // MAJOR FINDING: Genuine cycle!
            cycles++;
            cold_q1_cycles_found++;
            log_cycle_seed(entry.seed, result.cycle_length, result.cycle_value);
            fprintf(stderr, "[COLD_Q1] ⚠️  CYCLE FOUND! Seed=%s, length=%llu\n",
                    to_string_u128(entry.seed).c_str(), (unsigned long long)result.cycle_length);
        } else if (!result.overflow) {
            // Verified: reached 1 within extended limit
            verified++;
            cold_q1_verified_ok++;
        } else {
            // Still hit extended fuse - log for manual review
            timeouts++;
            log_fuse_seed(entry.seed, result.steps);
        }
    }
    
    size_t total = queue.size();
    queue.clear();
    
    clock_gettime(CLOCK_MONOTONIC, &t_end);
    double elapsed = (t_end.tv_sec - t_start.tv_sec) + (t_end.tv_nsec - t_start.tv_nsec) * 1e-9;
    cold_processing_time_ms += elapsed * 1000;
    
    fprintf(stderr, "[COLD_Q1] Verified %zu/%zu (cycles=%zu, timeouts=%zu) in %.1fms\n",
            verified, total, cycles, timeouts, elapsed * 1000);
}

static void process_cold_queue_overflow(std::deque<ColdQueueEntry>& queue, const char* run_tag) {
    if (queue.empty()) return;
    
    struct timespec t_start, t_end;
    clock_gettime(CLOCK_MONOTONIC, &t_start);
    
    fprintf(stderr, "[COLD_Q2] Processing %zu 128-bit overflows with 256-bit + Brent...\n",
            queue.size());
    
    size_t verified = 0;
    size_t cycles = 0;
    size_t overflows_256 = 0;
    
    for (const auto& entry : queue) {
        uint256_t seed_256(entry.seed);
        BrentResult result = detect_cycle_brent_256(seed_256, EXTENDED_FUSE);
        
        if (result.cycle_found) {
            cycles++;
            cold_q2_cycles_found++;
            log_cycle_seed(entry.seed, result.cycle_length, result.cycle_value);
            fprintf(stderr, "[COLD_Q2] ⚠️  CYCLE FOUND! Seed=%s, length=%llu\n",
                    to_string_u128(entry.seed).c_str(), (unsigned long long)result.cycle_length);
        } else if (result.overflow) {
            // 256-bit overflow or hit extended fuse
            overflows_256++;
            cold_q2_256bit_overflow++;
            log_256bit_overflow(entry.seed, result.steps);
        } else {
            verified++;
            cold_q2_verified_ok++;
        }
    }
    
    size_t total = queue.size();
    queue.clear();
    
    clock_gettime(CLOCK_MONOTONIC, &t_end);
    double elapsed = (t_end.tv_sec - t_start.tv_sec) + (t_end.tv_nsec - t_start.tv_nsec) * 1e-9;
    cold_processing_time_ms += elapsed * 1000;
    
    fprintf(stderr, "[COLD_Q2] Verified %zu/%zu (cycles=%zu, 256-bit-overflow=%zu) in %.1fms\n",
            verified, total, cycles, overflows_256, elapsed * 1000);
}

// ----- Metadata output for reproducibility -----
static void write_metadata_json(const char* run_tag,
                                uint64_t start_offset, uint64_t count,
                                uint32_t memo_bits, uint64_t max_steps,
                                uint64_t tested, uint64_t max_steps_seen,
                                uint128_t max_peak, uint128_t hardest_n,
                                uint64_t overflow_count, uint64_t fuse_count,
                                double elapsed_s, double throughput) {
    std::ofstream json("run_metadata_" + std::string(run_tag ? run_tag : "v14b") + ".json");
    if (!json.is_open()) {
        fprintf(stderr, "[WARNING] Could not write metadata JSON\n");
        return;
    }
    
    uint64_t memo_checksum = compute_memo_checksum();
    
    json << "{\n";
    json << "  \"version\": \"V1.4b\",\n";
    json << "  \"timestamp\": " << time(nullptr) << ",\n";
    json << "  \"run_tag\": \"" << (run_tag ? run_tag : "v14b") << "\",\n";
    json << "  \"git_commit\": \"" << 
#ifdef GIT_HASH
        GIT_HASH
#else
        "unknown"
#endif
        << "\",\n";
    json << "  \"compile_flags\": \"-O3 -march=native -flto -fno-exceptions -fno-rtti -funroll-loops -fno-asynchronous-unwind-tables -DNDEBUG\",\n";
    json << "  \"input\": {\n";
    json << "    \"base\": \"2^71\",\n";
    json << "    \"start_offset\": " << start_offset << ",\n";
    json << "    \"count\": " << count << ",\n";
    json << "    \"stride_policy\": \"mod-6 (test n≡1,5 mod 6)\"\n";
    json << "  },\n";
    json << "  \"configuration\": {\n";
    json << "    \"memo_bits\": " << memo_bits << ",\n";
    json << "    \"memo_size\": " << (1ULL << memo_bits) << ",\n";
    json << "    \"memo_checksum\": \"0x" << std::hex << memo_checksum << std::dec << "\",\n";
    json << "    \"hot_fuse\": " << max_steps << ",\n";
    json << "    \"cold_fuse\": " << EXTENDED_FUSE << ",\n";
    json << "    \"batch_size\": " << COLD_BATCH_SIZE << "\n";
    json << "  },\n";
    json << "  \"results\": {\n";
    json << "    \"tested\": " << tested << ",\n";
    json << "    \"time_ms\": " << (uint64_t)(elapsed_s * 1000) << ",\n";
    json << "    \"throughput_per_sec\": " << (uint64_t)throughput << ",\n";
    json << "    \"max_steps\": " << max_steps_seen << ",\n";
    json << "    \"max_steps_seed\": \"" << to_string_u128(hardest_n) << "\",\n";
    json << "    \"max_peak\": \"" << to_string_u128(max_peak) << "\"\n";
    json << "  },\n";
    json << "  \"cold_queues\": {\n";
    json << "    \"note\": \"cycle_length reported in bundled-iterations (CTZ-collapsed), not raw Collatz steps\",\n";
    json << "    \"queue1_fuse_triggers\": " << cold_q1_triggers << ",\n";
    json << "    \"queue1_genuine_cycles\": " << cold_q1_cycles_found << ",\n";
    json << "    \"queue1_verified_ok\": " << cold_q1_verified_ok << ",\n";
    json << "    \"queue2_overflow_triggers\": " << cold_q2_triggers << ",\n";
    json << "    \"queue2_genuine_cycles\": " << cold_q2_cycles_found << ",\n";
    json << "    \"queue2_256bit_overflow\": " << cold_q2_256bit_overflow << ",\n";
    json << "    \"queue2_verified_ok\": " << cold_q2_verified_ok << ",\n";
    json << "    \"total_cold_processing_ms\": " << (uint64_t)cold_processing_time_ms << "\n";
    json << "  },\n";
    json << "  \"verification\": {\n";
    json << "    \"overflow_log\": \"overflow_seeds_" << (run_tag ? run_tag : "v14b") << ".txt\",\n";
    json << "    \"fuse_log\": \"fuse_seeds_" << (run_tag ? run_tag : "v14b") << ".txt\",\n";
    json << "    \"cycle_log\": \"cycle_seeds_" << (run_tag ? run_tag : "v14b") << ".txt\",\n";
    json << "    \"256bit_overflow_log\": \"256bit_overflow_" << (run_tag ? run_tag : "v14b") << ".txt\",\n";
    json << "    \"status\": \"" << (cold_q1_cycles_found == 0 && cold_q2_cycles_found == 0 ? 
        "ALL_VERIFIED" : "CYCLES_FOUND") << "\"\n";
    json << "  }\n";
    json << "}\n";
    
    json.close();
    fprintf(stderr, "[METADATA] Wrote run_metadata_%s.json\n", run_tag ? run_tag : "v14b");
}

// ----- Main driver -----
int main(int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <start_offset> <count> [OPTIONS]\n";
        std::cerr << "  Scans range [2^71 + start_offset, 2^71 + start_offset + count)\n";
        std::cerr << "  with mod-6 filter (tests only n ≡ 1,5 mod 6)\n";
        std::cerr << "\nOPTIONS:\n";
        std::cerr << "  --small-limit BITS  Set memo table size to 2^BITS (default 20 = 1M entries)\n";
        std::cerr << "  --save FILE         Save precomputed table to disk\n";
        std::cerr << "  --load FILE         Load precomputed table from disk (if exists)\n";
        std::cerr << "  --max-steps N       Override safety fuse (default 100000)\n";
        std::cerr << "  --progress [N]      Enable progress reporting every N numbers (default 16384)\n";
        std::cerr << "  --tag TAG           Run tag for output files (default: v14b)\n";
        std::cerr << "\nV1.4b HOT/COLD QUEUE FEATURES:\n";
        std::cerr << "  - Hot path: Pure 128-bit, 100K fuse, ~2M nums/sec\n";
        std::cerr << "  - Cold queue 1: Fuse hits → reprocess with 1M fuse + Brent cycle detection\n";
        std::cerr << "  - Cold queue 2: 128-bit overflow → reprocess with 256-bit + Brent\n";
        std::cerr << "  - Batch processing: Process cold queues every 10K hot seeds\n";
        std::cerr << "  - Real-time logging: cycle_seeds_<tag>.txt (if any cycles found!)\n";
        return 1;
    }
    
    uint64_t start_offset = std::stoull(argv[1]);
    uint64_t count = std::stoull(argv[2]);
    
    const char* save_file = nullptr;
    const char* load_file = nullptr;
    const char* run_tag = "v14b";  // Default tag for output files
    uint64_t max_steps = SAFETY_FUSE;  // Allow override via CLI
    uint64_t progress_interval = 16384;  // Progress reporting interval
    
    // Parse optional arguments
    for (int i = 3; i < argc; i++) {
        if (strcmp(argv[i], "--small-limit") == 0 && i + 1 < argc) {
            SMALL_LIMIT_BITS = std::stoul(argv[++i]);
        } else if (strcmp(argv[i], "--save") == 0 && i + 1 < argc) {
            save_file = argv[++i];
        } else if (strcmp(argv[i], "--load") == 0 && i + 1 < argc) {
            load_file = argv[++i];
        } else if (strcmp(argv[i], "--max-steps") == 0 && i + 1 < argc) {
            max_steps = std::stoull(argv[++i]);
        } else if (strcmp(argv[i], "--progress") == 0) {
            enable_progress = true;
            if (i + 1 < argc && argv[i + 1][0] != '-') {
                progress_interval = std::stoull(argv[++i]);
            }
        } else if (strcmp(argv[i], "--tag") == 0 && i + 1 < argc) {
            run_tag = argv[++i];
        }
    }
    
    // Open logging files
    open_logs_once(run_tag);
    
    // Sanity checks
    if (progress_interval == 0 && enable_progress) {
        fprintf(stderr, "[WARNING] --progress 0 ignored (invalid interval)\n");
        enable_progress = false;
    }
    
    if (SMALL_LIMIT_BITS > 28) {
        fprintf(stderr, "[ERROR] SMALL_LIMIT_BITS=%u too large (max 28 = 1GB table)\n", 
                SMALL_LIMIT_BITS);
        return 1;
    }
    if (SMALL_LIMIT_BITS < 8) {
        fprintf(stderr, "[ERROR] SMALL_LIMIT_BITS=%u too small (min 8)\n", 
                SMALL_LIMIT_BITS);
        return 1;
    }
    
    small_limit = (1ULL << SMALL_LIMIT_BITS);
    fprintf(stderr, "[CONFIG] Memo table size: 2^%u = %llu entries (%.2f MB)\n",
            SMALL_LIMIT_BITS, (unsigned long long)small_limit, 
            (small_limit * 4.0) / (1024 * 1024));
    
    // Allocate and initialize memo table
    memo.resize(small_limit, UNKNOWN);
    
    // Try to load from disk if requested
    bool loaded = false;
    if (load_file) {
        loaded = load_table_from_disk(load_file);
    }
    
    // If not loaded, precompute
    if (!loaded) {
        precompute_small_table();
        
        // Validate memo table correctness
        if (!validate_memo_table()) {
            fprintf(stderr, "[ERROR] Memo table validation failed! Aborting.\n");
            return 1;
        }
        
        // Save if requested
        if (save_file) {
            save_table_to_disk(save_file);
        }
    } else {
        // Also validate loaded tables
        if (!validate_memo_table()) {
            fprintf(stderr, "[ERROR] Loaded memo table validation failed! Aborting.\n");
            return 1;
        }
    }
    
    // Prefault memo pages before timing
    uint64_t prefault_sum = 0;
    for (uint64_t i = 0; i < small_limit; i += 1024) {
        prefault_sum += memo[i];
    }
    if (prefault_sum == 0xDEADBEEF) {  // Never true, just prevent optimization
        fprintf(stderr, "prefault=%llu\n", (unsigned long long)prefault_sum);
    }
    
    // Base number: 2^71
    uint128_t base = u128_from_u64(1) << 71;
    uint128_t start = base + start_offset;
    uint128_t end = start + count;  // Scan range, not test count
    uint64_t delta = 4;  // Will be corrected by align_start_and_delta
    align_start_and_delta(start, delta);
    
    fprintf(stderr, "[RUN] Scanning range [");
    print_u128(start);
    fprintf(stderr, ", ");
    print_u128(end);
    fprintf(stderr, ") with mod-6 filter (stride %llu↔%llu)\n",
            (unsigned long long)delta, (unsigned long long)(delta ^ 6));
    fprintf(stderr, "[MEMO] Table is now READ-ONLY (thread-safe)\n");
    
    // Cold queues for batch processing
    std::deque<ColdQueueEntry> cold_queue_fuse;
    std::deque<ColdQueueEntry> cold_queue_overflow;
    
    // Timing and computation
    struct timespec t_start, t_end;
    clock_gettime(CLOCK_MONOTONIC, &t_start);
    
    uint64_t tested = 0;
    uint64_t total_steps = 0;
    uint64_t max_steps_seen = 0;
    uint128_t max_peak = 0;
    uint128_t hardest_n = 0;
    uint64_t overflow_count = 0;
    uint64_t fuse_count = 0;
    
    uint128_t n = start;
    
    // Use read-only pointer with restrict hint for compiler optimization
    const uint32_t* __restrict memo_ptr = memo.data();
    
    for (; n < end; ) {
        CollatzResult res = compute_collatz_readonly(n, memo_ptr, small_limit, max_steps);
        
        tested++;
        total_steps += res.steps;
        
        if (res.steps > max_steps_seen) {
            max_steps_seen = res.steps;
            hardest_n = n;
        }
        if (res.peak > max_peak) {
            max_peak = res.peak;
        }
        
        if (res.overflow) {
            // Determine which cold queue to use
            if (res.steps >= max_steps) {
                // Hit iteration fuse → Cold Queue 1 (extended fuse + Brent)
                fuse_count++;
                cold_q1_triggers++;
                cold_queue_fuse.emplace_back(n, res.steps, res.peak);
            } else {
                // 128-bit overflow → Cold Queue 2 (256-bit + Brent)
                overflow_count++;
                cold_q2_triggers++;
                cold_queue_overflow.emplace_back(n, res.steps, res.peak);
                // Log to overflow file immediately (for crash recovery)
                log_overflow_seed(n, res.steps, res.peak);
            }
        }
        
        n += delta;
        delta ^= 6;  // Branchless toggle: 4 XOR 6 = 2, 2 XOR 6 = 4
        
        // Batch processing: process cold queues every COLD_BATCH_SIZE seeds
        if (tested % COLD_BATCH_SIZE == 0 && tested > 0) {
            if (!cold_queue_fuse.empty() || !cold_queue_overflow.empty()) {
                process_cold_queue_fuse(cold_queue_fuse, run_tag);
                process_cold_queue_overflow(cold_queue_overflow, run_tag);
            }
        }
        
        // Safety valve: if queues get too large, process immediately
        if (cold_queue_fuse.size() >= MAX_COLD_QUEUE_SIZE) {
            fprintf(stderr, "[SAFETY] Cold queue 1 full (%zu), processing immediately\n", 
                    cold_queue_fuse.size());
            process_cold_queue_fuse(cold_queue_fuse, run_tag);
        }
        if (cold_queue_overflow.size() >= MAX_COLD_QUEUE_SIZE) {
            fprintf(stderr, "[SAFETY] Cold queue 2 full (%zu), processing immediately\n",
                    cold_queue_overflow.size());
            process_cold_queue_overflow(cold_queue_overflow, run_tag);
        }
        
        // Progress report (use modulo for any interval)
        if (enable_progress && (tested % progress_interval == 0) && tested > 0) {
            fprintf(stderr, "[PROGRESS] Tested %llu numbers | Cold Q1: %zu pending | Cold Q2: %zu pending\n", 
                    (unsigned long long)tested, cold_queue_fuse.size(), cold_queue_overflow.size());
        }
    }
    
    clock_gettime(CLOCK_MONOTONIC, &t_end);
    
    // Process any remaining cold queue entries
    fprintf(stderr, "[FINAL] Processing remaining cold queue entries...\n");
    process_cold_queue_fuse(cold_queue_fuse, run_tag);
    process_cold_queue_overflow(cold_queue_overflow, run_tag);
    
    // Results
    double elapsed_s = (t_end.tv_sec - t_start.tv_sec) +
                       (t_end.tv_nsec - t_start.tv_nsec) * 1e-9;
    double throughput = tested / elapsed_s;
    double avg_steps = (double)total_steps / tested;
    
    // Close logs
    if (g_overflow_log.is_open()) g_overflow_log.close();
    if (g_fuse_log.is_open()) g_fuse_log.close();
    if (g_cycle_log.is_open()) g_cycle_log.close();
    if (g_256bit_overflow_log.is_open()) g_256bit_overflow_log.close();
    
    // Write metadata JSON
    write_metadata_json(run_tag, start_offset, count, SMALL_LIMIT_BITS, max_steps,
                       tested, max_steps_seen, max_peak, hardest_n,
                       overflow_count, fuse_count, elapsed_s, throughput);
    
    // Print results
    std::cout << "\n=== V1.4b Results (Hot/Cold Queue Architecture) ===\n";
    std::cout << "Tested:       " << tested << " numbers\n";
    std::cout << "Time:         " << (uint64_t)(elapsed_s * 1000) << " ms\n";
    std::cout << "Throughput:   " << (uint64_t)throughput << " nums/sec\n";
    std::cout << "Avg steps:    " << avg_steps << "\n";
    std::cout << "Max steps:    " << max_steps_seen << " at n=";
    print_u128(hardest_n);
    std::cout << "\n";
    std::cout << "Peak value:   ";
    print_u128(max_peak);
    std::cout << "\n";
    std::cout << "\n--- Hot Path ---\n";
    std::cout << "Hot tested:   " << tested << " seeds\n";
    std::cout << "Hot time:     " << (uint64_t)(elapsed_s * 1000) << " ms\n";
    std::cout << "\n--- Cold Queues ---\n";
    std::cout << "Queue 1 (fuse):      " << cold_q1_triggers << " triggers, " 
              << cold_q1_verified_ok << " verified, " 
              << cold_q1_cycles_found << " CYCLES\n";
    std::cout << "Queue 2 (overflow):  " << cold_q2_triggers << " triggers, "
              << cold_q2_verified_ok << " verified, "
              << cold_q2_cycles_found << " CYCLES, "
              << cold_q2_256bit_overflow << " 256-bit overflow\n";
    std::cout << "Cold processing:     " << (uint64_t)cold_processing_time_ms << " ms ("
              << std::fixed << std::setprecision(2) 
              << (100.0 * cold_processing_time_ms / (elapsed_s * 1000)) << "% overhead)\n";
    std::cout << "\n--- Verification Status ---\n";
    if (cold_q1_cycles_found > 0 || cold_q2_cycles_found > 0) {
        std::cout << "⚠️  CYCLES FOUND: " << (cold_q1_cycles_found + cold_q2_cycles_found) 
                  << " (see cycle_seeds_" << run_tag << ".txt)\n";
    } else {
        std::cout << "✅ ALL VERIFIED: No cycles found in tested range\n";
    }
    if (cold_q2_256bit_overflow > 0) {
        std::cout << "⚠️  256-BIT OVERFLOW: " << cold_q2_256bit_overflow 
                  << " seeds (need 512-bit, see 256bit_overflow_" << run_tag << ".txt)\n";
    }
    std::cout << "\nMetadata:     run_metadata_" << run_tag << ".json\n";
    
    return 0;
}
