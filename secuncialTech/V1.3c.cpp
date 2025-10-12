// V1.3c: Precomputed read-only memo table for thread-safe parallelization
//
// KEY CHANGES FROM V1.3a:
// 1. Precompute phase: Fully populate 2^20 table BEFORE timing starts
//    - Use path-compression backfill during precompute
//    - Ensures all small trajectories are cached
// 2. Compute phase: Read-only memo access (no writes)
//    - Thread-safe, trivial to parallelize with OpenMP/CUDA
//    - Stable cache behavior
// 3. Optional disk save/load: Cache table across runs (steps_2p20.bin)
//
// BENEFITS:
// - Eliminates lazy-fill complexity in hot path
// - Thread-safe memo access (no locks needed)
// - Predictable cache behavior (no cold misses)
// - Foundation for V1.4 (OpenMP) and V1.5 (CUDA)
//
// EXPECTED PERFORMANCE: Similar to V1.3a (maybe slightly faster due to no write overhead)

#include <iostream>
#include <vector>
#include <cstdint>
#include <cstring>
#include <cstdio>
#include <ctime>
#include <sys/stat.h>

// ----- Basic 128-bit utilities -----
typedef __uint128_t uint128_t;

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

// ----- Constants -----
static constexpr uint64_t SAFETY_FUSE = 100000;
static constexpr uint32_t UNKNOWN = UINT32_MAX;
// Safe threshold for 3n+1: (UINT128_MAX - 1) / 3
static constexpr uint128_t MAX_SAFE = ((~(uint128_t)0) - 1) / 3;

// ----- Global memo table (read-only after precompute) -----
static uint32_t SMALL_LIMIT_BITS = 20;  // default 2^20 = 1M entries (4 MB) - mutable via --small-limit
static uint64_t small_limit = 0;
static std::vector<uint32_t> memo;

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
            100.0 * filled / small_limit);
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
    fprintf(stderr, "[SAVE] Wrote %lu entries to %s\n", small_limit, filename);
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
    fprintf(stderr, "[LOAD] Loaded %lu entries from %s\n", small_limit, filename);
    return true;
}

// ----- Collatz computation (READ-ONLY memo access) -----
static CollatzResult compute_collatz_readonly(uint128_t n, uint64_t max_steps = SAFETY_FUSE) {
    CollatzResult res = {0, n, false};
    
    if (n == 0) return res;
    if (n == 1) return res;
    
    uint128_t current = n;
    uint64_t steps = 0;
    
    while (current != 1) {
        // Early-exit check: READ-ONLY lookup
        if (current < small_limit) {
            uint32_t cached = memo[(uint64_t)current];
            if (cached != UNKNOWN) {
                res.steps = steps + cached;
                return res;
            }
        }
        
        // Safety fuse
        if (steps >= max_steps) {
            res.overflow = true;
            res.steps = steps;
            return res;
        }
        
        // Collatz step with CTZ (bundle odd step + even collapse)
        if ((current & 1) == 0) {
            int shift = ctz_u128(current);
            current >>= shift;
            steps += shift;
        } else {
            // Check for overflow BEFORE computing 3n+1
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

// ----- Main driver -----
int main(int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <start_offset> <count> [--small-limit BITS] [--save FILE] [--load FILE]\n";
        std::cerr << "  Tests numbers 2^71 + start_offset, +4, +6, +4, +6, ... (mod-6 pattern)\n";
        std::cerr << "  --small-limit BITS: Set memo table size to 2^BITS (default 20 = 1M entries)\n";
        std::cerr << "  --save FILE: Save precomputed table to disk\n";
        std::cerr << "  --load FILE: Load precomputed table from disk (if exists)\n";
        return 1;
    }
    
    uint64_t start_offset = std::stoull(argv[1]);
    uint64_t count = std::stoull(argv[2]);
    
    const char* save_file = nullptr;
    const char* load_file = nullptr;
    
    // Parse optional arguments
    for (int i = 3; i < argc; i++) {
        if (strcmp(argv[i], "--small-limit") == 0 && i + 1 < argc) {
            SMALL_LIMIT_BITS = std::stoul(argv[++i]);
        } else if (strcmp(argv[i], "--save") == 0 && i + 1 < argc) {
            save_file = argv[++i];
        } else if (strcmp(argv[i], "--load") == 0 && i + 1 < argc) {
            load_file = argv[++i];
        }
    }
    
    // Sanity check: prevent absurd table sizes
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
        fprintf(stderr, "prefault=%lu\n", prefault_sum);
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
    fprintf(stderr, ") with mod-6 filter (stride %lu↔%lu)\n", delta, delta^6);
    fprintf(stderr, "[MEMO] Table is now READ-ONLY (thread-safe)\n");
    
    // Timing and computation
    struct timespec t_start, t_end;
    clock_gettime(CLOCK_MONOTONIC, &t_start);
    
    uint64_t tested = 0;
    uint64_t total_steps = 0;
    uint64_t max_steps_seen = 0;
    uint128_t max_peak = 0;
    uint128_t hardest_n = 0;
    
    uint128_t n = start;
    
    for (; n < end; ) {
        CollatzResult res = compute_collatz_readonly(n);
        
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
            std::cout << "OVERFLOW at n=";
            print_u128(n);
            std::cout << " (steps=" << res.steps << ")\n";
            break;
        }
        
        n += delta;
        delta ^= 6;  // Branchless toggle: 4 XOR 6 = 2, 2 XOR 6 = 4
        
        // Progress report (power-of-2 for efficiency)
        if ((tested & ((1 << 14) - 1)) == 0 && tested > 0) {
            fprintf(stderr, "[PROGRESS] Tested %lu numbers\n", tested);
        }
    }
    
    clock_gettime(CLOCK_MONOTONIC, &t_end);
    
    // Results
    double elapsed_s = (t_end.tv_sec - t_start.tv_sec) +
                       (t_end.tv_nsec - t_start.tv_nsec) * 1e-9;
    double throughput = tested / elapsed_s;
    double avg_steps = (double)total_steps / tested;
    
    std::cout << "\n=== V1.3c Results (Precomputed Read-Only Memo) ===\n";
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
    
    return 0;
}
