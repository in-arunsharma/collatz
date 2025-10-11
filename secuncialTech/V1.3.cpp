/*
 * Collatz Conjecture - Version 1.3 (Early-Exit Memo Table)
 * 
 * Description: V1.2b + early-exit with small read-only memo table
 * Starting from: 2^71 and pushing tested limits
 * 
 * Algorithm optimizations:
 * - All V1.2b optimizations (CTZ, mod-6 filtering, XOR stride)
 * - Early-exit memo: when n < SMALL_LIM, use precomputed steps-to-1
 * - Lazy filling: compute once per small value, cache result
 * - Memory cost: 2^24 entries × 4 bytes ≈ 64 MB (fits in cache)
 * 
 * Key changes from V1.2b:
 * - Global memo table: steps_small[n] = steps from n to 1 (or -1 if unknown)
 * - Early exit when trajectory dips below SMALL_LIM
 * - CLI --small-limit to tune threshold (default 2^24)
 * - Tracks % of seeds that hit memo for telemetry
 * 
 * Expected: Multi-× speedup (most trajectories dip below 2^24 quickly)
 * Rationale: Large seeds dip early → skip long tail to 1
 */

#include <iostream>
#include <cstdint>
#include <chrono>
#include <iomanip>
#include <string>
#include <cstdlib>
#include <cstring>

// Type for handling large numbers (at least 2^71)
typedef __uint128_t uint128_t;

// Maximum value for uint128_t
constexpr uint128_t U128_MAX = ~(uint128_t)0;

// Early-exit memo table
static const uint32_t SMALL_LIM_DEFAULT = 1u << 24;  // 16,777,216 (64 MB)
static uint32_t SMALL_LIM = SMALL_LIM_DEFAULT;
static int32_t* steps_small = nullptr;  // -1 = unknown, else steps-to-1

// Hoist overflow check constant
static const uint128_t MAX_SAFE_N_FOR_3N1 = (U128_MAX - 1) / 3;

// Structure to hold results
struct CollatzResult {
    uint128_t number;
    uint64_t steps;
    bool reached_one;
    bool max_steps_exceeded;
    bool overflow_detected;
    uint128_t max_excursion;
    bool used_memo;  // NEW: track if we hit the memo table
};

// Initialize memo table
void init_small_table() {
    steps_small = (int32_t*)malloc(SMALL_LIM * sizeof(int32_t));
    if (!steps_small) {
        std::cerr << "ERROR: Failed to allocate memo table (" 
                  << (SMALL_LIM * sizeof(int32_t) / (1024*1024)) << " MB)" << std::endl;
        exit(1);
    }
    for (uint32_t i = 0; i < SMALL_LIM; ++i) {
        steps_small[i] = -1;
    }
    steps_small[1] = 0;  // Base case
}

void free_small_table() {
    free(steps_small);
    steps_small = nullptr;
}

// Count trailing zeros for 128-bit unsigned integers
static inline int ctz_u128(uint128_t x) {
    if (!x) return 128;
    uint64_t lo = (uint64_t)x;
    if (lo) return __builtin_ctzll(lo);
    uint64_t hi = (uint64_t)(x >> 64);
    return 64 + __builtin_ctzll(hi);
}

// Fast mod3 for uint128_t using 64-bit arithmetic
// Since 2^64 ≡ 1 (mod 3), we have (hi*2^64 + lo) % 3 == (hi + lo) % 3
static inline uint32_t mod3_u128(uint128_t x) {
    uint64_t lo = (uint64_t)x;
    uint64_t hi = (uint64_t)(x >> 64);
    uint32_t s = (uint32_t)((hi % 3) + (lo % 3));
    return (uint32_t)(s % 3);
}

// Function to print uint128_t numbers
void print_uint128(uint128_t n) {
    if (n == 0) {
        std::cout << "0";
        return;
    }
    
    char buffer[40];
    int i = 39;
    buffer[i] = '\0';
    
    while (n > 0 && i > 0) {
        buffer[--i] = '0' + (n % 10);
        n /= 10;
    }
    
    std::cout << &buffer[i];
}

// Accelerated Collatz with CTZ + early-exit memo + peak excursion tracking
CollatzResult compute_collatz(uint128_t n, uint64_t max_steps = 100000000) {
    CollatzResult r;
    r.number = n;
    r.steps = 0;
    r.reached_one = false;
    r.max_steps_exceeded = false;
    r.overflow_detected = false;
    r.max_excursion = n;
    r.used_memo = false;

    while (n != 1) {
        // Early-exit: if n < SMALL_LIM, check memo or compute+cache
        if (n < SMALL_LIM) {
            uint32_t idx = (uint32_t)n;
            int32_t memo = steps_small[idx];
            
            if (memo >= 0) {
                // Cache hit - use precomputed steps
                r.steps += (uint64_t)memo;
                r.reached_one = true;
                r.used_memo = true;
                return r;
            }
            
            // Cache miss - compute once for this small n
            // Stop when we hit another cached value or reach 1
            uint128_t m = n;
            uint64_t add = 0;
            
            while (m != 1 && m >= SMALL_LIM) {
                int z = ctz_u128(m);
                if (z) {
                    m >>= z;
                    add += z;
                    if (m > r.max_excursion) r.max_excursion = m;
                    continue;
                }
                
                // m is odd: do 3m+1 with overflow guard
                if (m > MAX_SAFE_N_FOR_3N1) {
                    r.overflow_detected = true;
                    return r;
                }
                
                uint128_t t = 3*m + 1;
                int k = ctz_u128(t);
                m = t >> k;
                add += 1 + k;
                if (m > r.max_excursion) r.max_excursion = m;
            }
            
            // Now m < SMALL_LIM or m == 1
            uint64_t tail = (m == 1) ? 0 : (uint64_t)steps_small[(uint32_t)m];
            uint64_t total = add + tail;
            
            // Memoize (cache for future lookups)
            if (total < INT32_MAX) {
                steps_small[idx] = (int32_t)total;
            }
            
            r.steps += total;
            r.reached_one = true;
            r.used_memo = true;
            return r;
        }
        
        // n >= SMALL_LIM: continue with accelerated kernel
        int z = ctz_u128(n);
        if (z) {
            n >>= z;
            r.steps += z;
            if (n > r.max_excursion) r.max_excursion = n;
            if (r.steps > max_steps) { r.max_steps_exceeded = true; return r; }
            continue;
        }

        // n is odd here: do 3n+1 with overflow guard
        if (n > MAX_SAFE_N_FOR_3N1) {
            r.overflow_detected = true;
            return r;
        }
        
        uint128_t m = 3*n + 1;

        // Strip all twos from m
        int k = ctz_u128(m);
        n = m >> k;

        // Track peak excursion
        if (n > r.max_excursion) r.max_excursion = n;

        // Count 1 step for (3n+1) and k steps for the k divisions by 2
        r.steps += 1 + k;
        if (r.steps > max_steps) { r.max_steps_exceeded = true; return r; }
    }

    r.reached_one = true;
    return r;
}

// Align starting value to first valid seed (odd & not multiple of 3)
// Also sets the correct first stride based on residue class
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

int main(int argc, char* argv[]) {
    std::cout << "=== Collatz Conjecture - Version 1.3 (Early-Exit Memo) ===" << std::endl;
    std::cout << "CTZ + mod-6 filtering + early-exit memo table" << std::endl;
    std::cout << std::endl;
    
    // Starting point: 2^71
    uint128_t start = (uint128_t)1 << 71;
    uint128_t end = start + 1000000;  // Test first million numbers from 2^71
    uint64_t max_steps = 100000000;    // Safety fuse
    uint32_t small_limit_pow = 24;     // Default: 2^24
    
    // Parse command line arguments
    // Format: ./V1.3 [offset] [count] [max_steps] [--small-limit N]
    for (int i = 1; i < argc; ++i) {
        if (strcmp(argv[i], "--small-limit") == 0 && i+1 < argc) {
            small_limit_pow = (uint32_t)std::stoul(argv[i+1]);
            ++i;
        } else if (i == 1) {
            uint64_t offset = std::stoull(argv[1]);
            start = ((uint128_t)1 << 71) + offset;
        } else if (i == 2) {
            uint64_t count = std::stoull(argv[2]);
            end = start + count;
        } else if (i == 3) {
            max_steps = std::stoull(argv[3]);
        }
    }
    
    // Set memo table size
    SMALL_LIM = 1u << small_limit_pow;
    
    std::cout << "Memo table: 2^" << small_limit_pow << " = " << SMALL_LIM 
              << " entries (" << (SMALL_LIM * sizeof(int32_t) / (1024*1024)) << " MB)" << std::endl;
    std::cout << "Testing range from 2^71 + offset" << std::endl;
    std::cout << "Start (requested): ";
    print_uint128(start);
    std::cout << std::endl;
    
    // Align start to valid seed (n ≡ 1,5 mod 6) and set correct first stride
    uint64_t delta = 4;  // Will be corrected by align_start_and_delta
    align_start_and_delta(start, delta);
    
    std::cout << "Start (aligned): ";
    print_uint128(start);
    std::cout << std::endl;
    std::cout << "End: ";
    print_uint128(end);
    std::cout << std::endl;
    std::cout << "Max steps safety fuse: " << max_steps << std::endl;
    std::cout << "Filter: Testing only n ≡ 1,5 (mod 6) - skips evens & multiples of 3" << std::endl;
    std::cout << std::endl;
    
    // Initialize memo table
    init_small_table();
    
    // Timing
    auto start_time = std::chrono::high_resolution_clock::now();
    
    uint64_t total_tested = 0;
    uint64_t longest_path_steps = 0;
    uint128_t number_with_max_steps = 0;
    uint64_t total_steps = 0;
    uint128_t global_max_excursion = 0;
    uint128_t number_with_max_excursion = 0;
    uint64_t memo_hits = 0;  // Count how many seeds hit the memo
    
    // Main computation loop with mod-6 filtering
    uint128_t n = start;
    
    for (; n < end; ) {
        CollatzResult result = compute_collatz(n, max_steps);
        
        if (result.overflow_detected) {
            std::cerr << "OVERFLOW detected at ";
            print_uint128(n);
            std::cerr << " after " << result.steps << " steps - stopping" << std::endl;
            break;
        }
        
        if (result.max_steps_exceeded) {
            std::cerr << "EXCEEDED MAX STEPS (" << max_steps << ") at ";
            print_uint128(n);
            std::cerr << " - stopping" << std::endl;
            break;
        }
        
        if (!result.reached_one) {
            std::cerr << "FAILED to reach 1 for ";
            print_uint128(n);
            std::cerr << std::endl;
            break;
        }
        
        total_tested++;
        total_steps += result.steps;
        if (result.used_memo) memo_hits++;
        
        if (result.steps > longest_path_steps) {
            longest_path_steps = result.steps;
            number_with_max_steps = n;
        }
        
        if (result.max_excursion > global_max_excursion) {
            global_max_excursion = result.max_excursion;
            number_with_max_excursion = n;
        }
        
        // Progress report every 10000 numbers (power-of-2 would be: & 0x3FFF == 0 for 16384)
        if (total_tested % 10000 == 0) {
            std::cout << "Progress: " << total_tested << " numbers tested" << std::endl;
        }
        
        // Next valid seed with branchless alternating stride
        n += delta;
        delta ^= 6;  // Branchless toggle: 4 XOR 6 = 2, 2 XOR 6 = 4
    }
    
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    
    // Free memo table
    free_small_table();
    
    // Guard against divide-by-zero
    if (total_tested == 0 || duration.count() == 0) {
        std::cout << "No numbers tested or zero elapsed time.\n";
        return 0;
    }
    
    // Results
    std::cout << std::endl;
    std::cout << "=== Results ===" << std::endl;
    std::cout << "Numbers tested: " << total_tested << std::endl;
    std::cout << "Memo hits: " << memo_hits << " (" << std::fixed << std::setprecision(1)
              << (memo_hits * 100.0 / total_tested) << "%)" << std::endl;
    std::cout << "Total time: " << duration.count() << " ms" << std::endl;
    std::cout << "Average time per number: " << std::fixed << std::setprecision(6) 
              << (double)duration.count() / total_tested << " ms" << std::endl;
    std::cout << "Numbers per second: " << std::fixed << std::setprecision(2)
              << (total_tested * 1000.0) / duration.count() << std::endl;
    std::cout << std::endl;
    std::cout << "Maximum steps: " << longest_path_steps << std::endl;
    std::cout << "Number with max steps: ";
    print_uint128(number_with_max_steps);
    std::cout << std::endl;
    std::cout << "Average steps: " << std::fixed << std::setprecision(2)
              << (double)total_steps / total_tested << std::endl;
    std::cout << std::endl;
    std::cout << "Maximum excursion: ";
    print_uint128(global_max_excursion);
    std::cout << std::endl;
    std::cout << "Number with max excursion: ";
    print_uint128(number_with_max_excursion);
    std::cout << std::endl;
    
    return 0;
}
