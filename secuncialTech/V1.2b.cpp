/*
 * Collatz Conjecture - Version 1.2b (Mod-6 Filtered Seeds)
 * 
 * Description: V1.2 + mod-6 filtering (test only n ≡ 1,5 mod 6)
 * Starting from: 2^71 and pushing tested limits
 * 
 * Algorithm optimizations:
 * - All V1.2 optimizations (CTZ, bitwise, accelerated odd step)
 * - Mod-6 filtering: skip even numbers and multiples of 3
 * - Branch-free iteration with alternating +4,+2 stride
 * - Peak excursion tracking for trajectory analysis
 * 
 * Key changes from V1.2:
 * - align_to_valid(): ensures start is odd and not divisible by 3
 * - Alternating stride (4,2,4,2,...) to hit only 1,5 (mod 6) residues
 * - No modulo in hot loop - pure addition and bit toggling
 * - max_excursion field to track peak value in trajectory
 * 
 * Expected: ~3× speedup over V1.2 (testing only ~33% of numbers)
 * Rationale: Even numbers and multiples of 3 are redundantly covered
 */

#include <iostream>
#include <cstdint>
#include <chrono>
#include <iomanip>
#include <string>

// Type for handling large numbers (at least 2^71)
typedef __uint128_t uint128_t;

// Maximum value for uint128_t
constexpr uint128_t U128_MAX = ~(uint128_t)0;

// Structure to hold results
struct CollatzResult {
    uint128_t number;
    uint64_t steps;
    bool reached_one;
    bool max_steps_exceeded;
    bool overflow_detected;
    uint128_t max_excursion;  // NEW: peak value reached in trajectory
};

// Count trailing zeros for 128-bit unsigned integers
static inline int ctz_u128(uint128_t x) {
    if (!x) return 128;
    uint64_t lo = (uint64_t)x;
    if (lo) return __builtin_ctzll(lo);
    uint64_t hi = (uint64_t)(x >> 64);
    return 64 + __builtin_ctzll(hi);
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

// Accelerated Collatz with CTZ + peak excursion tracking
// Note: Step limit is a safety fuse, NOT real cycle detection
CollatzResult compute_collatz(uint128_t n, uint64_t max_steps = 100000000) {
    CollatzResult r;
    r.number = n;
    r.steps = 0;
    r.reached_one = false;
    r.max_steps_exceeded = false;
    r.overflow_detected = false;
    r.max_excursion = n;  // Initialize with starting value

    while (n != 1) {
        // Collapse any even run in one shot
        int z = ctz_u128(n);
        if (z) {
            n >>= z;
            r.steps += z;
            if (n > r.max_excursion) r.max_excursion = n;
            if (r.steps > max_steps) { r.max_steps_exceeded = true; return r; }
            continue;
        }

        // n is odd here: do 3n+1 with overflow guard
        if (n > (U128_MAX - 1) / 3) { r.overflow_detected = true; return r; }
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

// Fast mod3 for uint128_t using 64-bit arithmetic
// Since 2^64 ≡ 1 (mod 3), we have (hi*2^64 + lo) % 3 == (hi + lo) % 3
static inline uint32_t mod3_u128(uint128_t x) {
    uint64_t lo = (uint64_t)x;
    uint64_t hi = (uint64_t)(x >> 64);
    uint32_t s = (uint32_t)((hi % 3) + (lo % 3));
    return (uint32_t)(s % 3);
}

// Align starting value to first valid seed (odd & not multiple of 3)
// Also sets the correct first stride based on residue class
// Result: n ≡ 1 or 5 (mod 6), delta is 4 or 2 accordingly
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
    std::cout << "=== Collatz Conjecture - Version 1.2b (Mod-6 Filtered) ===" << std::endl;
    std::cout << "CTZ optimization + mod-6 filtering (n ≡ 1,5 mod 6 only)" << std::endl;
    std::cout << std::endl;
    
    // Starting point: 2^71
    uint128_t start = (uint128_t)1 << 71;
    uint128_t end = start + 1000000;  // Test first million numbers from 2^71
    uint64_t max_steps = 100000000;    // Safety fuse (CLI configurable)
    
    // Allow command line arguments for range
    // NOTE: offset and count limited to uint64_t (2^64-1)
    // For larger values, use hex input or implement uint128_t parser
    if (argc >= 2) {
        // offset from 2^71 (uint64_t limit)
        uint64_t offset = std::stoull(argv[1]);
        start = ((uint128_t)1 << 71) + offset;
    }
    if (argc >= 3) {
        uint64_t count = std::stoull(argv[2]);
        end = start + count;
    }
    if (argc >= 4) {
        max_steps = std::stoull(argv[3]);
    }
    
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
    std::cout << "NOTE: CLI args limited to uint64_t (offset/count < 2^64)" << std::endl;
    std::cout << std::endl;
    
    // Timing
    auto start_time = std::chrono::high_resolution_clock::now();
    
    uint64_t total_tested = 0;
    uint64_t longest_path_steps = 0;
    uint128_t number_with_max_steps = 0;
    uint64_t total_steps = 0;
    uint128_t global_max_excursion = 0;
    uint128_t number_with_max_excursion = 0;
    
    // Main computation loop with mod-6 filtering
    // Iterate with alternating stride: +4, +2, +4, +2, ...
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
        
        if (result.steps > longest_path_steps) {
            longest_path_steps = result.steps;
            number_with_max_steps = n;
        }
        
        if (result.max_excursion > global_max_excursion) {
            global_max_excursion = result.max_excursion;
            number_with_max_excursion = n;
        }
        
        // Progress report every 10000 numbers
        if (total_tested % 10000 == 0) {
            std::cout << "Progress: " << total_tested << " numbers tested" << std::endl;
        }
        
        // Next valid seed with branchless alternating stride (no % in hot loop)
        n += delta;
        delta ^= 6;  // Branchless toggle: 4 XOR 6 = 2, 2 XOR 6 = 4
    }
    
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    
    // Guard against divide-by-zero
    if (total_tested == 0 || duration.count() == 0) {
        std::cout << "No numbers tested or zero elapsed time.\n";
        return 0;
    }
    
    // Results
    std::cout << std::endl;
    std::cout << "=== Results ===" << std::endl;
    std::cout << "Numbers tested: " << total_tested << std::endl;
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
