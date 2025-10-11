/*
 * Collatz Conjecture - Version 1.1 (Accelerated - CTZ + Bitwise)
 * 
 * Description: Collapse even runs using count-trailing-zeros + bitwise operations
 * Starting from: 2^71 and pushing tested limits
 * 
 * Algorithm optimizations:
 * - Collapse all consecutive divisions by 2 into single shift: n >> ctz(n)
 * - Use bitwise AND instead of modulo: n & 1
 * - Accelerated odd step: compute m = 3n+1, then strip all 2s from m
 * - Each iteration handles multiple original Collatz steps
 * 
 * Key changes from V1.0:
 * - ctz_u128() to count trailing zeros in 128-bit numbers
 * - Single shift operation instead of multiple divisions
 * - Reduces loop iterations dramatically
 * - Overflow detection before 3n+1 operation
 * 
 * Expected: Multi-× speedup over V1.0 baseline
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

// Accelerated Collatz with CTZ - NO PRINTS IN HOT PATH
// Note: Step limit is a safety fuse, NOT real cycle detection
CollatzResult compute_collatz(uint128_t n, uint64_t max_steps = 100000000) {
    CollatzResult r;
    r.number = n;
    r.steps = 0;
    r.reached_one = false;
    r.max_steps_exceeded = false;
    r.overflow_detected = false;

    while (n != 1) {
        // Collapse any even run in one shot
        int z = ctz_u128(n);
        if (z) {
            n >>= z;
            r.steps += z;
            if (r.steps > max_steps) { r.max_steps_exceeded = true; return r; }
            continue;
        }

        // n is odd here: do 3n+1 with overflow guard
        if (n > (U128_MAX - 1) / 3) { r.overflow_detected = true; return r; }
        uint128_t m = 3*n + 1;

        // Strip all twos from m
        int k = ctz_u128(m);
        n = m >> k;

        // Count 1 step for (3n+1) and k steps for the k divisions by 2
        r.steps += 1 + k;
        if (r.steps > max_steps) { r.max_steps_exceeded = true; return r; }
    }

    r.reached_one = true;
    return r;
}

int main(int argc, char* argv[]) {
    std::cout << "=== Collatz Conjecture - Version 1.1 (Accelerated) ===" << std::endl;
    std::cout << "CTZ optimization + bitwise operations" << std::endl;
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
    std::cout << "Start: ";
    print_uint128(start);
    std::cout << std::endl;
    std::cout << "End: ";
    print_uint128(end);
    std::cout << std::endl;
    std::cout << "Max steps safety fuse: " << max_steps << std::endl;
    std::cout << "NOTE: CLI args limited to uint64_t (offset/count < 2^64)" << std::endl;
    std::cout << std::endl;
    
    // Timing
    auto start_time = std::chrono::high_resolution_clock::now();
    
    uint64_t total_tested = 0;
    uint64_t longest_path_steps = 0;
    uint128_t number_with_max_steps = 0;
    uint64_t total_steps = 0;
    
    // Main computation loop
    for (uint128_t n = start; n < end; n++) {
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
        
        // Progress report every 10000 numbers
        if (total_tested % 10000 == 0) {
            std::cout << "Progress: " << total_tested << " numbers tested" << std::endl;
        }
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
    
    return 0;
}
