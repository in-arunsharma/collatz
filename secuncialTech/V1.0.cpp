/*
 * Collatz Conjecture - Version 1.0 (Baseline - CORRECTED)
 * 
 * Description: Basic sequential implementation with correctness fixes
 * Starting from: 2^71 and pushing tested limits
 * 
 * Algorithm:
 * - For each number n:
 *   - If n is even: n = n / 2
 *   - If n is odd:  n = 3n + 1 (with overflow check)
 * - Continue until n = 1
 * 
 * Correctness fixes:
 * - Overflow protection for 3n+1 operation
 * - Max steps is a safety fuse, NOT cycle detection
 * - CLI args documented as uint64_t limited
 * 
 * NOT YET IMPLEMENTED (future versions):
 * - Accelerated step: (3n+1)/2 for odd n (always produces even)
 * - Bitwise operations: n&1, n>>1 instead of %, /
 * - Skip even chains: n >> __builtin_ctz(n)
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

// Basic Collatz sequence computation with overflow protection
// Note: Step limit is a safety fuse, NOT real cycle detection
CollatzResult compute_collatz(uint128_t n, uint64_t max_steps = 100000000) {
    CollatzResult result;
    result.number = n;
    result.steps = 0;
    result.reached_one = false;
    result.max_steps_exceeded = false;
    result.overflow_detected = false;
    
    uint128_t original = n;
    
    while (n != 1) {
        // Collatz operation
        if (n % 2 == 0) {
            n = n / 2;  // Even: divide by 2
        } else {
            // Odd: check for overflow before 3*n+1
            if (n > (U128_MAX - 1) / 3) {
                std::cerr << "OVERFLOW detected at ";
                print_uint128(original);
                std::cerr << " after " << result.steps << " steps" << std::endl;
                result.overflow_detected = true;
                return result;
            }
            n = 3 * n + 1;  // Odd: multiply by 3 and add 1
        }
        
        result.steps++;
        
        // Safety fuse: prevent infinite loops (NOT real cycle detection)
        if (result.steps > max_steps) {
            std::cerr << "Warning: Exceeded max steps (" << max_steps << ") for ";
            print_uint128(original);
            std::cerr << " - use --max-steps to increase limit" << std::endl;
            result.max_steps_exceeded = true;
            return result;
        }
    }
    
    result.reached_one = true;
    return result;
}

int main(int argc, char* argv[]) {
    std::cout << "=== Collatz Conjecture - Version 1.0 (Baseline) ===" << std::endl;
    std::cout << "Sequential implementation with overflow protection" << std::endl;
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
            std::cerr << "OVERFLOW at ";
            print_uint128(n);
            std::cerr << " - stopping" << std::endl;
            break;
        }
        
        if (result.max_steps_exceeded) {
            std::cerr << "EXCEEDED MAX STEPS at ";
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
