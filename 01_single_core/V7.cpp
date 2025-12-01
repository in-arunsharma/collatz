#include <iostream>
#include <chrono>
#include <cstdint>
#include <cstdio>

using namespace std;
using namespace std::chrono;

__attribute__((always_inline))
inline int ctz128(__uint128_t n) 
{
    uint64_t low = (uint64_t)n;
    if (low != 0) return __builtin_ctzll(low);
    return 64 + __builtin_ctzll((uint64_t)(n >> 64));
}

__attribute__((always_inline))
inline __uint128_t collatz_algo(__uint128_t n, uint64_t& steps) 
{
    __uint128_t result;
    uint64_t low = (uint64_t)n;  // Hoist low bits for parity check
    
    if (low & 1) {
        result = 3 * n + 1;
        int zeros = ctz128(result);
        result >>= zeros;
        steps = zeros + 1;
        return result;
    } else {
        int zeros = ctz128(n);
        steps = zeros;
        n >>= zeros;
        return n;
    }
}

__attribute__((always_inline))
inline int64_t collatz_steps(__uint128_t n) {
    const __uint128_t original = n;
    __uint128_t current = n;
    uint64_t steps = 0;
    const uint64_t MAX_ITERATIONS = 100000;
    uint64_t iterations = 0;
    
    while (current != 1) {
        uint64_t n_steps = 0;
        current = collatz_algo(current, n_steps);
        steps += n_steps;
        
        if (current < original) {
            return steps;
        }
        
        if (++iterations > MAX_ITERATIONS) {
            return -1;
        }
    }
    return steps;
}

void print_uint128(__uint128_t x) {
    if (x > 9) print_uint128(x / 10);
    putchar('0' + x % 10);
}

int main() {
    const __uint128_t START = ((__uint128_t)1 << 71);
    const uint64_t COUNT = 1000000000;
    
    cout << "Collatz V7 - Compiler Hints" << endl;
    cout << "===========================" << endl;
    cout << "Start: 2^71 = "; print_uint128(START); cout << endl;
    cout << "Range: " << COUNT << " numbers" << endl;
    cout << "Optimizations: Mod-6 + __restrict__ + __attribute__((always_inline))" << endl << endl;
    
    uint64_t total_steps = 0;
    uint64_t max_steps = 0;
    __uint128_t max_steps_number = 0;
    uint64_t cycles_found = 0;
    uint64_t numbers_tested = 0;
    
    auto start_time = high_resolution_clock::now();
    
    // Find first odd number that is ≡ 1 or 5 (mod 6)
    __uint128_t n = START;
    if (!(n & 1)) n++;  // Make it odd
    
    // Adjust to first n ≡ 1 or 5 (mod 6)
    uint64_t mod6 = (uint64_t)(n % 6);
    if (mod6 == 3) n += 2;  // Skip from 3 to 5
    
    // Determine starting stride (4 if at 1, 2 if at 5)
    uint64_t stride = (mod6 == 1 || mod6 == 3) ? 4 : 2;
    
    __uint128_t end = START + COUNT;
    
    while (n < end) {
        int64_t steps = collatz_steps(n);
        
        if (steps == -1) {
            cycles_found++;
            cout << "CYCLE DETECTED at: "; print_uint128(n); cout << endl;
        } else {
            total_steps += steps;
            if ((uint64_t)steps > max_steps) {
                max_steps = steps;
                max_steps_number = n;
            }
        }
        
        numbers_tested++;
        n += stride;
        stride = 6 - stride;  // Alternate between 4 and 2
    }
    
    auto end_time = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(end_time - start_time);
    
    cout << "Results:" << endl;
    cout << "--------" << endl;
    cout << "Numbers tested: " << numbers_tested << endl;
    cout << "Total steps: " << total_steps << endl;
    cout << "Average steps: " << (numbers_tested > 0 ? total_steps / numbers_tested : 0) << endl;
    cout << "Max steps: " << max_steps << " at "; print_uint128(max_steps_number); cout << endl;
    cout << "Cycles found: " << cycles_found << endl;
    cout << "Time: " << duration.count() << " ms" << endl;
    cout << "Throughput: " << (duration.count() > 0 ? (COUNT * 1000ULL) / duration.count() : 0) 
         << " numbers/sec (range-based)" << endl;
    
    return 0;
}
