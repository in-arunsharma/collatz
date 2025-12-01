#include <iostream>
#include <chrono>
#include <cstdint>

using namespace std;
using namespace std::chrono;

inline int ctz128(__uint128_t n) {
    uint64_t low = (uint64_t)n;
    if (low != 0) return __builtin_ctzll(low);
    return 64 + __builtin_ctzll((uint64_t)(n >> 64));
}

inline __uint128_t collatz_algo(__uint128_t n, uint64_t& steps) {
    __uint128_t result;
    if (n & 1) {
		result = 3 * n + 1;
        int zeros = ctz128(result);
        result >>= zeros;
        steps = zeros + 1;
        return result;
        
    } else {
		int zeros = ctz128(n);
		steps = zeros;
        return n >>= zeros;
    }
}

int64_t collatz_steps(__uint128_t n) {
    __uint128_t slow = n;
    __uint128_t fast = n;
    uint64_t steps = 0;
    
    
    while (slow != 1) 
    {
		uint64_t n_steps = 0;
        slow = collatz_algo(slow, n_steps);
        steps += n_steps;
        
        if (fast != 1) fast = collatz_algo(fast, steps);
        if (fast != 1) fast = collatz_algo(fast, steps);
        
        if (slow == fast && slow != 1) {
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
    const uint64_t COUNT = 10000000;
    
    cout << "Collatz V2 - Floyd Cycle Detection + Bit Ops" << endl;
    cout << "============================================" << endl;
    cout << "Start: 2^71 = "; print_uint128(START); cout << endl;
    cout << "Testing: " << COUNT << " odd numbers" << endl << endl;
    
    uint64_t total_steps = 0;
    uint64_t max_steps = 0;
    __uint128_t max_steps_number = 0;
    uint64_t cycles_found = 0;
    
    auto start_time = high_resolution_clock::now();
    
    __uint128_t n = START;
    if (!(n & 1)) n++;
    
    for (uint64_t i = 0; i < COUNT/2; i++) {
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
        
        n += 2;
    }
    
    auto end_time = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(end_time - start_time);
    
    cout << "=== RESULTS ===" << endl;
    cout << "Time: " << duration.count() << " ms" << endl;
    cout << "Numbers/sec: " << (COUNT * 1000ULL / duration.count()) << endl;
    cout << "Total steps: " << total_steps << endl;
    cout << "Cycles found: " << cycles_found << endl;
    cout << "Max steps: " << max_steps << " for number: "; 
    print_uint128(max_steps_number); 
    cout << endl;
    
    return 0;
}
