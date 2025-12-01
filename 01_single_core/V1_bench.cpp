#include <iostream>
#include <set>
#include <chrono>
#include <cstdint>

using namespace std;
using namespace std::chrono;

// V1: Uses std::set for cycle detection (slow but simple)
int64_t collatz_steps(__uint128_t n) {
    set<__uint128_t> seen;
    int64_t steps = 0;
    
    while (n != 1) {
        if (seen.find(n) != seen.end()) {
            return -1;  // Cycle detected
        }
        seen.insert(n);
        
        if (n % 2 == 0) {
            n = n / 2;
        } else {
            n = 3 * n + 1;
        }
        steps++;
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
    
    cout << "Collatz V1 - std::set Cycle Detection (Baseline)" << endl;
    cout << "================================================" << endl;
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
