#include <iostream>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <omp.h>

using namespace std;
using namespace std::chrono;

// Cache line padding to prevent false sharing
struct alignas(64) ThreadLocalData {
    uint64_t total_steps;
    uint64_t max_steps;
    __uint128_t max_steps_number;
    uint64_t cycles_found;
    uint64_t numbers_tested;
    char padding[64 - (sizeof(uint64_t) * 4 + sizeof(__uint128_t)) % 64];
};

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
    uint64_t low = (uint64_t)n;
    
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
inline int64_t collatz_steps(__uint128_t n, __uint128_t start_range) {
    __uint128_t current = n;
    uint64_t steps = 0;
    const uint64_t MAX_ITERATIONS = 100000;
    uint64_t iterations = 0;
    
    while (current != 1) {
        uint64_t n_steps = 0;
        current = collatz_algo(current, n_steps);
        steps += n_steps;
        
        if (current < start_range) {
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
    const int NUM_THREADS = 12;
    
    cout << "Collatz Multi-threaded V5 - Cache-Aligned" << endl;
    cout << "==========================================" << endl;
    cout << "Start: 2^71 = "; print_uint128(START); cout << endl;
    cout << "Range: " << COUNT << " numbers" << endl;
    cout << "Threads: " << NUM_THREADS << endl;
    cout << "Optimizations: Mod-6 + Cache-aligned thread-local data" << endl << endl;
    
    omp_set_num_threads(NUM_THREADS);
    
    // Cache-aligned per-thread data to prevent false sharing
    ThreadLocalData thread_data[NUM_THREADS];
    for (int i = 0; i < NUM_THREADS; i++) {
        thread_data[i].total_steps = 0;
        thread_data[i].max_steps = 0;
        thread_data[i].max_steps_number = 0;
        thread_data[i].cycles_found = 0;
        thread_data[i].numbers_tested = 0;
    }
    
    auto start_time = high_resolution_clock::now();
    
    // Find first odd number that is ≡ 1 or 5 (mod 6)
    __uint128_t first_n = START;
    if (!(first_n & 1)) first_n++;
    uint64_t mod6 = (uint64_t)(first_n % 6);
    if (mod6 == 3) first_n += 2;
    
    __uint128_t end = START + COUNT;
    
    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        ThreadLocalData& local = thread_data[tid];
        
        #pragma omp for schedule(static) nowait
        for (uint64_t i = 0; i < COUNT / 6; i++) {
            __uint128_t n = first_n + i * 6;
            
            if (n >= end) continue;
            
            // Process n ≡ 1 (mod 6)
            int64_t steps = collatz_steps(n, first_n);
            if (steps == -1) {
                local.cycles_found++;
            } else {
                local.total_steps += steps;
                if ((uint64_t)steps > local.max_steps) {
                    local.max_steps = steps;
                    local.max_steps_number = n;
                }
            }
            local.numbers_tested++;
            
            // Process n ≡ 5 (mod 6)
            __uint128_t n2 = n + 4;
            if (n2 < end) {
                steps = collatz_steps(n2, first_n);
                if (steps == -1) {
                    local.cycles_found++;
                } else {
                    local.total_steps += steps;
                    if ((uint64_t)steps > local.max_steps) {
                        local.max_steps = steps;
                        local.max_steps_number = n2;
                    }
                }
                local.numbers_tested++;
            }
        }
    }
    
    // Sequential reduction after parallel work (no lock contention)
    uint64_t total_steps = 0;
    uint64_t max_steps = 0;
    __uint128_t max_steps_number = 0;
    uint64_t cycles_found = 0;
    uint64_t numbers_tested = 0;
    
    for (int i = 0; i < NUM_THREADS; i++) {
        total_steps += thread_data[i].total_steps;
        numbers_tested += thread_data[i].numbers_tested;
        cycles_found += thread_data[i].cycles_found;
        if (thread_data[i].max_steps > max_steps) {
            max_steps = thread_data[i].max_steps;
            max_steps_number = thread_data[i].max_steps_number;
        }
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
    cout << "Speedup vs single-core: " << (3636.0 / duration.count()) << "x" << endl;
    
    return 0;
}
