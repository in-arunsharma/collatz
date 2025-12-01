#include <iostream>
#include <chrono>
#include <cstdint>
#include <cuda_runtime.h>

using namespace std;
using namespace std::chrono;

// CUDA device function - count trailing zeros for 128-bit
__device__ int ctz128(unsigned long long low, unsigned long long high) 
{
    if (low != 0) return __ffsll(low) - 1;  // __ffsll returns 1-based index
    if (high != 0) return 64 + __ffsll(high) - 1;
    return 128;
}

// CUDA kernel - each thread processes one number
__global__ void collatz_kernel(
    unsigned long long start_low,
    unsigned long long start_high,
    unsigned long long count,
    unsigned long long* step_counts,
    int* cycle_flags
) {
    unsigned long long idx = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned long long stride = blockDim.x * gridDim.x;
    
    // START is already aligned to mod-6 filtering
    // We need to find first odd n ≡ 1 or 5 (mod 6)
    unsigned long long first_low = start_low;
    unsigned long long first_high = start_high;
    
    // Make it odd
    if ((first_low & 1) == 0) {
        first_low++;
        if (first_low == 0) first_high++;
    }
    
    // Adjust to first n ≡ 1 or 5 (mod 6)
    // For 2^71, it's even, so +1 makes it 2^71+1 which is ≡ 1 (mod 6)
    unsigned long long mod6 = first_low % 6;
    if (mod6 == 3) {
        first_low += 2;
        if (first_low < 2) first_high++;
    }
    
    // Each thread processes multiple number pairs with stride
    for (unsigned long long i = idx; i < count / 6; i += stride) {
        // Calculate n = first_n + i * 6 (mod-6 filtering)
        unsigned long long n_low = first_low + i * 6;
        unsigned long long n_high = first_high;
        
        // Handle overflow
        if (n_low < first_low) n_high++;
        
        // Process n ≡ 1 (mod 6)
        unsigned long long curr_low = n_low;
        unsigned long long curr_high = n_high;
        unsigned long long steps = 0;
        unsigned long long iterations = 0;
        const unsigned long long MAX_ITER = 100000;
        
        while (!(curr_low == 1 && curr_high == 0)) {
            // Check if below start range (early termination)
            if (curr_high < start_high || (curr_high == start_high && curr_low < start_low)) {
                break;
            }
            
            if (curr_low & 1) {
                // 3n + 1
                unsigned long long temp_low = curr_low * 3;
                unsigned long long temp_high = curr_high * 3;
                
                // Add carry from low to high
                if (temp_low < curr_low) temp_high++;
                
                // Add 1
                temp_low++;
                if (temp_low == 0) temp_high++;
                
                // Count trailing zeros and shift
                int zeros = ctz128(temp_low, temp_high);
                
                if (zeros < 64) {
                    curr_low = (temp_low >> zeros) | (temp_high << (64 - zeros));
                    curr_high = temp_high >> zeros;
                } else {
                    curr_low = temp_high >> (zeros - 64);
                    curr_high = 0;
                }
                
                steps += zeros + 1;
            } else {
                // n / 2^k
                int zeros = ctz128(curr_low, curr_high);
                
                if (zeros < 64) {
                    curr_low = (curr_low >> zeros) | (curr_high << (64 - zeros));
                    curr_high = curr_high >> zeros;
                } else {
                    curr_low = curr_high >> (zeros - 64);
                    curr_high = 0;
                }
                
                steps += zeros;
            }
            
            if (++iterations > MAX_ITER) {
                cycle_flags[i * 2] = 1;
                steps = 0;
                break;
            }
        }
        
        step_counts[i * 2] = steps;
        
        // Process n ≡ 5 (mod 6) - n + 4
        n_low += 4;
        if (n_low < 4) n_high++;
        
        curr_low = n_low;
        curr_high = n_high;
        steps = 0;
        iterations = 0;
        
        while (!(curr_low == 1 && curr_high == 0)) {
            if (curr_high < start_high || (curr_high == start_high && curr_low < start_low)) {
                break;
            }
            
            if (curr_low & 1) {
                unsigned long long temp_low = curr_low * 3;
                unsigned long long temp_high = curr_high * 3;
                if (temp_low < curr_low) temp_high++;
                temp_low++;
                if (temp_low == 0) temp_high++;
                
                int zeros = ctz128(temp_low, temp_high);
                if (zeros < 64) {
                    curr_low = (temp_low >> zeros) | (temp_high << (64 - zeros));
                    curr_high = temp_high >> zeros;
                } else {
                    curr_low = temp_high >> (zeros - 64);
                    curr_high = 0;
                }
                steps += zeros + 1;
            } else {
                int zeros = ctz128(curr_low, curr_high);
                if (zeros < 64) {
                    curr_low = (curr_low >> zeros) | (curr_high << (64 - zeros));
                    curr_high = curr_high >> zeros;
                } else {
                    curr_low = curr_high >> (zeros - 64);
                    curr_high = 0;
                }
                steps += zeros;
            }
            
            if (++iterations > MAX_ITER) {
                cycle_flags[i * 2 + 1] = 1;
                steps = 0;
                break;
            }
        }
        
        step_counts[i * 2 + 1] = steps;
    }
}

void print_uint128(unsigned long long low, unsigned long long high) {
    if (high > 0) {
        cout << high << low;  // Simplified printing
    } else {
        cout << low;
    }
}

int main() {
    // 2^71 = (2^7 << 64) | 0 = high=128, low=0
    const unsigned long long START_LOW = 0;
    const unsigned long long START_HIGH = 128;  // 2^7, since 2^71 = 2^7 * 2^64
    const unsigned long long COUNT = 1000000000;
    
    cout << "Collatz GPU V1 - CUDA" << endl;
    cout << "=====================" << endl;
    cout << "Start: 2^71" << endl;
    cout << "Range: " << COUNT << " numbers" << endl;
    cout << "Optimizations: Mod-6 + CUDA parallel" << endl << endl;
    
    // Allocate device memory
    unsigned long long* d_steps;
    int* d_cycles;
    unsigned long long num_pairs = COUNT / 6 + 1;
    
    cudaMalloc(&d_steps, num_pairs * 2 * sizeof(unsigned long long));
    cudaMalloc(&d_cycles, num_pairs * 2 * sizeof(int));
    cudaMemset(d_cycles, 0, num_pairs * 2 * sizeof(int));
    
    // Launch configuration: 256 threads per block, enough blocks to cover all work
    int threadsPerBlock = 256;
    int numBlocks = (num_pairs + threadsPerBlock - 1) / threadsPerBlock;
    
    cout << "GPU Configuration:" << endl;
    cout << "Blocks: " << numBlocks << ", Threads/Block: " << threadsPerBlock << endl << endl;
    
    auto start_time = high_resolution_clock::now();
    
    // Launch kernel
    collatz_kernel<<<numBlocks, threadsPerBlock>>>(
        START_LOW, START_HIGH, COUNT, d_steps, d_cycles
    );
    
    // Wait for GPU to finish
    cudaDeviceSynchronize();
    
    auto end_time = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(end_time - start_time);
    
    // Copy results back
    unsigned long long* h_steps = new unsigned long long[num_pairs * 2];
    int* h_cycles = new int[num_pairs * 2];
    
    cudaMemcpy(h_steps, d_steps, num_pairs * 2 * sizeof(unsigned long long), cudaMemcpyDeviceToHost);
    cudaMemcpy(h_cycles, d_cycles, num_pairs * 2 * sizeof(int), cudaMemcpyDeviceToHost);
    
    // Calculate statistics
    unsigned long long total_steps = 0;
    unsigned long long max_steps = 0;
    unsigned long long cycles_found = 0;
    unsigned long long numbers_tested = 0;
    
    for (unsigned long long i = 0; i < num_pairs * 2; i++) {
        if (h_cycles[i]) {
            cycles_found++;
        } else {
            total_steps += h_steps[i];
            if (h_steps[i] > max_steps) max_steps = h_steps[i];
            numbers_tested++;
        }
    }
    
    cout << "Results:" << endl;
    cout << "--------" << endl;
    cout << "Numbers tested: " << numbers_tested << endl;
    cout << "Total steps: " << total_steps << endl;
    cout << "Average steps: " << (numbers_tested > 0 ? total_steps / numbers_tested : 0) << endl;
    cout << "Max steps: " << max_steps << endl;
    cout << "Cycles found: " << cycles_found << endl;
    cout << "Time: " << duration.count() << " ms" << endl;
    cout << "Throughput: " << (duration.count() > 0 ? (COUNT * 1000ULL) / duration.count() : 0) 
         << " numbers/sec (range-based)" << endl;
    cout << "Speedup vs single-core: " << (3636.0 / duration.count()) << "x" << endl;
    cout << "Speedup vs 12-core CPU: " << (469.0 / duration.count()) << "x" << endl;
    
    // Cleanup
    delete[] h_steps;
    delete[] h_cycles;
    cudaFree(d_steps);
    cudaFree(d_cycles);
    
    return 0;
}
