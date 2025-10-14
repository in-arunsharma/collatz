// Advanced GPP Optimization Ideas - Potential 20-40% improvement

// 1. MEMORY PREFETCHING OPTIMIZATION
__attribute__((always_inline)) static inline
CollatzResult compute_collatz_prefetch(uint128_t n, const uint32_t* __restrict memo_ptr,
                                      uint64_t memo_limit, uint64_t max_steps = HOT_PATH_FUSE) {
    CollatzResult res = {0, n, false};
    if (n == 0 || n == 1) return res;
    
    uint128_t current = n;
    uint64_t steps = 0;
    
    while (current != 1) {
        if (steps >= max_steps) {
            res.steps = steps;
            res.overflow = true;
            return res;
        }
        
        if (current < memo_limit) {
            // PREFETCH next likely memory location
            uint128_t next_guess = (current & 1) ? (3 * current + 1) >> 1 : current >> 1;
            if (next_guess < memo_limit) {
                __builtin_prefetch(&memo_ptr[next_guess], 0, 1);  // Prefetch for read
            }
            
            uint32_t cached = memo_ptr[(uint64_t)current];
            if (cached != UNKNOWN) {
                res.steps = steps + cached;
                return res;
            }
        }
        
        if (res.peak < current) res.peak = current;
        
        // OPTIMIZED: Combine even/odd processing
        if ((current & 1) == 0) {
            int shift = __builtin_ctzll((uint64_t)current);  // Use 64-bit CTZ for speed
            if (shift == 0 && current > UINT64_MAX) {
                shift = __builtin_ctzll((uint64_t)(current >> 64)) + 64;
            }
            current >>= shift;
            steps += shift;
        } else {
            if (current > MAX_SAFE) {
                res.steps = steps;
                res.overflow = true;
                return res;
            }
            uint128_t t = 3 * current + 1;
            int k = __builtin_ctzll((uint64_t)t);  // Use 64-bit CTZ first
            if (k == 0 && t > UINT64_MAX) {
                k = __builtin_ctzll((uint64_t)(t >> 64)) + 64;
            }
            current = t >> k;
            steps += 1 + k;
        }
    }
    
    res.steps = steps;
    return res;
}

// 2. BATCH PROCESSING OPTIMIZATION
void process_batch_optimized(uint128_t* seeds, int batch_size, 
                            const uint32_t* __restrict memo_ptr,
                            ThreadStats& local) {
    // Process seeds in batches to improve cache locality
    for (int i = 0; i < batch_size; i += 8) {  // Process 8 at a time
        int end = std::min(i + 8, batch_size);
        
        // Prefetch all seeds in batch
        for (int j = i; j < end; j++) {
            if (seeds[j] < MEMO_TABLE_SIZE) {
                __builtin_prefetch(&memo_ptr[seeds[j]], 0, 1);
            }
        }
        
        // Process the batch
        for (int j = i; j < end; j++) {
            CollatzResult res = compute_collatz_prefetch(seeds[j], memo_ptr, MEMO_TABLE_SIZE);
            // Update local stats...
        }
    }
}

// 3. COMPILER OPTIMIZATION FLAGS
/*
Advanced compiler flags to try:
-O3 -march=sapphirerapids -mtune=sapphirerapids
-flto -ffast-math -funroll-loops
-fprefetch-loop-arrays
-fprofile-generate (for PGO)
-mavx512f -mavx512cd -mavx512bw -mavx512dq -mavx512vl
*/