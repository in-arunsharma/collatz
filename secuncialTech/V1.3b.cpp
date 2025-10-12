/*
 * Collatz Conjecture - Version 1.3b (4-Lane Interleaved Walker)
 * 
 * Description: V1.3a + 4-lane ILP optimization (still single-thread)
 * Starting from: 2^71 and pushing tested limits
 * 
 * Algorithm optimizations (same as V1.3a):
 * - All V1.3a optimizations (path compression, prefaulting, etc.)
 * - Early-exit memo with path-compression fill
 * 
 * ILP improvement (NEW in V1.3b):
 * - Process 4 independent seeds in lockstep (interleaved execution)
 * - Each lane: independent state (n, steps, max_excursion)
 * - Reduces CPU dependency stalls by providing parallel work
 * - When lane finishes, pull next seed and continue
 * - Still single-thread, no synchronization needed
 * 
 * Expected: 5-20% single-thread boost from better ILP utilization
 * Rationale: Modern CPUs can execute independent instructions in parallel
 */

#include <iostream>
#include <cstdint>
#include <chrono>
#include <iomanip>
#include <string>
#include <cstdlib>
#include <cstring>
#include <algorithm>
#include <array>

// Type for handling large numbers (at least 2^71)
typedef __uint128_t uint128_t;

// Maximum value for uint128_t
constexpr uint128_t U128_MAX = ~(uint128_t)0;

// Early-exit memo table with uint32_t and sentinel
static const uint32_t SMALL_LIM_DEFAULT_POW = 20;  // 2^20 = 4MB (cache-friendly)
static uint32_t SMALL_LIM = 1u << SMALL_LIM_DEFAULT_POW;
constexpr uint32_t UNKNOWN = 0xFFFFFFFFu;  // Sentinel for unknown entries
static uint32_t* steps_small = nullptr;  // steps-to-1; UNKNOWN if not set

// Hoist overflow check constant
static const uint128_t MAX_SAFE_N_FOR_3N1 = (U128_MAX - 1) / 3;

// Per-lane state (avoids false sharing when we go multi-threaded later)
struct LaneState {
    uint128_t n;
    uint64_t steps;
    uint128_t max_excursion;
    bool finished;
    bool used_memo;
};

// Structure to hold final results
struct CollatzResult {
    uint128_t number;
    uint64_t steps;
    bool reached_one;
    bool max_steps_exceeded;
    bool overflow_detected;
    uint128_t max_excursion;
    bool used_memo;
};

// Initialize memo table
void init_small_table() {
    steps_small = (uint32_t*)malloc(SMALL_LIM * sizeof(uint32_t));
    if (!steps_small) {
        std::cerr << "ERROR: Failed to allocate memo table (" 
                  << (SMALL_LIM * sizeof(uint32_t) / (1024*1024)) << " MB)" << std::endl;
        exit(1);
    }
    std::fill(steps_small, steps_small + SMALL_LIM, UNKNOWN);
    steps_small[1] = 0;  // Base case
}

// Prefault table pages to avoid page-fault storms during timing
void prefault_small_table() {
    const size_t page = 4096;
    volatile uint8_t sink = 0;
    for (size_t i = 0; i < (size_t)SMALL_LIM * sizeof(uint32_t); i += page) {
        sink ^= *((uint8_t*)steps_small + i);
    }
    (void)sink;  // Suppress unused warning
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

// Forward declaration
static CollatzResult compute_collatz_full(uint128_t n, uint64_t max_steps);

// Advance one lane by one "bundle" (either collapse evens OR do 3n+1+collapse)
// Returns: true if lane should continue, false if finished/memo-hit
static inline bool advance_lane_one_step(LaneState& lane, uint64_t max_steps) {
    if (lane.n == 1) {
        lane.finished = true;
        return false;
    }
    
    // Check if we've entered small domain and can use memo
    if (lane.n < SMALL_LIM) {
        uint32_t idx = (uint32_t)lane.n;
        uint32_t memo = steps_small[idx];
        
        if (memo != UNKNOWN) {
            // Memo hit - add cached steps and finish
            lane.steps += (uint64_t)memo;
            lane.finished = true;
            lane.used_memo = true;
            return false;
        }
        
        // Memo miss - we'll compute and cache (path compression happens in full computation)
        // For simplicity in interleaved mode, compute to completion when hitting small domain
        CollatzResult full = compute_collatz_full(lane.n, max_steps);
        if (full.overflow_detected || full.max_steps_exceeded) {
            lane.finished = true;
            return false;
        }
        lane.steps += full.steps;
        if (full.max_excursion > lane.max_excursion) {
            lane.max_excursion = full.max_excursion;
        }
        lane.finished = true;
        lane.used_memo = full.used_memo;
        return false;
    }
    
    // n >= SMALL_LIM: do one bundle
    int z = ctz_u128(lane.n);
    if (z) {
        // Collapse even run
        lane.n >>= z;
        lane.steps += z;
        if (lane.n > lane.max_excursion) lane.max_excursion = lane.n;
    } else {
        // n is odd: do 3n+1 and collapse result
        if (lane.n > MAX_SAFE_N_FOR_3N1) {
            lane.finished = true;
            return false;
        }
        
        uint128_t m = 3 * lane.n + 1;
        int k = ctz_u128(m);
        lane.n = m >> k;
        lane.steps += 1 + k;
        if (lane.n > lane.max_excursion) lane.max_excursion = lane.n;
    }
    
    // Check safety fuse
    if (lane.steps > max_steps) {
        lane.finished = true;
        return false;
    }
    
    return true;  // Continue
}

// Full computation for when we hit small domain (with path compression)
static CollatzResult compute_collatz_full(uint128_t n, uint64_t max_steps __attribute__((unused))) {
    CollatzResult r;
    r.number = n;
    r.steps = 0;
    r.reached_one = false;
    r.max_steps_exceeded = false;
    r.overflow_detected = false;
    r.max_excursion = n;
    r.used_memo = false;
    
    // This is the same logic as V1.3a for small-domain handling
    if (n < SMALL_LIM) {
        uint32_t idx = (uint32_t)n;
        uint32_t memo = steps_small[idx];
        
        if (memo != UNKNOWN) {
            r.steps = (uint64_t)memo;
            r.reached_one = true;
            r.used_memo = true;
            return r;
        }
        
        // Compute with path compression (simplified for now)
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
        while (m != 1 && m < SMALL_LIM) {
            uint32_t sm = (uint32_t)m;
            if (steps_small[sm] != UNKNOWN) {
                add += (uint64_t)steps_small[sm];
                break;
            }
            
            int z = ctz_u128(m);
            if (z) {
                m >>= z;
                add += z;
                continue;
            }
            
            uint128_t t = 3*m + 1;
            int k = ctz_u128(t);
            m = t >> k;
            add += 1 + k;
        }
        
        // Cache result
        if (add <= UINT32_MAX && steps_small[idx] == UNKNOWN) {
            steps_small[idx] = (uint32_t)add;
        }
        
        r.steps = add;
        r.reached_one = true;
        r.used_memo = true;
        return r;
    }
    
    // Should not reach here in interleaved mode
    r.reached_one = false;
    return r;
}

// Align starting value to first valid seed (odd & not multiple of 3)
static inline void align_start_and_delta(uint128_t &n, uint64_t &delta) {
    if ((uint64_t)n % 2 == 0) { n += 1; }
    
    uint32_t r3 = mod3_u128(n);
    if (r3 == 0) { n += 2; r3 = 2; }
    
    delta = (r3 == 1) ? 4 : 2;
}

// Get next valid seed with mod-6 filtering
static inline uint128_t next_seed(uint128_t n, uint64_t& delta) {
    uint128_t next = n + delta;
    delta ^= 6;  // Toggle 4 <-> 2
    return next;
}

int main(int argc, char* argv[]) {
    std::cout << "=== Collatz Conjecture - Version 1.3b (4-Lane ILP) ===" << std::endl;
    std::cout << "Interleaved walker for ILP boost (still single-thread)" << std::endl;
    std::cout << std::endl;
    
    // Starting point: 2^71
    uint128_t start = (uint128_t)1 << 71;
    uint128_t end = start + 1000000;
    uint64_t max_steps = 100000000;
    uint32_t small_limit_pow = SMALL_LIM_DEFAULT_POW;
    
    // Parse command line arguments
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
    
    SMALL_LIM = 1u << small_limit_pow;
    
    std::cout << "Memo table: 2^" << small_limit_pow << " = " << SMALL_LIM 
              << " entries (" << (SMALL_LIM * sizeof(uint32_t) / (1024*1024)) << " MB)" << std::endl;
    std::cout << "Lanes: 4 (interleaved execution for ILP)" << std::endl;
    std::cout << "Testing range from 2^71 + offset" << std::endl;
    std::cout << "Start (requested): ";
    print_uint128(start);
    std::cout << std::endl;
    
    // Align start to valid seed
    uint64_t delta = 4;
    align_start_and_delta(start, delta);
    
    std::cout << "Start (aligned): ";
    print_uint128(start);
    std::cout << std::endl;
    std::cout << "End: ";
    print_uint128(end);
    std::cout << std::endl;
    std::cout << "Filter: Testing only n ≡ 1,5 (mod 6)" << std::endl;
    std::cout << std::endl;
    
    // Initialize and prefault memo table
    std::cout << "Initializing memo table..." << std::flush;
    init_small_table();
    std::cout << " done." << std::endl;
    std::cout << "Prefaulting pages..." << std::flush;
    prefault_small_table();
    std::cout << " done." << std::endl;
    std::cout << std::endl;
    
    // Timing starts AFTER prefaulting
    auto start_time = std::chrono::high_resolution_clock::now();
    
    // Global aggregates
    uint64_t total_tested = 0;
    uint64_t longest_path_steps = 0;
    uint128_t number_with_max_steps = 0;
    uint64_t total_steps = 0;
    uint128_t global_max_excursion = 0;
    uint128_t number_with_max_excursion = 0;
    uint64_t memo_hits = 0;
    
    // Power-of-2 progress reporting
    const uint64_t progress_every = 16384;
    const uint64_t progress_mask = progress_every - 1;
    
    // 4-lane state
    constexpr int NUM_LANES = 4;
    std::array<LaneState, NUM_LANES> lanes;
    std::array<uint128_t, NUM_LANES> lane_start_numbers;
    
    // Initialize lanes with first 4 seeds
    uint128_t current_seed = start;
    for (int i = 0; i < NUM_LANES; ++i) {
        if (current_seed < end) {
            lane_start_numbers[i] = current_seed;
            lanes[i].n = current_seed;
            lanes[i].steps = 0;
            lanes[i].max_excursion = current_seed;
            lanes[i].finished = false;
            lanes[i].used_memo = false;
            current_seed = next_seed(current_seed, delta);
        } else {
            lanes[i].finished = true;
        }
    }
    
    // Main interleaved loop
    while (true) {
        bool any_active = false;
        
        // Advance each lane by one step (interleaved execution)
        for (int i = 0; i < NUM_LANES; ++i) {
            if (!lanes[i].finished) {
                any_active = true;
                advance_lane_one_step(lanes[i], max_steps);
            }
        }
        
        // Check for finished lanes and emit results
        for (int i = 0; i < NUM_LANES; ++i) {
            if (lanes[i].finished && lane_start_numbers[i] != 0) {
                // Emit this lane's result
                total_tested++;
                total_steps += lanes[i].steps;
                if (lanes[i].used_memo) memo_hits++;
                
                if (lanes[i].steps > longest_path_steps) {
                    longest_path_steps = lanes[i].steps;
                    number_with_max_steps = lane_start_numbers[i];
                }
                
                if (lanes[i].max_excursion > global_max_excursion) {
                    global_max_excursion = lanes[i].max_excursion;
                    number_with_max_excursion = lane_start_numbers[i];
                }
                
                // Progress report
                if ((total_tested & progress_mask) == 0) {
                    std::cout << "Progress: " << total_tested << " numbers tested" << std::endl;
                }
                
                // Load next seed into this lane
                if (current_seed < end) {
                    lane_start_numbers[i] = current_seed;
                    lanes[i].n = current_seed;
                    lanes[i].steps = 0;
                    lanes[i].max_excursion = current_seed;
                    lanes[i].finished = false;
                    lanes[i].used_memo = false;
                    current_seed = next_seed(current_seed, delta);
                } else {
                    // No more seeds - mark lane as permanently done
                    lane_start_numbers[i] = 0;
                }
            }
        }
        
        if (!any_active) break;
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
