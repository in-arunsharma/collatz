// V1.5-mpi-lean: MPI + OpenMP wrapper around the proven V1.5 hot loop
// - Keep the fast 128-bit hot path and read-only memo table
// - No cold-queue work in the timed section (pure hot run for max throughput)
// - One MPI rank per node recommended; each rank uses all cores via OpenMP
//
// Build (GCC):   mpicxx -O3 -march=native -flto -fopenmp -DNDEBUG -o V1.5-mpi-lean V15_MPI_Lean.cpp
// Build (ICX):   mpiicx -O3 -xHost -qopenmp -flto -DNDEBUG -o V1.5-mpi-lean V15_MPI_Lean.cpp
//
// Run (2 nodes, 1 rank/node; SLURM): srun --ntasks-per-node=1 --cpus-per-task=112 --cpu-bind=cores ./V1.5-mpi-lean 0 1000000000
// Env (recommended):
//   OMP_PLACES=cores  OMP_PROC_BIND=spread  I_MPI_PIN=1  I_MPI_PIN_DOMAIN=omp  KMP_BLOCKTIME=0
//
// Usage: ./V1.5-mpi-lean <start_offset> <count> [--threads N] [--small-limit BITS]
//
// Notes:
// - Uses mod-6 filter (only n ≡ 1,5 mod 6).
// - Memo table default 2^20 (4 MB), read-only during compute.
// - Prints aggregate throughput and per-rank rate on rank 0.

#include <mpi.h>
#include <omp.h>
#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <string>
#include <iostream>
#include <iomanip>
#include <ctime>

// ----- Basic 128-bit -----
typedef __uint128_t uint128_t;

// ----- Globals (read-only after precompute) -----
static uint32_t SMALL_LIMIT_BITS = 20;            // 2^20 = 1,048,576 entries
static uint64_t small_limit = 0;
static std::vector<uint32_t> memo;
static constexpr uint32_t UNKNOWN = UINT32_MAX;
static constexpr uint64_t HOT_FUSE = 100000;
static constexpr uint128_t MAX_SAFE_128 = ((~(uint128_t)0) - 1) / 3;

// ----- 128-bit helpers -----
static inline int ctz_u128(uint128_t x) {
    if (x == 0) return 128;
    uint64_t lo = (uint64_t)x;
    if (lo) return __builtin_ctzll(lo);
    uint64_t hi = (uint64_t)(x >> 64);
    return 64 + __builtin_ctzll(hi);
}
static inline uint128_t u128_from_u64(uint64_t x){ return (uint128_t)x; }

static inline void print_u128(uint128_t v) {
    if (v == 0) { std::cout << '0'; return; }
    char buf[64]; int p=0;
    while (v){ buf[p++] = '0' + (v % 10); v/=10; }
    while (p--) std::cout << buf[p];
}
static inline std::string to_string_u128(uint128_t v){
    if (!v) return "0";
    std::string s; s.reserve(40);
    while (v){ s.push_back('0' + (v%10)); v/=10; }
    std::reverse(s.begin(), s.end());
    return s;
}

// ----- Memo precompute (single-threaded, small & fast) -----
static void precompute_small_table() {
    memo.assign(small_limit, UNKNOWN);
    memo[0] = UNKNOWN; memo[1] = 0;

    for (uint64_t n = 2; n < small_limit; ++n) {
        if (memo[n] != UNKNOWN) continue;

        std::vector<uint64_t> path_val;
        std::vector<uint64_t> path_steps;
        uint128_t cur = (uint128_t)n;
        uint64_t steps = 0;
        bool completed = false;

        while (true) {
            if (cur == 1) { completed = true; break; }
            if (cur < small_limit && memo[(uint64_t)cur] != UNKNOWN) {
                steps += memo[(uint64_t)cur];
                completed = true; break;
            }
            if (steps >= HOT_FUSE) break;

            if (cur < small_limit) {
                path_val.push_back((uint64_t)cur);
                path_steps.push_back(steps);
            }

            if (((uint64_t)cur & 1ULL) == 0) {
                int sh = ctz_u128(cur);
                cur >>= sh; steps += sh;
            } else {
                if (cur > MAX_SAFE_128) break;
                uint128_t t = 3*cur + 1;
                int k = ctz_u128(t);
                cur = t >> k; steps += 1 + k;
            }
        }
        if (!completed) continue;
        for (int i=(int)path_val.size()-1; i>=0; --i) {
            uint64_t v = path_val[i];
            if (memo[v]==UNKNOWN) memo[v] = (uint32_t)(steps - path_steps[i]);
        }
    }
}

static bool validate_memo_table() {
    auto has=[&](uint64_t i){ return i < small_limit; };
    if (has(1) && memo[1]!=0) { std::fprintf(stderr,"[VALIDATE] memo[1]=%u\n",memo[1]); return false; }
    if (has(2) && memo[2]!=1) { std::fprintf(stderr,"[VALIDATE] memo[2]=%u\n",memo[2]); return false; }
    if (has(3) && memo[3]!=7) { std::fprintf(stderr,"[VALIDATE] memo[3]=%u\n",memo[3]); return false; }
    return true;
}

// ----- Hot kernel (read-only memo; no logging/queues) -----
struct CollatzResult { uint64_t steps; uint128_t peak; bool overflow; };

__attribute__((always_inline)) static inline
CollatzResult compute_collatz_readonly(uint128_t n,
                                       const uint32_t* __restrict memo_ptr,
                                       uint64_t memo_lim,
                                       uint64_t fuse = HOT_FUSE)
{
    CollatzResult r{0,n,false};
    if (n<=1) return r;
    uint128_t cur = n; uint64_t steps=0;
    while (cur != 1) {
        if (steps >= fuse) { r.steps=steps; return r; }
        if (cur < memo_lim) {
            uint32_t m = memo_ptr[(uint64_t)cur];
            if (m != UNKNOWN) { r.steps = steps + m; return r; }
        }
        if (cur > r.peak) r.peak = cur;

        if (((uint64_t)cur & 1ULL)==0) {
            int sh = ctz_u128(cur);
            cur >>= sh; steps += sh;
        } else {
            if (cur > MAX_SAFE_128) { r.overflow=true; r.steps=steps; return r; }
            uint128_t t = 3*cur + 1;
            int k = ctz_u128(t);
            cur = t >> k; steps += 1 + k;
        }
    }
    r.steps = steps; return r;
}

// ----- mod-6 helpers -----
static inline uint32_t mod3_u128(uint128_t n) {
    uint64_t lo = (uint64_t)n, hi = (uint64_t)(n>>64);
    return (lo%3 + hi%3) % 3; // because 2^64 ≡ 1 (mod 3)
}

static inline void align_start_and_delta(uint128_t& n, uint64_t& delta) {
    if (((uint64_t)n & 1ULL) == 0) n += 1;
    uint32_t r3 = mod3_u128(n);
    if (r3 == 0) n += 2;
    delta = (mod3_u128(n) == 1) ? 4 : 2; // start stride (then toggle 2<->4)
}

static inline uint128_t seed_from_index(uint128_t start_aligned, uint64_t idx) {
    // Deterministic mapping for n ≡ 1,5 (mod 6): pairs {+0,+4}, {+6,+10}, ...
    // index -> start + 6*(idx/2) + (idx%2 ? 4 : 0)
    return start_aligned + (uint128_t)( (idx >> 1) * 6ULL ) + ( (idx & 1ULL) ? 4ULL : 0ULL );
}

// Conservative candidate count (we still guard with n<end inside the loop)
static inline uint64_t approx_candidates(uint128_t start, uint128_t end) {
    if (end <= start) return 0;
    uint128_t range = end - start;
    return (uint64_t)((range * 2) / 3); // about 2/3 are 1 or 5 mod 6
}

// ----- rank-local compute -----
struct Stats {
    uint64_t tested=0, total_steps=0, max_steps=0;
    uint128_t hardest=0, max_peak=0;
};

static void run_rank(uint128_t global_base, uint64_t start_off, uint64_t count,
                     Stats& out)
{
    // Build rank's [start, end) range
    uint128_t start = global_base + start_off;
    uint128_t end   = start + count;

    // Align start to {odd, not %3==0} and pick first stride (2 or 4)
    uint64_t first_delta = 4;
    align_start_and_delta(start, first_delta);

    // Estimate candidate count; actual loop guards n<end
    uint64_t candidates = approx_candidates(start, end);

    const uint32_t* __restrict memo_ptr = memo.data();

    // Time the hot loop only
    double t0 = MPI_Wtime();

    // OpenMP parallel over abstract indices; compute n on the fly
    // Use static scheduling for very regular work; chunk size big to cut overhead.
    #pragma omp parallel
    {
        uint64_t local_tested=0, local_steps=0, local_max_steps=0;
        uint128_t local_hardest=0, local_max_peak=0;

        #pragma omp for schedule(static, 4096)
        for (uint64_t idx=0; idx<candidates; ++idx) {
            uint128_t n = seed_from_index(start, idx);
            if (n >= end) continue;

            CollatzResult r = compute_collatz_readonly(n, memo_ptr, small_limit, HOT_FUSE);
            local_tested++;
            local_steps += r.steps;
            if (r.steps > local_max_steps) { local_max_steps = r.steps; local_hardest = n; }
            if (r.peak  > local_max_peak)   local_max_peak   = r.peak;
        }

        #pragma omp critical
        {
            out.tested      += local_tested;
            out.total_steps += local_steps;
            if (local_max_steps > out.max_steps) { out.max_steps = local_max_steps; out.hardest = local_hardest; }
            if (local_max_peak  > out.max_peak)  { out.max_peak  = local_max_peak; }
        }
    }

    double t1 = MPI_Wtime();
    // Pack elapsed in total_steps upper bits? No: we’ll reduce with MAX(separate time)
    // We just print per-rank on rank 0 using reductions below.
    (void)t0; (void)t1;
}

// ----- main -----
int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int rank=0, size=1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (argc < 3) {
        if (rank==0) {
            std::fprintf(stderr,"Usage: %s <start_offset> <count> [--threads N] [--small-limit BITS]\n", argv[0]);
        }
        MPI_Finalize(); return 1;
    }

    // Parse CLI
    uint64_t global_start_offset = std::strtoull(argv[1], nullptr, 10);
    uint64_t global_count        = std::strtoull(argv[2], nullptr, 10);
    int user_threads = 0;
    for (int i=3;i<argc;i++){
        if (!std::strcmp(argv[i],"--threads") && i+1<argc) user_threads = std::atoi(argv[++i]);
        else if (!std::strcmp(argv[i],"--small-limit") && i+1<argc) SMALL_LIMIT_BITS = (uint32_t)std::atoi(argv[++i]);
    }

    // Threads
    if (user_threads>0) omp_set_num_threads(user_threads);
    int omp_threads = omp_get_max_threads();

    // Memo
    small_limit = (1ULL << SMALL_LIMIT_BITS);
    if (rank==0) std::fprintf(stderr, "[CONFIG] ranks=%d  threads/rank=%d  memo=2^%u (%llu entries, %.2f MB)\n",
                               size, omp_threads, SMALL_LIMIT_BITS,
                               (unsigned long long)small_limit, small_limit*4.0/(1024*1024));
    memo.resize(small_limit, UNKNOWN);
    precompute_small_table();
    if (rank==0 && !validate_memo_table()) { std::fprintf(stderr,"[ERROR] memo self-test failed\n"); MPI_Abort(MPI_COMM_WORLD,2); }

    // Split work evenly by count
    uint64_t count_per_rank = global_count / size;
    uint64_t my_start_off = global_start_offset + (uint64_t)rank * count_per_rank;
    uint64_t my_count     = (rank == size-1) ? (global_count - (uint64_t)rank * count_per_rank)
                                             : count_per_rank;

    // Base number: 2^71
    uint128_t base = u128_from_u64(1) << 71;

    // Barrier, then time
    MPI_Barrier(MPI_COMM_WORLD);
    double t0 = MPI_Wtime();

    Stats s{};
    run_rank(base, my_start_off, my_count, s);

    double t1 = MPI_Wtime();
    double local_time = t1 - t0;

    // Reduce
    uint64_t R_tested=0, R_steps=0, R_max_steps=0;
    double   R_max_time=0.0;

    MPI_Reduce(&s.tested,    &R_tested,    1, MPI_UINT64_T, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&s.total_steps,&R_steps,     1, MPI_UINT64_T, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&s.max_steps, &R_max_steps, 1, MPI_UINT64_T, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&local_time,  &R_max_time,  1, MPI_DOUBLE,   MPI_MAX, 0, MPI_COMM_WORLD);

    // (Optional) per-rank quick print for debugging
    // if (rank==0) std::fprintf(stderr,"Per-rank time (max): %.3f s\n", R_max_time);

    if (rank==0) {
        double thr = (R_max_time>0) ? ( (double)R_tested / R_max_time ) : 0.0;
        double avg_steps = (R_tested>0) ? ((double)R_steps / (double)R_tested) : 0.0;

        std::cout << "\n=== V1.5-mpi (lean) Results ===\n";
        std::cout << "MPI ranks:    " << size << "\n";
        std::cout << "Threads/rank: " << omp_threads << "\n";
        std::cout << "Tested:       " << R_tested << "\n";
        std::cout << "Wall time:    " << (uint64_t)(R_max_time*1000) << " ms\n";
        std::cout << "Throughput:   " << (uint64_t)thr << " nums/sec\n";
        std::cout << "Per-rank:     " << (uint64_t)(thr/size) << " nums/sec\n";
        std::cout << "Avg steps:    " << std::fixed << std::setprecision(2) << avg_steps << "\n";
        std::cout << "Max steps:    " << R_max_steps << "\n";
        std::cout << "\nHint: If Per-rank << 1.3e8, recheck compiler (GCC/ICX) and CPU pinning.\n";
    }

    MPI_Finalize();
    return 0;
}
