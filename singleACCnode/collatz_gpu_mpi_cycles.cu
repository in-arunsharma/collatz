// collatz_gpu_mpi_cycles.cu
// MPI + CUDA Collatz scanner with end-of-run cycle printing on rank 0.
//
// Build:
//   nvcc -O3 -arch=sm_90 --device-int128 collatz_gpu_mpi_cycles.cu -o collatz_gpu_mpi_cycles -lmpi
//
// Notes:
// - One MPI rank per GPU is the simplest setup (srun -n <num_gpus_total>).
// - The kernel scans candidates (skipping evens & multiples of 3) and pushes flagged seeds
//   (128-bit overflow or fuse-hit) to a device queue.
// - After the kernel, the host verifies flags with Brent(128), and if needed Brent(256) on CPU.
// - Only rank 0 prints a final summary; with --print-cycles it also prints each cycle line.
//
// CLI:
//   ./collatz_gpu_mpi_cycles <start_offset> <count> [--small-limit <bits>] [--tag <name>] [--no-print-cycles]
//
//   <start_offset>, <count> : range within [2^71 + start_offset, 2^71 + start_offset + count)
//   --small-limit <bits>    : memo table size as 2^bits (default: 20 → ~4MB)
//   --tag <name>            : tag (currently only used in headers; logging to files removed for clarity)
//   --no-print-cycles       : do not print the individual cycle lines (summary still prints)
//
// -------------------------------------------------------------------------

#include <mpi.h>
#include <cuda_runtime.h>

#include <cstdio>
#include <cstdint>
#include <cstring>
#include <vector>
#include <string>
#include <algorithm>
#include <iostream>
#include <cassert>

// ---------------- Common typedefs & constants ----------------
using u128 = unsigned __int128;

static constexpr uint64_t SAFETY_FUSE   = 100000;   // fast path fuse
static constexpr uint64_t EXTENDED_FUSE = 1000000;  // verification fuse
static constexpr uint32_t UNKNOWN       = UINT32_MAX;

// ---------------- Helpers: u128 (host & device) ----------------
__host__ __device__ static inline u128 u128_from_u64(uint64_t v){ return (u128)v; }

// Replace your current ctz_u128 with this pair:

// 64-bit trailing-zero count that works in both host & device builds
__host__ __device__ static inline int ctz64_u(uint64_t x) {
#ifdef __CUDA_ARCH__
    // device: __ffsll returns 1 + index of LSB set; subtract 1
    return __ffsll((long long)x) - 1;
#else
    // host: built-in works here
    return __builtin_ctzll(x);
#endif
}

__host__ __device__ static inline int ctz_u128(u128 x) {
    if (x == 0) return 128;
    uint64_t lo = (uint64_t)x;
    if (lo) return ctz64_u(lo);
    uint64_t hi = (uint64_t)(x >> 64);
    return 64 + ctz64_u(hi);
}
__host__ __device__ static inline u128 max_safe_u128() {
    // ((~(u128)0) - 1) / 3
    return ( ( (u128)~(u128)0 ) - 1 ) / 3;
}
// Convert u128 to decimal string (host only)
// Note: for cycle printing only; not performance critical
static std::string to_string_u128(u128 v) {
    if (v == 0) return "0";
    std::string s; s.reserve(40);
    while (v) {
        int d = (int)(v % 10);
        s.push_back('0' + d);
        v /= 10;
    }
    std::reverse(s.begin(), s.end());
    return s;
}


// ---------------- Device memo table ----------------
struct DevMemo {
    const uint32_t* tbl;
    uint64_t limit;
};

__device__ __forceinline__ uint32_t memo_lookup(const DevMemo& M, u128 x){
    // If x fits in 64 bits and < limit, check memo
    uint64_t v = (uint64_t)x;
    if ((u128)v == x && v < M.limit){
        uint32_t m = M.tbl[v];
        return m;
    }
    return UNKNOWN;
}

// ---------------- Device: fast-path compute (128-bit) ----------------
__device__ __forceinline__
void compute_collatz_readonly_gpu(u128 n, const DevMemo& M,
                                  uint64_t max_steps,
                                  uint64_t& steps_out,
                                  bool& overflow_out,
                                  u128& peak_out)
{
    steps_out = 0;
    overflow_out = false;
    peak_out = n;
    if (n <= 1) return;

    const u128 MAX_SAFE = max_safe_u128();
    u128 cur = n;

    while (cur != 1) {
        if (steps_out >= max_steps) return;

        if (cur < M.limit) {
            uint32_t cached = memo_lookup(M, cur);
            if (cached != UNKNOWN) {
                steps_out += cached;
                return;
            }
        }
        if (cur > peak_out) peak_out = cur;

        if ( ((uint64_t)cur & 1ULL) == 0ULL ){
            int sh = ctz_u128(cur);
            cur >>= sh;
            steps_out += sh;
        } else {
            if (cur > MAX_SAFE) { overflow_out = true; return; }
            u128 t = 3*cur + 1;
            int k = ctz_u128(t);
            cur = t >> k;
            steps_out += 1 + k;
        }
    }
}

// ---------------- Device: flagged queue ----------------
struct FlaggedSeed {
    u128     n;
    uint64_t steps_before;
    uint8_t  reason;  // 1=overflow, 2=fuse, 3=both (unlikely)
};

__device__ unsigned int g_q_count = 0;

__device__ void push_flagged(FlaggedSeed* q, unsigned int q_cap,
                             u128 n, uint64_t steps, uint8_t reason)
{
    unsigned int pos = atomicAdd(&g_q_count, 1);
    if (pos < q_cap){
        q[pos].n = n;
        q[pos].steps_before = steps;
        q[pos].reason = reason;
    }
    // else: drop (rare); could set a global overflow flag if desired
}

// ---------------- Device: kernel ----------------
struct ThreadAcc {
    uint64_t tested;
    uint64_t steps_sum;
    uint64_t overflow_cnt;
    uint64_t fuse_cnt;
};

__global__ void collatz_kernel(u128 start, u128 end, uint64_t first_delta,
                               DevMemo M,
                               FlaggedSeed* __restrict flag_q, unsigned int q_cap,
                               ThreadAcc* __restrict acc_out)
{
    const uint64_t tid = blockIdx.x * (uint64_t)blockDim.x + threadIdx.x;
    const uint64_t stride = gridDim.x * (uint64_t)blockDim.x;

    // Walk the 2/4 delta sequence to position this thread at its first candidate
    u128 n = start;
    uint64_t d = first_delta;
    for (uint64_t i=0;i<tid;i++){
        n += d;
        d ^= 6ULL; // 2 ↔ 4
    }

    // Accumulators
    uint64_t tested=0, steps_sum=0, over_c=0, fuse_c=0;

    // Grid-stride loop: for each "round", advance 'stride' positions in the 2/4 toggle
    for (;;) {
        if (n >= end) break;

        // compute
        uint64_t steps=0;
        bool overflow=false;
        u128 peak;
        compute_collatz_readonly_gpu(n, M, SAFETY_FUSE, steps, overflow, peak);

        tested++;
        steps_sum += steps;

        uint8_t reason = 0;
        if (overflow){ over_c++; reason |= 1; }
        if (steps >= SAFETY_FUSE){ fuse_c++; reason |= 2; }
        if (reason){
            push_flagged(flag_q, q_cap, n, steps, reason);
        }

        // advance by 'stride' toggles
        for (uint64_t s=0; s<stride; s++){
            n += d;
            d ^= 6ULL;
        }
    }

    // Write per-thread accumulators
    acc_out[tid].tested = tested;
    acc_out[tid].steps_sum = steps_sum;
    acc_out[tid].overflow_cnt = over_c;
    acc_out[tid].fuse_cnt = fuse_c;
}

// ---------------- Host: small memo table ----------------
// Minimal but correct base entries. You can graft a full precompute if desired.
static void build_min_memo(std::vector<uint32_t>& memo){
    if (memo.size() == 0) return;
    std::fill(memo.begin(), memo.end(), UNKNOWN);
    memo[0] = UNKNOWN;
    if (memo.size() > 1) memo[1] = 0;
    if (memo.size() > 2) memo[2] = 1;
    if (memo.size() > 3) memo[3] = 7;
}

// ---------------- Host: align start & delta (skip evens and multiples of 3) ----------------
static inline uint32_t mod3_u128(u128 n) {
    uint64_t lo = (uint64_t)n;
    uint64_t hi = (uint64_t)(n >> 64);
    uint32_t r_lo = lo % 3;
    uint32_t r_hi = hi % 3;
    return (r_lo + r_hi) % 3;
}

static inline void align_start_and_delta(u128 &n, uint64_t &delta) {
    if ( ((uint64_t)n & 1ULL) == 0ULL ) n += 1;
    uint32_t r3 = mod3_u128(n);
    if (r3 == 0) n += 2;
    r3 = mod3_u128(n);
    delta = (r3 == 1) ? 4 : 2;
}

// ---------------- Host: Brent(128) & Brent(256) ----------------
struct Brent128Res {
    bool cycle_found;
    uint64_t steps;
    uint64_t cycle_length;
    u128 meet_value;
    bool overflow;
};

static Brent128Res detect_cycle_brent_128_host(u128 n, uint64_t max_steps) {
    Brent128Res r{false,0,0,0,false};
    u128 tortoise = n, hare = n;
    uint64_t power = 1, lambda = 1, steps = 0;
    const u128 MAX_SAFE = max_safe_u128();

    while (steps < max_steps) {
        if (hare == 1) { r.steps = steps; return r; }

        if ( ((uint64_t)hare & 1ULL) == 0ULL ){
            int sh = ctz_u128(hare);
            hare >>= sh; steps += sh;
        } else {
            if (hare > MAX_SAFE){ r.overflow=true; r.steps=steps; return r; }
            u128 t = 3*hare + 1;
            int k = ctz_u128(t);
            hare = t >> k; steps += 1 + k;
        }
        if (tortoise == hare){
            r.cycle_found = true;
            r.cycle_length = lambda;
            r.steps = steps;
            r.meet_value = hare;
            return r;
        }
        if (lambda == power){
            tortoise = hare;
            power <<= 1;
            lambda = 0;
        }
        lambda++;
    }
    r.steps = steps; r.overflow = true;
    return r;
}

// Lightweight uint256 for host Brent(256)
struct u256 {
    uint64_t v[4]; // little-endian limbs: v[0] is lowest

    u256(){ v[0]=v[1]=v[2]=v[3]=0; }
    explicit u256(u128 x){
        v[0] = (uint64_t)x;
        v[1] = (uint64_t)(x >> 64);
        v[2] = v[3] = 0;
    }
    bool is_zero() const { return (v[0]|v[1]|v[2]|v[3])==0; }
    bool is_one()  const { return v[0]==1 && v[1]==0 && v[2]==0 && v[3]==0; }
    bool is_even() const { return (v[0] & 1ULL)==0; }
};

static inline int ctz_u256(const u256& a){
    if (a.v[0]) return __builtin_ctzll(a.v[0]);
    if (a.v[1]) return 64 + __builtin_ctzll(a.v[1]);
    if (a.v[2]) return 128 + __builtin_ctzll(a.v[2]);
    if (a.v[3]) return 192 + __builtin_ctzll(a.v[3]);
    return 256;
}
static inline bool eq_u256(const u256& a, const u256& b){
    return a.v[0]==b.v[0] && a.v[1]==b.v[1] && a.v[2]==b.v[2] && a.v[3]==b.v[3];
}
static inline bool gt_u256(const u256& a, const u256& b){
    for (int i=3;i>=0;i--){
        if (a.v[i] > b.v[i]) return true;
        if (a.v[i] < b.v[i]) return false;
    }
    return false;
}
static inline u256 shr_u256(const u256& a, int sh){
    if (sh<=0) return a;
    if (sh>=256) return u256();
    u256 r;
    int limb = sh/64, bits = sh%64;
    for (int i=0;i<4-limb;i++){
        uint64_t low = a.v[i+limb] >> bits;
        uint64_t hi  = 0;
        if (bits && i+limb+1<4) hi = a.v[i+limb+1] << (64-bits);
        r.v[i] = low | hi;
    }
    for (int i=4-limb;i<4;i++) r.v[i]=0;
    return r;
}
static inline u256 add_u256_u64(const u256& a, uint64_t b){
    u256 r=a;
    __uint128_t s = ( (__uint128_t)r.v[0] + b );
    r.v[0] = (uint64_t)s;
    uint64_t c = (uint64_t)(s>>64);
    for (int i=1;i<4 && c;i++){
        s = ( (__uint128_t)r.v[i] + c );
        r.v[i] = (uint64_t)s;
        c = (uint64_t)(s>>64);
    }
    return r;
}
static inline u256 mul_u256_u64(const u256& a, uint64_t m){
    u256 r;
    __uint128_t c=0;
    for (int i=0;i<4;i++){
        __uint128_t p = (__uint128_t)a.v[i]*m + c;
        r.v[i] = (uint64_t)p;
        c = p>>64;
    }
    return r;
}
static inline bool would_overflow_3n_plus_1_u256(const u256& a){
    // Keep same MAX_SAFE_256 as user code (~0x5555... pattern)
    static const u256 LIMIT = []{
        u256 t; t.v[0]=0x5555555555555555ULL;
                t.v[1]=0x5555555555555555ULL;
                t.v[2]=0x5555555555555555ULL;
                t.v[3]=0x5555555555555555ULL;
        return t;
    }();
    return gt_u256(a, LIMIT) || eq_u256(a, LIMIT);
}

struct Brent256Res {
    bool cycle_found;
    uint64_t steps;
    uint64_t cycle_length;
    std::string cycle_value_dec; // for printing/log; derived at detection
    bool overflow;
};

static std::string u256_to_dec(u256 x){
    if (x.is_zero()) return "0";
    // long division by 10
    std::string out;
    while (!x.is_zero()){
        __uint128_t rem=0;
        for (int i=3;i>=0;i--){
            __uint128_t cur = (rem<<64) | x.v[i];
            x.v[i] = (uint64_t)(cur / 10);
            rem    = (uint64_t)(cur % 10);
        }
        out.push_back('0' + (int)rem);
    }
    std::reverse(out.begin(), out.end());
    return out;
}

static Brent256Res detect_cycle_brent_256_host(u256 n, uint64_t max_steps) {
    Brent256Res r{false,0,0,"",false};
    u256 tortoise = n, hare = n;
    uint64_t power=1, lambda=1, steps=0;

    while (steps < max_steps) {
        if (hare.is_one()) { r.steps=steps; return r; }

        if (hare.is_even()){
            int sh = ctz_u256(hare);
            hare = shr_u256(hare, sh);
            steps += sh;
        } else {
            if (would_overflow_3n_plus_1_u256(hare)){ r.overflow=true; r.steps=steps; return r; }
            u256 t = add_u256_u64(mul_u256_u64(hare,3), 1);
            int k = ctz_u256(t);
            hare = shr_u256(t, k);
            steps += 1 + k;
        }

        if (eq_u256(tortoise, hare)){
            r.cycle_found = true;
            r.cycle_length = lambda;
            r.steps = steps;
            r.cycle_value_dec = u256_to_dec(hare);
            return r;
        }

        if (lambda == power){
            tortoise = hare;
            power <<= 1;
            lambda = 0;
        }
        lambda++;
    }
    r.steps=steps; r.overflow=true;
    return r;
}

// ---------------- Host: align, counting, and printing ----------------
static void align_task_range(u128 base, uint64_t start_offset, uint64_t count,
                             u128& start, u128& end, uint64_t& delta)
{
    start = base + (u128)start_offset;
    end   = start + (u128)count;
    align_start_and_delta(start, delta);
}

// Cycle record for final printing
struct CycleRec {
    u128 seed;
    uint64_t len;
    // for 128-bit cycles: meet value in u128; for 256-bit we'll store decimal string
    bool is256;
    u128 meet128;
    std::string meet256_dec;
};

// ---------------- Main ----------------
int main(int argc, char** argv){
    MPI_Init(&argc, &argv);
    int rank=0, size=1;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);

    if (argc < 3){
        if (rank==0){
            fprintf(stderr,"Usage: %s <start_offset> <count> [--small-limit <bits>] [--tag <name>] [--no-print-cycles]\n", argv[0]);
        }
        MPI_Finalize(); return 1;
    }

    uint64_t start_offset = std::stoull(argv[1]);
    uint64_t count        = std::stoull(argv[2]);
    uint32_t SMALL_LIMIT_BITS = 20;
    const char* run_tag   = "v15cuda";
    bool print_cycles = true;

    for (int i=3;i<argc;i++){
        if (!strcmp(argv[i],"--small-limit") && i+1<argc) SMALL_LIMIT_BITS = std::stoul(argv[++i]);
        else if (!strcmp(argv[i],"--tag") && i+1<argc)     run_tag = argv[++i];
        else if (!strcmp(argv[i],"--no-print-cycles"))     print_cycles = false;
    }

    // pick device per-rank
    int dev_count=1;
    cudaGetDeviceCount(&dev_count);
    cudaSetDevice(rank % std::max(1,dev_count));

    // Build minimal memo and upload
    const uint64_t small_limit = (1ULL << SMALL_LIMIT_BITS);
    std::vector<uint32_t> memo(small_limit);
    build_min_memo(memo);

    uint32_t* d_memo=nullptr;
    cudaMalloc(&d_memo, small_limit*sizeof(uint32_t));
    cudaMemcpy(d_memo, memo.data(), small_limit*sizeof(uint32_t), cudaMemcpyHostToDevice);
    DevMemo M{d_memo, small_limit};

    // MPI slice
    uint64_t count_per = count/size;
    uint64_t start_r   = start_offset + (uint64_t)rank * count_per;
    uint64_t count_r   = (rank==size-1) ? (count - (uint64_t)rank*count_per) : count_per;

    // Task range near 2^71
    const u128 base = (u128)1 << 71;
    u128 start_u128, end_u128;
    uint64_t delta=4;
    align_task_range(base, start_r, count_r, start_u128, end_u128, delta);

    if (rank==0){
        fprintf(stderr,"[MPI] size=%d\n", size);
    }

    // Launch params (tune as needed)
    const int BLOCK=256;
    const int GRID = 2048; // ~524,288 threads
    const size_t NTHREADS = (size_t)BLOCK * GRID;

    // Device accumulators
    ThreadAcc* d_acc=nullptr;
    cudaMalloc(&d_acc, NTHREADS*sizeof(ThreadAcc));

    // Device flag queue
    const unsigned int QCAP = 2'000'000; // capacity for flagged seeds per rank
    FlaggedSeed* d_q=nullptr;
    cudaMalloc(&d_q, QCAP*sizeof(FlaggedSeed));
    // reset queue counter
    cudaMemset(&g_q_count, 0, sizeof(unsigned int));

    // Barrier & timing
    MPI_Barrier(MPI_COMM_WORLD);
    double t0 = MPI_Wtime();

    // Launch kernel
    collatz_kernel<<<GRID,BLOCK>>>(start_u128, end_u128, delta, M, d_q, QCAP, d_acc);
    cudaDeviceSynchronize();

    // Pull back accumulators
    std::vector<ThreadAcc> acc(NTHREADS);
    cudaMemcpy(acc.data(), d_acc, NTHREADS*sizeof(ThreadAcc), cudaMemcpyDeviceToHost);

    // Pull back queue size and contents
    unsigned int h_q_count = 0;
    cudaMemcpyFromSymbol(&h_q_count, g_q_count, sizeof(unsigned int));
    if (h_q_count > QCAP) h_q_count = QCAP; // safety
    std::vector<FlaggedSeed> flagged(h_q_count);
    if (h_q_count) {
        cudaMemcpy(flagged.data(), d_q, h_q_count*sizeof(FlaggedSeed), cudaMemcpyDeviceToHost);
    }

    // Local aggregates
    unsigned long long tested_local=0, steps_local=0, over_local=0, fuse_local=0;
    for (const auto& a : acc){
        tested_local   += a.tested;
        steps_local    += a.steps_sum;
        over_local     += a.overflow_cnt;
        fuse_local     += a.fuse_cnt;
    }

    // Verification on host: Brent(128) and Brent(256) as needed
    std::vector<CycleRec> cycles_local;
    cycles_local.reserve(1024);

    for (const auto& fs : flagged){
        // First, Brent(128)
        Brent128Res r128 = detect_cycle_brent_128_host(fs.n, EXTENDED_FUSE);
        if (r128.cycle_found){
            CycleRec c;
            c.seed = fs.n;
            c.len  = r128.cycle_length;
            c.is256 = false;
            c.meet128 = r128.meet_value;
            cycles_local.push_back(c);
            continue;
        }
        if (r128.overflow){
            // Promote to 256
            u256 seed256(fs.n);
            Brent256Res r256 = detect_cycle_brent_256_host(seed256, EXTENDED_FUSE);
            if (r256.cycle_found){
                CycleRec c;
                c.seed = fs.n;
                c.len  = r256.cycle_length;
                c.is256 = true;
                c.meet256_dec = r256.cycle_value_dec;
                cycles_local.push_back(c);
            }
            // else: either overflowed 256 or timed out → no cycle confirmed
        }
        // else: verified OK within 128 without cycle (fuse false positive)
    }

    // Timing
    MPI_Barrier(MPI_COMM_WORLD);
    double elapsed = MPI_Wtime() - t0;

    // MPI reductions
    double elapsed_max=0.0;
    MPI_Reduce(&elapsed, &elapsed_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    unsigned long long tested_sum=0, steps_sum=0, over_sum=0, fuse_sum=0;
    MPI_Reduce(&tested_local, &tested_sum, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&steps_local,  &steps_sum,  1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&over_local,   &over_sum,   1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&fuse_local,   &fuse_sum,   1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

    // Gather cycles to rank 0 for printing once
    // We'll serialize as plain text lines to keep it simple.
    // First, build local string
    std::string local_cycles_blob;
    {
        for (const auto& c : cycles_local){
            local_cycles_blob += (c.is256 ? "256 " : "128 ");
            local_cycles_blob += to_string_u128(c.seed);
            local_cycles_blob += " ";
            local_cycles_blob += std::to_string(c.len);
            local_cycles_blob += " ";
            if (c.is256) {
                local_cycles_blob += c.meet256_dec;
            } else {
                local_cycles_blob += to_string_u128(c.meet128);
            }
            local_cycles_blob += "\n";
        }
    }
    // MPI gather counts then data
    int local_len = (int)local_cycles_blob.size();
    std::vector<int> all_lens;
    if (rank==0) all_lens.resize(size);
    MPI_Gather(&local_len, 1, MPI_INT, rank==0?all_lens.data():nullptr, 1, MPI_INT, 0, MPI_COMM_WORLD);

    std::string all_blob;
    if (rank==0){
        size_t total=0;
        for (int i=0;i<size;i++) total += all_lens[i];
        all_blob.resize(total);
    }

    // displacements
    std::vector<int> displs;
    if (rank==0){
        displs.resize(size);
        int off=0;
        for (int i=0;i<size;i++){ displs[i]=off; off += all_lens[i]; }
    }

    MPI_Gatherv(local_cycles_blob.data(), local_len, MPI_CHAR,
                rank==0?all_blob.data():nullptr,
                rank==0?all_lens.data():nullptr,
                rank==0?displs.data():nullptr,
                MPI_CHAR, 0, MPI_COMM_WORLD);

    // Rank 0 prints once
    if (rank==0){
        double thr = tested_sum / std::max(1e-9, elapsed_max);
        double avg_steps = tested_sum ? (double)steps_sum / (double)tested_sum : 0.0;

        std::cout << "\n=== V1.5-cuda (MPI+CUDA) Summary ===\n";
        std::cout << "Tag:               " << run_tag << "\n";
        std::cout << "MPI ranks:         " << size << "\n";
        std::cout << "Tested:            " << tested_sum << " numbers\n";
        std::cout << "Time (max rank):   " << (uint64_t)(elapsed_max*1000) << " ms\n";
        std::cout << "Throughput:        " << (uint64_t)thr << " nums/sec\n";
        std::cout << "Avg steps:         " << avg_steps << "\n";
        std::cout << "Overflows (fast):  " << over_sum << "\n";
        std::cout << "Fuse hits (fast):  " << fuse_sum << "\n";

        // Cycle print at the end (single place)
        if (!all_blob.empty()){
            // Count lines
            size_t cnt = std::count(all_blob.begin(), all_blob.end(), '\n');
            std::cout << "Cycles found:      " << cnt << "\n";
            if (print_cycles){
                std::cout << "--- CYCLES (format: <bits> <seed> <length> <meet_value>) ---\n";
                std::cout << all_blob;
            }
        } else {
            std::cout << "Cycles found:      0\n";
        }
        std::cout << "=== SUCCESS ===\n";
    }

    // Cleanup
    cudaFree(d_q);
    cudaFree(d_acc);
    cudaFree(d_memo);

    MPI_Finalize();
    return 0;
}
