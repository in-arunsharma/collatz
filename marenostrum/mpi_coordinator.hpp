// mpi_coordinator.hpp - MPI Work Distribution for Collatz on MareNostrum 5
// Handles embarrassingly parallel work partitioning across nodes

#ifndef MPI_COORDINATOR_HPP
#define MPI_COORDINATOR_HPP

#include <cstdint>
#include <string>
#include <vector>
#include <cstdio>
#include <mpi.h>

namespace mn5 {

// Work assignment for a single MPI rank
struct WorkAssignment {
    int rank;                  // MPI rank ID
    __uint128_t start;         // Starting seed (2^71 + offset)
    __uint128_t end;           // Ending seed (exclusive)
    uint64_t total_seeds;      // Total seeds to process
    uint64_t total_candidates; // After mod-6 filter
    
    WorkAssignment() : rank(0), start(0), end(0), total_seeds(0), total_candidates(0) {}
};

// Global results aggregation
struct GlobalResults {
    uint64_t total_tested;
    uint64_t total_time_ms;
    uint64_t max_steps;
    __uint128_t max_steps_seed;
    __uint128_t max_peak;
    
    uint64_t total_overflow;
    uint64_t total_fuse;
    uint64_t total_cycles;
    
    GlobalResults() : total_tested(0), total_time_ms(0), max_steps(0),
                      max_steps_seed(0), max_peak(0), total_overflow(0),
                      total_fuse(0), total_cycles(0) {}
};

// Per-rank results (what workers send back to master)
struct RankResults {
    int rank;
    uint64_t tested;
    uint64_t time_ms;
    uint64_t max_steps;
    __uint128_t max_steps_seed;
    __uint128_t max_peak;
    
    uint64_t overflow_count;
    uint64_t fuse_count;
    uint64_t cycle_count;
    
    double throughput;  // nums/sec
};

// MPI Coordinator class
class MPICoordinator {
private:
    int rank_;
    int size_;
    bool is_master_;
    
public:
    MPICoordinator() : rank_(0), size_(1), is_master_(true) {
        MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
        MPI_Comm_size(MPI_COMM_WORLD, &size_);
        is_master_ = (rank_ == 0);
    }
    
    int rank() const { return rank_; }
    int size() const { return size_; }
    bool is_master() const { return is_master_; }
    bool is_worker() const { return !is_master_; }
    
    // Partition work across all ranks
    // Master (rank 0) coordinates, ranks 1..N-1 are workers
    std::vector<WorkAssignment> partition_work(
        __uint128_t base_start,
        uint64_t total_count
    ) {
        std::vector<WorkAssignment> assignments;
        
        // Number of worker ranks (exclude master rank 0)
        int num_workers = size_ - 1;
        if (num_workers <= 0) {
            fprintf(stderr, "[ERROR] Need at least 2 MPI ranks (1 master + 1 worker)\n");
            return assignments;
        }
        
        // Each worker gets equal chunk
        uint64_t chunk_size = total_count / num_workers;
        uint64_t remainder = total_count % num_workers;
        
        uint64_t current_offset = 0;
        
        for (int worker_rank = 1; worker_rank < size_; ++worker_rank) {
            WorkAssignment wa;
            wa.rank = worker_rank;
            
            // Give remainder to first few workers
            uint64_t this_chunk = chunk_size;
            if (worker_rank - 1 < (int)remainder) {
                this_chunk++;
            }
            
            wa.start = base_start + current_offset;
            wa.end = wa.start + this_chunk;
            wa.total_seeds = this_chunk;
            
            // Mod-6 filter: only test n ≡ 1,5 (mod 6)
            // Actual candidates ≈ total_seeds / 3
            wa.total_candidates = (this_chunk + 2) / 3; // Conservative estimate
            
            assignments.push_back(wa);
            current_offset += this_chunk;
        }
        
        return assignments;
    }
    
    // Get work assignment for current rank
    WorkAssignment get_my_work(const std::vector<WorkAssignment>& assignments) {
        for (const auto& wa : assignments) {
            if (wa.rank == rank_) {
                return wa;
            }
        }
        return WorkAssignment(); // Empty if not found
    }
    
    // Master: Broadcast work assignments to all workers
    void broadcast_work(const WorkAssignment& work) {
        // Pack work assignment into buffer
        struct WorkPacket {
            uint64_t start_lo, start_hi;
            uint64_t end_lo, end_hi;
            uint64_t total_seeds;
            uint64_t total_candidates;
        } packet;
        
        packet.start_lo = (uint64_t)work.start;
        packet.start_hi = (uint64_t)(work.start >> 64);
        packet.end_lo = (uint64_t)work.end;
        packet.end_hi = (uint64_t)(work.end >> 64);
        packet.total_seeds = work.total_seeds;
        packet.total_candidates = work.total_candidates;
        
        // Send to specific rank
        if (is_master_ && work.rank != 0) {
            MPI_Send(&packet, sizeof(packet), MPI_BYTE, work.rank, 0, MPI_COMM_WORLD);
        }
    }
    
    // Worker: Receive work assignment from master
    WorkAssignment receive_work() {
        struct WorkPacket {
            uint64_t start_lo, start_hi;
            uint64_t end_lo, end_hi;
            uint64_t total_seeds;
            uint64_t total_candidates;
        } packet;
        
        MPI_Status status;
        MPI_Recv(&packet, sizeof(packet), MPI_BYTE, 0, 0, MPI_COMM_WORLD, &status);
        
        WorkAssignment work;
        work.rank = rank_;
        work.start = ((__uint128_t)packet.start_hi << 64) | packet.start_lo;
        work.end = ((__uint128_t)packet.end_hi << 64) | packet.end_lo;
        work.total_seeds = packet.total_seeds;
        work.total_candidates = packet.total_candidates;
        
        return work;
    }
    
    // Worker: Send results back to master
    void send_results(const RankResults& results) {
        struct ResultPacket {
            int rank;
            uint64_t tested;
            uint64_t time_ms;
            uint64_t max_steps;
            uint64_t max_steps_seed_lo, max_steps_seed_hi;
            uint64_t max_peak_lo, max_peak_hi;
            uint64_t overflow_count;
            uint64_t fuse_count;
            uint64_t cycle_count;
            double throughput;
        } packet;
        
        packet.rank = results.rank;
        packet.tested = results.tested;
        packet.time_ms = results.time_ms;
        packet.max_steps = results.max_steps;
        packet.max_steps_seed_lo = (uint64_t)results.max_steps_seed;
        packet.max_steps_seed_hi = (uint64_t)(results.max_steps_seed >> 64);
        packet.max_peak_lo = (uint64_t)results.max_peak;
        packet.max_peak_hi = (uint64_t)(results.max_peak >> 64);
        packet.overflow_count = results.overflow_count;
        packet.fuse_count = results.fuse_count;
        packet.cycle_count = results.cycle_count;
        packet.throughput = results.throughput;
        
        MPI_Send(&packet, sizeof(packet), MPI_BYTE, 0, 1, MPI_COMM_WORLD);
    }
    
    // Master: Receive results from a worker
    RankResults receive_results(int from_rank) {
        struct ResultPacket {
            int rank;
            uint64_t tested;
            uint64_t time_ms;
            uint64_t max_steps;
            uint64_t max_steps_seed_lo, max_steps_seed_hi;
            uint64_t max_peak_lo, max_peak_hi;
            uint64_t overflow_count;
            uint64_t fuse_count;
            uint64_t cycle_count;
            double throughput;
        } packet;
        
        MPI_Status status;
        MPI_Recv(&packet, sizeof(packet), MPI_BYTE, from_rank, 1, MPI_COMM_WORLD, &status);
        
        RankResults results;
        results.rank = packet.rank;
        results.tested = packet.tested;
        results.time_ms = packet.time_ms;
        results.max_steps = packet.max_steps;
        results.max_steps_seed = ((__uint128_t)packet.max_steps_seed_hi << 64) | packet.max_steps_seed_lo;
        results.max_peak = ((__uint128_t)packet.max_peak_hi << 64) | packet.max_peak_lo;
        results.overflow_count = packet.overflow_count;
        results.fuse_count = packet.fuse_count;
        results.cycle_count = packet.cycle_count;
        results.throughput = packet.throughput;
        
        return results;
    }
    
    // Master: Aggregate results from all workers
    GlobalResults aggregate_results(const std::vector<RankResults>& all_results) {
        GlobalResults global;
        
        for (const auto& r : all_results) {
            global.total_tested += r.tested;
            global.total_time_ms += r.time_ms;
            global.total_overflow += r.overflow_count;
            global.total_fuse += r.fuse_count;
            global.total_cycles += r.cycle_count;
            
            if (r.max_steps > global.max_steps) {
                global.max_steps = r.max_steps;
                global.max_steps_seed = r.max_steps_seed;
            }
            
            if (r.max_peak > global.max_peak) {
                global.max_peak = r.max_peak;
            }
        }
        
        return global;
    }
    
    // Print work distribution summary (master only)
    void print_work_distribution(const std::vector<WorkAssignment>& assignments) {
        if (!is_master_) return;
        
        fprintf(stderr, "\n=== MPI Work Distribution ===\n");
        fprintf(stderr, "Total MPI ranks: %d (1 master + %d workers)\n", 
                size_, size_ - 1);
        fprintf(stderr, "\nWorker Assignments:\n");
        
        for (const auto& wa : assignments) {
            // Convert uint128_t to string for printing (simplified)
            fprintf(stderr, "  Rank %d: %lu seeds (~%lu candidates after mod-6)\n",
                    wa.rank, wa.total_seeds, wa.total_candidates);
        }
        fprintf(stderr, "\n");
    }
    
    // Print results summary (master only)
    void print_global_results(const GlobalResults& results, double wall_time_sec) {
        if (!is_master_) return;
        
        fprintf(stderr, "\n=== Global Results Summary ===\n");
        fprintf(stderr, "Total tested:     %lu numbers\n", results.total_tested);
        fprintf(stderr, "Wall time:        %.2f seconds\n", wall_time_sec);
        fprintf(stderr, "Total CPU time:   %.2f seconds\n", results.total_time_ms / 1000.0);
        fprintf(stderr, "Throughput:       %.2f M nums/sec\n", 
                results.total_tested / wall_time_sec / 1e6);
        fprintf(stderr, "Max steps:        %lu\n", results.max_steps);
        fprintf(stderr, "Total overflow:   %lu\n", results.total_overflow);
        fprintf(stderr, "Total fuse hits:  %lu\n", results.total_fuse);
        fprintf(stderr, "Total cycles:     %lu\n", results.total_cycles);
        fprintf(stderr, "\n");
    }
    
    // Barrier synchronization
    void barrier() {
        MPI_Barrier(MPI_COMM_WORLD);
    }
};

} // namespace mn5

#endif // MPI_COORDINATOR_HPP
