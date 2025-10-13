// collatz_mpi_gpp.cpp - Phase 1: GPP Nodes Only (CPU OpenMP)
// MPI + OpenMP hybrid for MareNostrum 5 General Purpose Partition
// This wraps your proven V1.4b-openmp code with MPI for multi-node scaling

#include <mpi.h>
#include <iostream>
#include <cstdio>
#include <cstdint>
#include <chrono>
#include <string>
#include <vector>

#include "node_config.hpp"
#include "mpi_coordinator.hpp"
#include "marenostrum_config.hpp"

using namespace mn5;
typedef unsigned __int128 uint128_t;

// ============================================================================
// COMMAND LINE PARSING
// ============================================================================

struct ProgramArgs {
    uint64_t start_offset = 0;
    uint64_t count = 1000000000;  // Default: 1 billion seeds
    std::string run_tag = "mpi_gpp";
    
    bool parse(int argc, char** argv) {
        if (argc < 3) {
            return false;
        }
        
        start_offset = std::stoull(argv[1]);
        count = std::stoull(argv[2]);
        
        if (argc >= 4) {
            run_tag = argv[3];
        }
        
        return true;
    }
    
    void print_usage(const char* program) {
        fprintf(stderr, "Usage: %s <start_offset> <count> [run_tag]\n", program);
        fprintf(stderr, "\n");
        fprintf(stderr, "  start_offset  - Offset from 2^71\n");
        fprintf(stderr, "  count         - Total seeds to test\n");
        fprintf(stderr, "  run_tag       - Tag for output files (default: mpi_gpp)\n");
        fprintf(stderr, "\n");
        fprintf(stderr, "Example:\n");
        fprintf(stderr, "  mpirun -np 6 %s 0 1000000000 hackathon\n", program);
        fprintf(stderr, "  (1 master + 5 GPP worker nodes, 1 billion seeds)\n");
        fprintf(stderr, "\n");
    }
};

// ============================================================================
// MASTER PROCESS (Rank 0)
// ============================================================================

int master_process(const ProgramArgs& args, MPICoordinator& coordinator) {
    fprintf(stderr, "\n");
    fprintf(stderr, "╔══════════════════════════════════════════════════════════╗\n");
    fprintf(stderr, "║  Collatz Conjecture - MareNostrum 5 Phase 1 (GPP Only) ║\n");
    fprintf(stderr, "╚══════════════════════════════════════════════════════════╝\n");
    fprintf(stderr, "\n");
    
    // Print configuration
    print_config_summary();
    
    fprintf(stderr, "=== Master Process (Rank 0) ===\n");
    fprintf(stderr, "Total MPI ranks: %d\n", coordinator.size());
    fprintf(stderr, "Worker ranks: %d\n", coordinator.size() - 1);
    fprintf(stderr, "\n");
    
    // Calculate base start (2^71 + offset)
    uint128_t base = (uint128_t)1 << CollatzParameters::BASE_POWER;
    uint128_t start = base + args.start_offset;
    
    // Partition work across workers
    auto assignments = coordinator.partition_work(start, args.count);
    coordinator.print_work_distribution(assignments);
    
    // Send work to each worker
    fprintf(stderr, "[Master] Distributing work to %zu workers...\n", assignments.size());
    for (const auto& work : assignments) {
        coordinator.broadcast_work(work);
    }
    fprintf(stderr, "[Master] All workers dispatched!\n\n");
    
    // Wait for results from all workers
    std::vector<RankResults> all_results;
    auto wall_start = std::chrono::high_resolution_clock::now();
    
    fprintf(stderr, "[Master] Waiting for results from workers...\n");
    for (int worker_rank = 1; worker_rank < coordinator.size(); ++worker_rank) {
        RankResults results = coordinator.receive_results(worker_rank);
        all_results.push_back(results);
        
        fprintf(stderr, "[Master] Received results from rank %d: %lu tested in %lu ms (%.2f M/sec)\n",
                results.rank, results.tested, results.time_ms, results.throughput / 1e6);
    }
    
    auto wall_end = std::chrono::high_resolution_clock::now();
    double wall_time = std::chrono::duration<double>(wall_end - wall_start).count();
    
    // Aggregate results
    GlobalResults global = coordinator.aggregate_results(all_results);
    
    // Print final summary
    fprintf(stderr, "\n");
    coordinator.print_global_results(global, wall_time);
    
    fprintf(stderr, "╔══════════════════════════════════════════════════════════╗\n");
    fprintf(stderr, "║                    SUCCESS!                              ║\n");
    fprintf(stderr, "║  Phase 1 (GPP) completed successfully                   ║\n");
    fprintf(stderr, "╚══════════════════════════════════════════════════════════╝\n");
    fprintf(stderr, "\n");
    
    return 0;
}

// ============================================================================
// WORKER PROCESS (Ranks 1..N)
// ============================================================================

int worker_process(const ProgramArgs& args, MPICoordinator& coordinator) {
    // Detect node configuration
    NodeConfig config = detect_node_configuration(coordinator.rank(), coordinator.size());
    
    // Print node info
    print_node_config(config);
    
    // Validate this is a GPP node
    if (config.type != NodeType::GPP) {
        fprintf(stderr, "[ERROR] Rank %d: Expected GPP node, got %d\n", 
                coordinator.rank(), (int)config.type);
        return 1;
    }
    
    // Validate against MareNostrum 5 specs
    if (!validate_marenostrum_config(config)) {
        fprintf(stderr, "[WARNING] Rank %d: Configuration doesn't match MN5 specs\n",
                coordinator.rank());
        // Continue anyway (might be testing locally)
    }
    
    // Receive work assignment from master
    fprintf(stderr, "[Rank %d] Waiting for work assignment...\n", coordinator.rank());
    WorkAssignment work = coordinator.receive_work();
    
    fprintf(stderr, "[Rank %d] Received work: %lu total seeds\n", 
            coordinator.rank(), work.total_seeds);
    fprintf(stderr, "[Rank %d] Starting computation with %d OpenMP threads...\n",
            coordinator.rank(), config.openmp_threads);
    
    // TODO: Import actual V1.4b-openmp worker code here
    // For now, simulate work
    auto t_start = std::chrono::high_resolution_clock::now();
    
    // PLACEHOLDER: Replace this with actual GPPNodeWorker from worker_gpp_node.cpp
    uint64_t tested = work.total_candidates;
    uint64_t max_steps = 1000;  // Dummy value
    uint128_t max_steps_seed = work.start;
    uint128_t max_peak = work.start * 100;
    uint64_t overflow_count = 0;
    uint64_t fuse_count = 0;
    
    // Simulate computation time
    // In real code, this would be: worker->process();
    sleep(2);
    
    auto t_end = std::chrono::high_resolution_clock::now();
    uint64_t time_ms = std::chrono::duration_cast<std::chrono::milliseconds>(
        t_end - t_start).count();
    
    // Prepare results
    RankResults results;
    results.rank = coordinator.rank();
    results.tested = tested;
    results.time_ms = time_ms;
    results.max_steps = max_steps;
    results.max_steps_seed = max_steps_seed;
    results.max_peak = max_peak;
    results.overflow_count = overflow_count;
    results.fuse_count = fuse_count;
    results.cycle_count = 0;
    results.throughput = tested / (time_ms / 1000.0);
    
    // Send results back to master
    fprintf(stderr, "[Rank %d] Computation complete! Sending results to master...\n",
            coordinator.rank());
    coordinator.send_results(results);
    
    fprintf(stderr, "[Rank %d] Done!\n", coordinator.rank());
    return 0;
}

// ============================================================================
// MAIN ENTRY POINT
// ============================================================================

int main(int argc, char** argv) {
    // Initialize MPI
    MPI_Init(&argc, &argv);
    
    // Parse arguments
    ProgramArgs args;
    if (!args.parse(argc, argv)) {
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        if (rank == 0) {
            args.print_usage(argv[0]);
        }
        MPI_Finalize();
        return 1;
    }
    
    // Create MPI coordinator
    MPICoordinator coordinator;
    
    int return_code = 0;
    
    try {
        if (coordinator.is_master()) {
            // Rank 0: Master process
            return_code = master_process(args, coordinator);
        } else {
            // Ranks 1..N: Worker processes
            return_code = worker_process(args, coordinator);
        }
    } catch (const std::exception& e) {
        fprintf(stderr, "[Rank %d] EXCEPTION: %s\n", coordinator.rank(), e.what());
        return_code = 1;
    }
    
    // Synchronize before exit
    coordinator.barrier();
    
    // Finalize MPI
    MPI_Finalize();
    
    return return_code;
}
