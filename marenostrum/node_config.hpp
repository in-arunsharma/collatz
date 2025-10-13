// node_config.hpp - Hardware Detection and Configuration for MareNostrum 5
// Detects node type (GPP vs ACC) and configures optimal settings

#ifndef NODE_CONFIG_HPP
#define NODE_CONFIG_HPP

#include <cstdint>
#include <string>
#include <fstream>
#include <sstream>
#include <cstdlib>

#ifdef __CUDACC__
#include <cuda_runtime.h>
#endif

namespace mn5 {

// Node type enumeration
enum class NodeType {
    UNKNOWN,
    GPP,    // General Purpose Partition (CPU-only)
    ACC     // Accelerated Partition (GPU + CPU)
};

// Hardware configuration for a node
struct NodeConfig {
    NodeType type;
    int mpi_rank;
    int mpi_size;
    
    // CPU configuration
    int cpu_cores;           // Total physical cores
    int cpu_threads;         // Hyperthreads (if enabled)
    double cpu_freq_ghz;     // Base frequency
    uint64_t memory_gb;      // System memory
    
    // GPU configuration (ACC nodes only)
    int gpu_count;
    int gpu_sm_count;        // Streaming multiprocessors per GPU
    uint64_t gpu_memory_gb;  // Memory per GPU
    
    // Optimal thread counts
    int openmp_threads;      // Recommended OpenMP thread count
    
    // Node identification
    std::string hostname;
    std::string cpu_model;
    std::string gpu_model;
    
    NodeConfig() : type(NodeType::UNKNOWN), mpi_rank(0), mpi_size(1),
                   cpu_cores(0), cpu_threads(0), cpu_freq_ghz(0.0),
                   memory_gb(0), gpu_count(0), gpu_sm_count(0),
                   gpu_memory_gb(0), openmp_threads(0) {}
};

// Detect CUDA-capable GPUs
inline int detect_cuda_devices() {
#ifdef __CUDACC__
    int device_count = 0;
    cudaError_t err = cudaGetDeviceCount(&device_count);
    if (err != cudaSuccess) {
        return 0;
    }
    return device_count;
#else
    return 0;
#endif
}

// Get GPU information
inline bool get_gpu_info(int device_id, std::string& model, int& sm_count, uint64_t& memory_gb) {
#ifdef __CUDACC__
    cudaDeviceProp prop;
    if (cudaGetDeviceProperties(&prop, device_id) == cudaSuccess) {
        model = prop.name;
        sm_count = prop.multiProcessorCount;
        memory_gb = prop.totalGlobalMem / (1024ULL * 1024 * 1024);
        return true;
    }
#endif
    return false;
}

// Read CPU info from /proc/cpuinfo
inline bool get_cpu_info(int& cores, std::string& model, double& freq_ghz) {
    std::ifstream cpuinfo("/proc/cpuinfo");
    if (!cpuinfo.is_open()) return false;
    
    std::string line;
    cores = 0;
    bool found_model = false;
    
    while (std::getline(cpuinfo, line)) {
        if (line.find("processor") == 0) {
            cores++;
        } else if (!found_model && line.find("model name") == 0) {
            size_t colon = line.find(':');
            if (colon != std::string::npos) {
                model = line.substr(colon + 2);
                found_model = true;
                
                // Try to extract frequency from model name
                size_t at_pos = model.find('@');
                if (at_pos != std::string::npos) {
                    std::string freq_str = model.substr(at_pos + 1);
                    freq_ghz = std::stod(freq_str);
                }
            }
        }
    }
    
    return cores > 0;
}

// Read total system memory from /proc/meminfo
inline uint64_t get_system_memory_gb() {
    std::ifstream meminfo("/proc/meminfo");
    if (!meminfo.is_open()) return 0;
    
    std::string line;
    while (std::getline(meminfo, line)) {
        if (line.find("MemTotal") == 0) {
            std::istringstream iss(line);
            std::string label;
            uint64_t mem_kb;
            iss >> label >> mem_kb;
            return mem_kb / (1024 * 1024); // Convert KB to GB
        }
    }
    return 0;
}

// Get hostname
inline std::string get_hostname() {
    char hostname[256];
    if (gethostname(hostname, sizeof(hostname)) == 0) {
        return std::string(hostname);
    }
    return "unknown";
}

// Detect node type based on hardware
inline NodeType detect_node_type(int cpu_cores, int gpu_count) {
    // MareNostrum 5 specific detection:
    // GPP nodes: 112 cores (2× 56-core Sapphire Rapids 8480+)
    // ACC nodes: 80 cores (2× 40-core Sapphire Rapids 8460Y+) + 4 GPUs
    
    if (gpu_count >= 4) {
        // Has GPUs -> ACC node
        return NodeType::ACC;
    } else if (cpu_cores >= 100) {
        // High core count, no GPUs -> GPP node
        return NodeType::GPP;
    }
    
    return NodeType::UNKNOWN;
}

// Determine optimal OpenMP thread count
inline int get_optimal_thread_count(NodeType type, int cpu_cores) {
    // On MareNostrum 5:
    // GPP: Use all 112 cores
    // ACC: Use 80 cores (leave some headroom for GPU management)
    
    const char* omp_num_threads = std::getenv("OMP_NUM_THREADS");
    if (omp_num_threads) {
        return std::atoi(omp_num_threads);
    }
    
    switch (type) {
        case NodeType::GPP:
            return cpu_cores; // Use all cores
        case NodeType::ACC:
            return cpu_cores; // Use all cores (GPU kernels run independently)
        default:
            return cpu_cores > 0 ? cpu_cores : 1;
    }
}

// Main configuration detection function
inline NodeConfig detect_node_configuration(int mpi_rank, int mpi_size) {
    NodeConfig config;
    config.mpi_rank = mpi_rank;
    config.mpi_size = mpi_size;
    config.hostname = get_hostname();
    
    // Detect CPU
    config.cpu_cores = 0;
    if (!get_cpu_info(config.cpu_cores, config.cpu_model, config.cpu_freq_ghz)) {
        // Fallback: try to get from environment
        const char* omp_max = std::getenv("OMP_NUM_THREADS");
        if (omp_max) {
            config.cpu_cores = std::atoi(omp_max);
        } else {
            config.cpu_cores = 1; // Safe fallback
        }
    }
    
    // Detect system memory
    config.memory_gb = get_system_memory_gb();
    
    // Detect GPUs
    config.gpu_count = detect_cuda_devices();
    if (config.gpu_count > 0) {
        // Get info from first GPU (assume all GPUs are identical)
        get_gpu_info(0, config.gpu_model, config.gpu_sm_count, config.gpu_memory_gb);
    }
    
    // Determine node type
    config.type = detect_node_type(config.cpu_cores, config.gpu_count);
    
    // Set optimal thread count
    config.openmp_threads = get_optimal_thread_count(config.type, config.cpu_cores);
    
    return config;
}

// Print configuration (for debugging/logging)
inline void print_node_config(const NodeConfig& config, FILE* out = stderr) {
    const char* type_str = "UNKNOWN";
    switch (config.type) {
        case NodeType::GPP: type_str = "GPP (CPU-only)"; break;
        case NodeType::ACC: type_str = "ACC (GPU+CPU)"; break;
        default: break;
    }
    
    fprintf(out, "[MPI Rank %d/%d] Node Configuration:\n", 
            config.mpi_rank, config.mpi_size);
    fprintf(out, "  Hostname:     %s\n", config.hostname.c_str());
    fprintf(out, "  Node Type:    %s\n", type_str);
    fprintf(out, "  CPU Model:    %s\n", config.cpu_model.c_str());
    fprintf(out, "  CPU Cores:    %d\n", config.cpu_cores);
    fprintf(out, "  CPU Freq:     %.1f GHz\n", config.cpu_freq_ghz);
    fprintf(out, "  System RAM:   %lu GB\n", config.memory_gb);
    
    if (config.gpu_count > 0) {
        fprintf(out, "  GPU Count:    %d\n", config.gpu_count);
        fprintf(out, "  GPU Model:    %s\n", config.gpu_model.c_str());
        fprintf(out, "  GPU Memory:   %lu GB per GPU\n", config.gpu_memory_gb);
        fprintf(out, "  GPU SMs:      %d per GPU\n", config.gpu_sm_count);
    }
    
    fprintf(out, "  OpenMP Threads: %d\n", config.openmp_threads);
    fprintf(out, "\n");
}

// Validate configuration against expected MareNostrum 5 specs
inline bool validate_marenostrum_config(const NodeConfig& config) {
    switch (config.type) {
        case NodeType::GPP:
            // Expect: 112 cores, 256-1024 GB RAM, no GPUs
            if (config.cpu_cores < 100 || config.cpu_cores > 120) {
                fprintf(stderr, "[WARNING] GPP node has unexpected core count: %d (expected ~112)\n", 
                        config.cpu_cores);
                return false;
            }
            if (config.gpu_count > 0) {
                fprintf(stderr, "[WARNING] GPP node should not have GPUs\n");
                return false;
            }
            break;
            
        case NodeType::ACC:
            // Expect: 80 cores, 512 GB RAM, 4 Hopper GPUs
            if (config.cpu_cores < 70 || config.cpu_cores > 90) {
                fprintf(stderr, "[WARNING] ACC node has unexpected core count: %d (expected ~80)\n", 
                        config.cpu_cores);
                return false;
            }
            if (config.gpu_count != 4) {
                fprintf(stderr, "[WARNING] ACC node has %d GPUs (expected 4)\n", 
                        config.gpu_count);
                return false;
            }
            break;
            
        default:
            fprintf(stderr, "[ERROR] Unknown node type!\n");
            return false;
    }
    
    return true;
}

} // namespace mn5

#endif // NODE_CONFIG_HPP
