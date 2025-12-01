#!/bin/bash
# Compile and run CUDA version
# Requires: NVIDIA GPU with CUDA toolkit installed

echo "Compiling CUDA version..."
nvcc -O3 -arch=sm_60 -o V1_cuda V1_cuda.cu

if [ $? -eq 0 ]; then
    echo "Compilation successful!"
    echo "Running on GPU..."
    ./V1_cuda
else
    echo "Compilation failed. Make sure CUDA toolkit is installed."
    echo "Install with: sudo apt install nvidia-cuda-toolkit"
fi
