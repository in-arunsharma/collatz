module purge
module load intel impi/2021.10.0
which mpiicpc && mpiicpc --version   # should be icpc classic wrapper

mpiicpc -Ofast -xCORE-AVX512 -qopenmp -flto -std=c++17 \
        -o V1.5-mpi-lean V1.5-mpi-lean.cpp

ldd ./V1.5-mpi-lean | egrep 'iomp|mpi' || true