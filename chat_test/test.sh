module purge
module load GCC/12.3.0 impi/2021.10.0
export I_MPI_CXX=g++                      # make mpicxx call g++
mpicxx -O3 -march=native -flto -fopenmp -DNDEBUG \
  -o V1.5-mpi-lean V15_MPI_Lean.cpp

which mpicxx
mpicxx -show