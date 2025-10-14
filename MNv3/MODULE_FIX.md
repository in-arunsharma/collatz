## 🔧 Module Name Fix for MareNostrum 5

### Problem
`gcc/12.2.0` doesn't exist on MareNostrum 5

### Solution
Run this on MareNostrum to find available modules:

```bash
./find_modules.sh
```

Or manually:
```bash
module avail gcc
module avail mpi
```

### Updated Scripts
All scripts now try `gcc/13.1.0` and `openmpi/4.1.5` (adjust if needed).

### Build Options

**Option 1: OpenMPI + GCC (RECOMMENDED)**
```bash
module purge
module load gcc/13.1.0 openmpi/4.1.5
mpicxx -O3 -march=native -flto -fopenmp ... -o V1.5-mpi-lean V1.5-mpi-lean.cpp
```

**Option 2: Intel MPI + GCC backend**
```bash
module purge
module load gcc impi/2021.10.0
export I_MPI_CXX=g++
mpicxx -O3 -march=native -flto -fopenmp ... -o V1.5-mpi-lean V1.5-mpi-lean.cpp
```

**Option 3: Use whatever GCC is available**
```bash
# Find GCC versions
module spider gcc

# Load any version (13.x, 12.x, 11.x all work)
module load gcc/<version>

# Then load matching MPI
module load openmpi  # or impi with I_MPI_CXX=g++
```

### Commands on MareNostrum

```bash
cd /gpfs/projects/nct_352/nct01225/collatz/MNv3/

# Find available modules
./find_modules.sh

# If gcc/13.1.0 exists, build:
./build_v15_mpi_lean.sh

# If NOT, edit build_v15_mpi_lean.sh with correct version, then:
./build_v15_mpi_lean.sh

# Then test:
sbatch test_v15mpi_1node_sanity.slurm
```

### Key Point
**ANY GCC version works** - just avoid icpc (Intel Classic Compiler)!

The important part is:
1. Use GCC (any version 9+)
2. NOT icpc (Intel Classic)
3. With `-march=native -O3 -flto`
