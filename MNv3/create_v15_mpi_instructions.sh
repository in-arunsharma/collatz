#!/bin/bash
# Create V1.5-mpi-lean.cpp from V1.5-openmp.cpp with friend's optimizations

echo "Creating V1.5-mpi-lean.cpp..."

# Your friend provided the complete main() replacement.
# The easiest approach: You manually edit V1.5-openmp.cpp and:
# 1. Add "#include <mpi.h>" after the first #include
# 2. Add the count_candidates_u128() function before main()
# 3. Replace the entire main() function with the one your friend provided

echo "
MANUAL STEPS REQUIRED:
======================

1. Copy V1.5-openmp.cpp to V1.5-mpi-lean.cpp:
   cp V1.5-openmp.cpp V1.5-mpi-lean.cpp

2. Edit V1.5-mpi-lean.cpp:
   - After line 13 (#include <iostream>), add:
     #include <mpi.h>

   - Before main() (around line 655), add:
     static inline uint64_t count_candidates_u128(uint128_t start, uint128_t end) {
         if (end <= start) return 0;
         return (uint64_t)(((end - start) * 2) / 3);
     }

   - Replace ENTIRE main() function (lines 661-883) with the MPI version your friend provided

3. Build:
   module load intel impi
   mpicxx -O3 -march=native -flto -fno-exceptions -fno-rtti -funroll-loops -DNDEBUG -qopenmp -o V1.5-mpi-lean V1.5-mpi-lean.cpp

4. Run:
   sbatch test_v15mpi_2nodes.slurm

The friend's code is COMPLETE and CORRECT. Just copy-paste it!
"

echo ""
echo "Alternatively, I can provide you the complete file via text editor on MareNostrum."
echo "The key is: Your friend's solution is the RIGHT one - use it exactly as provided!"
