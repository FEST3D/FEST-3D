#ifdef __GFORTRAN__
    use mpi
#endif    
    implicit none
#ifdef __INTEL_COMPILER
    include "mpif.h"
#endif
