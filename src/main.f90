program main
#ifdef __GFORTRAN__
 use mpi
#endif
  use global_vars
  use solver

  implicit none
#ifdef __INTEL_COMPILER
  include "mpif.h"
#endif

  integer :: ierr
  call MPI_INIT(ierr)

  call setup_solver()

!  do while (.not. converged())
  do while (.true.)
     !print *, 'Iteration', iter
!     call flush()
     call step()
     if (current_iter == max_iters) then
        exit
     end if
  end do
  call destroy_solver()
  call MPI_FINALIZE(ierr)
end program main
