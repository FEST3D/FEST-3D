program main
  use mpi
  use global_vars
  use solver

  implicit none
!  include "mpif.h"

  integer :: ierr
  call MPI_INIT(ierr)

  call setup_solver()

  do while (.not. converged())
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
