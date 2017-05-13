program main
#ifdef __GFORTRAN__
 use mpi
#endif
  use global_vars
  use solver
  use start_finish, only:  start_run
  use start_finish, only: finish_run

  implicit none
#ifdef __INTEL_COMPILER
  include "mpif.h"
#endif

  call start_run()

!  do while (.not. converged())
  do while (.true.)
     call step()
     if (current_iter == max_iters) then
        exit
     end if
     if (want_to_stop==1) max_iters=current_iter+1
  end do

  call finish_run()

end program main
