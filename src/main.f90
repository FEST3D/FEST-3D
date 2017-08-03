program main
  !-------------------------------------------------
  ! 170801 - Jatinder Pal Singh Sandhu
  ! change : - explicit call with use module instead of using whole module
  !          - include error and mpi file
  ! 170803 - jatinder Pal Singh Sandhu
  ! change : - new name to step -> iterate_one_more_time_step
  !------------------------------------------------

  use global_vars,  only: process_id
  use global_vars,  only: max_iters
  use global_vars,  only: current_iter
  use global_vars,  only: want_to_stop
  use global_vars,  only: tolerance
  use global_vars,  only: resnorm
  use solver     ,  only: iterate_one_more_time_step
  use convergence,  only: converged
  use start_finish, only:  start_run
  use start_finish, only: finish_run

#include "error.inc"
#include "mpi.inc"

!--------Start---------!
  call start_run()

  do while ((current_iter <= max_iters) .and. (.not. converged()))
     call iterate_one_more_time_step()
  end do

  call finish_run()
!--------Stop---------!

end program main
