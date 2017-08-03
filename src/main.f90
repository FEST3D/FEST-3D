program main
  !-------------------------------------------------
  ! 170801 - Jatinder Pal Singh Sandhu
  ! change : - explicit call with use module instead of using whole module
  !          - include error and mpi file
  !------------------------------------------------

  use global_vars,  only: process_id
  use global_vars,  only: max_iters
  use global_vars,  only: current_iter
  use global_vars,  only: want_to_stop
  use global_vars,  only: tolerance
  use global_vars,  only: resnorm
  use solver     ,  only: step
  use solver     ,  only: converged
  use start_finish, only:  start_run
  use start_finish, only: finish_run

#include "error.inc"
#include "mpi.inc"

!--------Start---------!
  call start_run()

!  do while (.not. converged())
  resnorm=1
  do while ((current_iter <= max_iters) .and. (resnorm>tolerance))
     call step()! new name "iterate"
  end do

  call finish_run()
!--------Stop---------!

end program main
