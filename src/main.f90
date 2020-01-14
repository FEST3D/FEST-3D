program main
  !< Main program
  !-------------------------------------------------
  use solver     ,  only: iterate_one_more_time_step
  use solver     ,  only: control
  use convergence,  only: converged
  use solver, only:  start_run
  use solver, only: finish_run

!--------Start---------!
  call start_run()

  do while ((control%current_iter <= control%max_iters) .and. (.not. converged(control)) .and. (.not. control%Halt))
     call iterate_one_more_time_step()
  end do

  call finish_run()
!--------Stop---------!

end program main
