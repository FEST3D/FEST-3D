program main
  use mpi
  use solver

  implicit none

  integer :: ierr
  call MPI_INIT(ierr)

  call setup_solver()

  do while (.not. converged())
     !print *, 'Iteration', iter
     call flush()
     call step()
     if (iter == max_iters) then
        exit
     end if
  end do
  call destroy_solver()
  call MPI_FINALIZE(ierr)
end program main
