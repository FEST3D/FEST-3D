module start_finish
  use fclose, only: close_all_files
  use solver, only: destroy_solver
  use solver, only: setup_solver
#ifdef __GFORTRAN__
  use mpi
#endif    
  implicit none
#ifdef __INTEL_COMPILER
  include "mpif.h"
#endif
  private

  public :: abort_run
  public :: finish_run
  public :: start_run

  contains

    subroutine abort_run()
      implicit none
      integer :: ierr

      call close_all_files()
      call destroy_solver()
      call MPI_FINALIZE(ierr)
      stop

    end subroutine abort_run

    subroutine finish_run()
      implicit none
      integer :: ierr

      call close_all_files()
      call destroy_solver()
      call MPI_FINALIZE(ierr)

    end subroutine finish_run

    subroutine start_run()
      implicit none
      integer :: ierr

      call MPI_INIT(ierr)
      call setup_solver()

    end subroutine start_run

end module start_finish
