module mass_error
  !----------------------------------------------------
  ! this module contains subroutine that 
  !----------------------------------------------------

  use global,      only: RESNORM_FILE_UNIT

  use global_vars, only: imx
  use global_vars, only: jmx
  use global_vars, only: kmx
  use global_vars, only: density
  use global_vars, only: x_speed
  use global_vars, only: y_speed
  use global_vars, only: z_speed
  use global_vars, only: volume
  use global_vars, only: xA
  use global_vars, only: yA
  use global_vars, only: zA
  use global_vars, only: xnx
  use global_vars, only: xny
  use global_vars, only: xnz
  use global_vars, only: ynx
  use global_vars, only: yny
  use global_vars, only: ynz
  use global_vars, only: znx
  use global_vars, only: zny
  use global_vars, only: znz
  use global_vars, only: merror
  use global_vars, only: face_names
  use global_vars, only: id
  use global_vars, only: current_iter

  use utils,      only: dmsg
  use utils,      only: dealloc
  use utils,      only: alloc
  use layout,     only: process_id
  use layout,     only: total_process
#ifdef __GFORTRAN__
  use mpi
#endif

  implicit none
#ifdef __INTEL_COMPILER
  include "mpif.h"
#endif
  private

  integer :: face_num
  real, dimension(:), allocatable :: recv_buf
  real :: dmass=0.


  public :: calculate_mass_error

  contains
    subroutine calculate_mass_error()
      implicit none
      dmass=0.
      if(process_id==0) call alloc(recv_buf, 0, total_process-1)
      dmass=merror
      call mpi_gather_mass_error()
      if(process_id==0) call dealloc(recv_buf)
    end subroutine calculate_mass_error

    subroutine mpi_gather_mass_error()
      implicit none
      integer :: ierr
    
        call MPI_Gather(dmass, 1, MPI_DOUBLE_PRECISION,&
                        recv_buf , 1, MPI_DOUBLE_PRECISION,&
                        0,MPI_COMM_WORLD, ierr)
      if(process_id==0) merror=abs(sum(recv_buf))

    end subroutine mpi_gather_mass_error

      

end module mass_error
