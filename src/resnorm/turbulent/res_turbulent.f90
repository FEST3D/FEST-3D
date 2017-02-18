module Res_turbulent
  !----------------------------------------------------
  ! this module contains subroutine that 
  ! 1. check if time for resnorm dump is arrived
  ! 2. calculate resnorm
  ! 3. send those resnorm to processor number 0
  ! 4. Recalulate resnorm based on information 
  !    availble from all processors
  ! 5. Append the data to resnorm file
  !----------------------------------------------------

  use global_vars, only: density_inf
  use global_vars, only: tk_inf
  use global_vars, only: tw_inf
  use global_vars, only: TKE_residue
  use global_vars, only: omega_residue
  use global_vars, only:  turb_resnorm
  use global_vars, only:   TKE_resnorm
  use global_vars, only: omega_resnorm
  use global_vars, only:   TKE_resnorm_0
  use global_vars, only:  turb_resnorm_0
  use global_vars, only: omega_resnorm_0
  use global_vars, only:   TKE_resnorm_0s
  use global_vars, only:  turb_resnorm_0s
  use global_vars, only: omega_resnorm_0s
  use global_vars, only:   TKE_resnorm_d1
  use global_vars, only:  turb_resnorm_d1
  use global_vars, only: omega_resnorm_d1
  use global_vars, only: current_iter
  use global_vars, only: res_write_interval
  use global_vars, only: turbulence

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
  integer, parameter                :: n = 6    !total number of turbulent variables
  real                              :: speed_inf
  real, dimension(n)                :: res_send_buf
  real, dimension(:),   allocatable :: root_res_recv_buf
  real, dimension(:,:), allocatable :: global_resnorm

  public :: compute_turbulent_resnorm

  contains

    subroutine compute_turbulent_resnorm()
      implicit none

        call compute_block_resnorm()
        if (current_iter == 1) then
          call store_intial_resnorm()
        end if
        if (process_id == 0) then
          call alloc(root_res_recv_buf, 1,total_process*n, &
              errmsg='Error: Unable to allocate memory to root_res_recv_buf')
          root_res_recv_buf = 0.
        end if
        call send_resnorm_to_process_0()
        if (process_id == 0) then
          call alloc(global_resnorm, 1,total_process, 1,n, &
              errmsg='Error: Unable to allocate memory to root_res_recv_buf')
          call recv_resnorm_to_process_0()
          call recalculate_collective_resnorm()
          call dealloc(global_resnorm)
          call dealloc(root_res_recv_buf)
        end if

    end subroutine compute_turbulent_resnorm

    subroutine compute_block_resnorm()

        implicit none
        
        call dmsg(1, 'res_turbuent', 'compute_block_norm')

        select case (turbulence)
          
          case ('sst')
            call R_TKE()
            call R_omega()

          case DEFAULT
            call dmsg(5, "res_turbulence", "compute_block_resnorm", &
                       "ERROR: Turbulence model not recognised")
            STOP

          end select

        turb_resnorm =     (                      &
                            TKE_resnorm    + &
                            omega_resnorm    &
                           )

    end subroutine compute_block_resnorm


    subroutine R_TKE()
      implicit none

      TKE_resnorm = sum(                               &
                        (TKE_residue(:, :, :)/         &
                        (density_inf * tk_inf)) ** 2   &
                       )                               

    end subroutine R_TKE

    subroutine R_omega()
      implicit none

      omega_resnorm = sum(                              &
                          (omega_residue(:, :, :)/      &
                          (density_inf * tw_inf)) ** 2  &
                         )                              

    end subroutine R_omega

    subroutine store_intial_resnorm()
      implicit none
  
         turb_resnorm_0  =  turb_resnorm
          TKE_resnorm_0  =   TKE_resnorm
        omega_resnorm_0  = omega_resnorm

    end subroutine store_intial_resnorm

    subroutine send_resnorm_to_process_0()
      implicit none
      integer :: ierr

      res_send_buf(1) =  turb_resnorm
      res_send_buf(2) =   TKE_resnorm
      res_send_buf(3) = omega_resnorm
      res_send_buf(4) =  turb_resnorm_0
      res_send_buf(5) =   TKE_resnorm_0
      res_send_buf(6) = omega_resnorm_0

      call MPI_Gather(res_send_buf, n, MPI_DOUBLE_PRECISION, &
        root_res_recv_buf, n, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

    end subroutine send_resnorm_to_process_0

    subroutine recv_resnorm_to_process_0()
      implicit none
      integer :: id
      do id = 0,total_process-1
        global_resnorm(id+1, 1) = root_res_recv_buf(1+n*id)
        global_resnorm(id+1, 2) = root_res_recv_buf(2+n*id)
        global_resnorm(id+1, 3) = root_res_recv_buf(3+n*id)
        global_resnorm(id+1, 4) = root_res_recv_buf(4+n*id)
        global_resnorm(id+1, 5) = root_res_recv_buf(5+n*id)
        global_resnorm(id+1, 6) = root_res_recv_buf(6+n*id)
      end do
    end  subroutine recv_resnorm_to_process_0

    subroutine recalculate_collective_resnorm()
      implicit none
!      real  ::  turb_resnorm_0s
!      real  ::   TKE_resnorm_0s
!      real  :: omega_resnorm_0s

         turb_resnorm     = SUM(global_resnorm(: ,1))
          TKE_resnorm     = SUM(global_resnorm(: ,2))
        omega_resnorm     = SUM(global_resnorm(: ,3))
         turb_resnorm_0s  = SUM(global_resnorm(: ,4))
          TKE_resnorm_0s  = SUM(global_resnorm(: ,5))
        omega_resnorm_0s  = SUM(global_resnorm(: ,6))

         turb_resnorm     = SQRT( turb_resnorm  )
          TKE_resnorm     = SQRT(  TKE_resnorm  )
        omega_resnorm     = SQRT(omega_resnorm  )
         turb_resnorm_0s  = SQRT( turb_resnorm_0s)
          TKE_resnorm_0s  = SQRT(  TKE_resnorm_0s)
        omega_resnorm_0s  = SQRT(omega_resnorm_0s)

         turb_resnorm_d1  =  turb_resnorm / turb_resnorm_0s
          TKE_resnorm_d1  =   TKE_resnorm /  TKE_resnorm_0s
        omega_resnorm_d1  = omega_resnorm /omega_resnorm_0s

    end subroutine recalculate_collective_resnorm

end module res_turbulent

