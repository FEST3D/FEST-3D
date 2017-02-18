module resnorm_
  !----------------------------------------------------
  ! this module contains subroutine that 
  ! 1. check if time for resnorm dump is arrived
  ! 2. calculate resnorm
  ! 3. send those resnorm to processor number 0
  ! 4. Recalulate resnorm based on information 
  !    availble from all processors
  ! 5. Append the data to resnorm file
  !----------------------------------------------------

  use global,      only: RESNORM_FILE_UNIT
  use global,      only: RES_CONTROL_FILE_UNIT
  use global,      only: res_control_file
  use global,      only: resnorm_file
  use global,      only: LONG_BUFFER_LENGTH
  use global,      only: STRING_BUFFER_LENGTH
  use global,       only: resnorm_number

  use global_vars, only: gm
  use global_vars, only: density_inf
  use global_vars, only: x_speed_inf
  use global_vars, only: y_speed_inf
  use global_vars, only: z_speed_inf
  use global_vars, only: pressure_inf
  use global_vars, only: tk_inf
  use global_vars, only: tw_inf
  use global_vars, only:   mass_residue
  use global_vars, only:  x_mom_residue
  use global_vars, only:  y_mom_residue
  use global_vars, only:  z_mom_residue
  use global_vars, only: energy_residue
  use global_vars, only:    TKE_residue
  use global_vars, only:  omega_residue
  use global_vars, only:        resnorm
  use global_vars, only:    vis_resnorm
  use global_vars, only:   turb_resnorm
  use global_vars, only:   cont_resnorm     
  use global_vars, only:  x_mom_resnorm    
  use global_vars, only:  y_mom_resnorm    
  use global_vars, only:  z_mom_resnorm    
  use global_vars, only: energy_resnorm   
  use global_vars, only:    TKE_resnorm    
  use global_vars, only:  omega_resnorm    
  use global_vars, only:        resnorm_d1
  use global_vars, only:    vis_resnorm_d1
  use global_vars, only:   turb_resnorm_d1 
  use global_vars, only:   cont_resnorm_d1
  use global_vars, only:  x_mom_resnorm_d1
  use global_vars, only:  y_mom_resnorm_d1
  use global_vars, only:  z_mom_resnorm_d1
  use global_vars, only: energy_resnorm_d1
  use global_vars, only:    TKE_resnorm_d1
  use global_vars, only:  omega_resnorm_d1
  use global_vars, only:        resnorm_0
  use global_vars, only:    vis_resnorm_0
  use global_vars, only:   turb_resnorm_0 
  use global_vars, only:   cont_resnorm_0
  use global_vars, only:  x_mom_resnorm_0
  use global_vars, only:  y_mom_resnorm_0
  use global_vars, only:  z_mom_resnorm_0
  use global_vars, only: energy_resnorm_0
  use global_vars, only:    TKE_resnorm_0
  use global_vars, only:  omega_resnorm_0
  use global_vars, only:        resnorm_0s
  use global_vars, only:    vis_resnorm_0s
  use global_vars, only:   turb_resnorm_0s 
  use global_vars, only:   cont_resnorm_0s
  use global_vars, only:  x_mom_resnorm_0s
  use global_vars, only:  y_mom_resnorm_0s
  use global_vars, only:  z_mom_resnorm_0s
  use global_vars, only: energy_resnorm_0s
  use global_vars, only:    TKE_resnorm_0s
  use global_vars, only:  omega_resnorm_0s
  use global_vars, only: current_iter
  use global_vars, only: res_write_interval
  use global_vars, only: write_percision

  use global_vars, only: turbulence

  use utils,      only: dmsg
  use utils,      only: dealloc
  use utils,      only: alloc
  use layout,     only: process_id
  use layout,     only: total_process
  use string

  use res_viscous,  only: compute_viscous_resnorm
  use res_turbulent, only: compute_turbulent_resnorm

#ifdef __GFORTRAN
  use mpi
#endif

  implicit none
#ifdef __INTEL_COMPILER
  include "mpif.h"
#endif
  private

  integer  :: write_num
  integer, dimension(resnorm_number, 2) :: write_permission
  real,    dimension(:), allocatable, target :: write_data
  character(len=STRING_BUFFER_LENGTH), dimension(:), allocatable :: write_name

  public :: write_resnorm
  public :: setup_resnorm
  public :: destroy_resnorm

  contains

    subroutine write_resnorm()
      implicit none
      character(len=20) :: frm
!      character(len=*), parameter :: delimiter = "    "

      call dmsg(1, 'resnorm_', 'write_resnorm')
      call compute_resnorm()

      if (process_id == 0) then
        write(frm, '(A, I0, A)') "(I0, ", write_num, "e20.10E2, 4x)"
        write(RESNORM_FILE_UNIT, frm) current_iter, write_data(1:write_num) 
      end if
    end subroutine write_resnorm

    subroutine read_permission_to_write()
      implicit none 
      integer :: i

      call dmsg(1, 'resnorm_', 'read_permision_to_write')

      open(RES_CONTROL_FILE_UNIT, FILE=res_control_file, STATUS='old', ACTION='read')

      read(RES_CONTROL_FILE_UNIT,*)
      do i = 1,resnorm_number
        read(RES_CONTROL_FILE_UNIT,*) write_permission(i, 1:2)
      end do

      close(RES_CONTROL_FILE_UNIT)

    end subroutine read_permission_to_write

    subroutine compute_resnorm()
      implicit none
      call dmsg(1, 'resnorm_', 'compute_resnorm')


      call compute_viscous_resnorm()
      if (turbulence /= 'none') then
        call compute_turbulent_resnorm()
      end if
     
      if (process_id ==0 ) then
      resnorm     = sqrt(                        &
                           vis_resnorm    **2  + &
                          turb_resnorm    **2    &
                        )
      resnorm_0s   = sqrt(                        &
                           vis_resnorm_0s  **2  + &
                          turb_resnorm_0s  **2    &
                        )
      resnorm_d1  = resnorm / resnorm_0s

      end if

    end subroutine compute_resnorm

    subroutine setup_resnorm()
      implicit none
      integer :: i
      integer :: forward_count=0
      integer :: back_count = resnorm_number*2
      character(len=LONG_BUFFER_LENGTH) :: header

      call dmsg(1, 'resnorm_', 'setup_resnorm')

      call read_permission_to_write()

      write_num = sum(write_permission)

      
      call alloc(write_data, 1,resnorm_number*2,&
                errmsg='Unable to allocate memory to write_data in resnorm module')
      allocate(write_name(1:write_num))

      if (write_permission(1,1) == 1) then
        forward_count = forward_count + 1
        write_name(forward_count) = 'resnorm '
        resnorm  => write_data(forward_count)
      else
        resnorm  => write_data(back_count)
        back_count = back_count - 1
      end if

      if (write_permission(1,2) == 1) then
        forward_count = forward_count + 1
        write_name(forward_count) = 'resnorm_norm1 '
        resnorm_d1  => write_data(forward_count)
      else
        resnorm_d1  => write_data(back_count)
        back_count = back_count - 1
      end if

      if (write_permission(2,1) == 1) then
        forward_count = forward_count + 1
        write_name(forward_count) = 'vis_resnorm '
        vis_resnorm  => write_data(forward_count)
      else
        vis_resnorm  => write_data(back_count)
        back_count = back_count - 1
      end if

      if (write_permission(2,2) == 1) then
        forward_count = forward_count + 1
        write_name(forward_count) = 'vis_resnorm_norm1 '
        vis_resnorm_d1  => write_data(forward_count)
      else
        vis_resnorm_d1  => write_data(back_count)
        back_count = back_count - 1
      end if

      if (write_permission(3,1) == 1) then
        forward_count = forward_count + 1
        write_name(forward_count) = 'turb_resnorm '
        turb_resnorm  => write_data(forward_count)
      else
        turb_resnorm  => write_data(back_count)
        back_count = back_count - 1
      end if

      if (write_permission(3,2) == 1) then
        forward_count = forward_count + 1
        write_name(forward_count) = 'turb_resnorm_norm1 '
        turb_resnorm_d1  => write_data(forward_count)
      else
        turb_resnorm_d1  => write_data(back_count)
        back_count = back_count - 1
      end if

      if (write_permission(4,1) == 1) then
        forward_count = forward_count + 1
        write_name(forward_count) = 'cont_resnorm '
        cont_resnorm  => write_data(forward_count)
      else
        cont_resnorm  => write_data(back_count)
        back_count = back_count - 1
      end if

      if (write_permission(4,2) == 1) then
        forward_count = forward_count + 1
        write_name(forward_count) = 'cont_resnorm_norm1 '
        cont_resnorm_d1  => write_data(forward_count)
      else
        cont_resnorm_d1  => write_data(back_count)
        back_count = back_count - 1
      end if

      if (write_permission(5,1) == 1) then
        forward_count = forward_count + 1
        write_name(forward_count) = 'x_mom_resnorm '
        x_mom_resnorm  => write_data(forward_count)
      else
        x_mom_resnorm  => write_data(back_count)
        back_count = back_count - 1
      end if

      if (write_permission(5,2) == 1) then
        forward_count = forward_count + 1
        write_name(forward_count) = 'x_mom_resnorm_norm1 '
        x_mom_resnorm_d1  => write_data(forward_count)
      else
        x_mom_resnorm_d1  => write_data(back_count)
        back_count = back_count - 1
      end if

      if (write_permission(6,1) == 1) then
        forward_count = forward_count + 1
        write_name(forward_count) = 'y_mom_resnorm '
        y_mom_resnorm  => write_data(forward_count)
      else
        y_mom_resnorm  => write_data(back_count)
        back_count = back_count - 1
      end if

      if (write_permission(6,2) == 1) then
        forward_count = forward_count + 1
        write_name(forward_count) = 'y_mom_resnorm_norm1 '
        y_mom_resnorm_d1  => write_data(forward_count)
      else
        y_mom_resnorm_d1  => write_data(back_count)
        back_count = back_count - 1
      end if

      if (write_permission(7,1) == 1) then
        forward_count = forward_count + 1
        write_name(forward_count) = 'z_mom_resnorm '
        z_mom_resnorm  => write_data(forward_count)
      else
        z_mom_resnorm  => write_data(back_count)
        back_count = back_count - 1
      end if

      if (write_permission(7,2) == 1) then
        forward_count = forward_count + 1
        write_name(forward_count) = 'z_mom_resnorm_norm1 '
        z_mom_resnorm_d1  => write_data(forward_count)
      else
        z_mom_resnorm_d1  => write_data(back_count)
        back_count = back_count - 1
      end if

      if (write_permission(8,1) == 1) then
        forward_count = forward_count + 1
        write_name(forward_count) = 'energy_resnorm '
        energy_resnorm  => write_data(forward_count)
      else
        energy_resnorm  => write_data(back_count)
        back_count = back_count - 1
      end if

      if (write_permission(8,2) == 1) then
        forward_count = forward_count + 1
        write_name(forward_count) = 'energy_resnorm_norm1 '
        energy_resnorm_d1  => write_data(forward_count)
      else
        energy_resnorm_d1  => write_data(back_count)
        back_count = back_count - 1
      end if

      if (write_permission(9,1) == 1) then
        forward_count = forward_count + 1
        write_name(forward_count) = 'TKE_resnorm '
        TKE_resnorm  => write_data(forward_count)
      else
        TKE_resnorm  => write_data(back_count)
        back_count = back_count - 1
      end if

      if (write_permission(9,2) == 1) then
        forward_count = forward_count + 1
        write_name(forward_count) = 'TKE_resnorm_norm1 '
        TKE_resnorm_d1  => write_data(forward_count)
      else
        TKE_resnorm_d1  => write_data(back_count)
        back_count = back_count - 1
      end if

      if (write_permission(10,1) == 1) then
        forward_count = forward_count + 1
        write_name(forward_count) = 'omega_resnorm '
        omega_resnorm  => write_data(forward_count)
      else
        omega_resnorm  => write_data(back_count)
        back_count = back_count - 1
      end if

      if (write_permission(10,2) == 1) then
        forward_count = forward_count + 1
        write_name(forward_count) = 'omega_resnorm_norm1 '
        omega_resnorm_d1  => write_data(forward_count)
      else
        omega_resnorm_d1  => write_data(back_count)
        back_count = back_count - 1
      end if

      header = " iter "
      do i = 1, write_num
        write(header,'(A)') trim(header)//"  "//trim(write_name(i))
      end do
      
      if (process_id == 0) then
      open(RESNORM_FILE_UNIT, file=resnorm_file, ACTION='write')
      write(RESNORM_FILE_UNIT, '(A)') header
      end if


      ! intializing whole resnorm to 0.
      write_data = 0.
!             resnorm_0   = 0.
!         vis_resnorm_0   = 0.
!        turb_resnorm_0   = 0.
!        cont_resnorm_0   = 0.
!       x_mom_resnorm_0   = 0.
!       y_mom_resnorm_0   = 0.
!       z_mom_resnorm_0   = 0.
!      energy_resnorm_0   = 0.
!         TKE_resnorm_0   = 0.
!       omega_resnorm_0   = 0.


      end subroutine setup_resnorm

      subroutine destroy_resnorm()
        
        call dmsg(1, 'resnorm_', 'destroy_resnorm')
        deallocate(write_name)
        call dealloc(write_data)
        if (process_id == 0) then
          close(RESNORM_FILE_UNIT)
        end if

      end subroutine destroy_resnorm



end module resnorm_

