module resnorm
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
  use global,      only: LONG_BUFFER_STRING

  use global_vars, only: gm
  use global_vars, only: density_inf
  use global_vars, only: x_speed_inf
  use global_vars, only: y_speed_inf
  use global_vars, only: z_speed_inf
  use global_vars, only: pressure_inf
  use global_vars, only: tk_inf
  use global_vars, only: tw_inf
  use global_vars, only: mass_residue
  use global_vars, only: x_mom_residue
  use global_vars, only: y_mom_residue
  use global_vars, only: z_mom_residue
  use global_vars, only: energy_residue
  use global_vars, only: TKE_residue
  use global_vars, only: omega_residue
  use global_vars, only: vis_resnorm
  use global_vars, only: vis_resnorm_0
  use global_vars, only: vis_resnorm_abs
  use global_vars, only: cont_resnorm
  use global_vars, only: cont_resnorm_0
  use global_vars, only: x_mom_resnorm
  use global_vars, only: x_mom_resnorm_0
  use global_vars, only: y_mom_resnorm
  use global_vars, only: y_mom_resnorm_0
  use global_vars, only: z_mom_resnorm
  use global_vars, only: z_mom_resnorm_0
  use global_vars, only: energy_resnorm
  use global_vars, only: energy_resnorm_0
  use global_vars, only: TKE_resnorm
  use global_vars, only: TKE_resnorm_0
  use global_vars, only: omega_resnorm
  use global_vars, only: omega_resnorm_0
  use global_vars, only: current_iter
  use global_vars, only: res_write_interval

  use global_vars, only: resnorm_number

  use utils,      only: dmsg
  use utils,      only: dealloc
  use utils,      only: alloc
  use layout,     only: process_id
  use layout,     only: total_process

  use res_viscous  only: compute_viscous_resnorm
  use res_turbuent only: compute_turbulent_resnorm


  use mpi

  implicit none
  private

  public :: write_resnorm
  integer  :: write_num
  integer, dimension(resnorm_number), allocatable :: write_permission
  real, dimension(:), allocatable, pointer :: write_data
  character(len=STRING_BUFFER_LENGTH), dimension(:), allocatable :: write_name

  contains

    subroutine write_resnorm()
      implicit none
      character(len=*), dimension(write_num) :: delimiter="    "

      call compute_resnorm()

       write(RESNORM_FILE_UNIT, '(f0.16, A)') (/write_data(i), delimiter(i), i=1,write_num/) 
    end subroutine write_resnorm

    subroutine read_permission_to_write()
      implicit none 
      integer :: i

      open(RES_CONTROL_FILE_UNIT, FILE=res_control_file, STATUS='old', ACTION='read')

      do i = 1,resnorm_number
        read(RES_CONTROL_FILE_UNIT,*) write_permission(i)
      end do

      close(RES_CONTROL_FILE_UNIT)

    end subroutine read_permission_to_write

    subroutine compute_resnorm()

      resnorm           = 0.
      resnorm_abs       = 0.
      vis_resnorm       = 0.
      vis_resnorm_abs   = 0.
      turb_resnorm      = 0.
      turb_resnorm_abs  = 0.
      cont_resnorm      = 0.
      x_mom_resnorm     = 0.
      y_mom_resnorm     = 0.
      z_mom_resnorm     = 0.
      energy_resnorm    = 0.
      TKE_resnorm       = 0.
      omega_resnorm     = 0.

      call compute_viscous_resnorm()
      if (turbulence /= 'none') then
        call compute_turbulent_resnorm()
      end if
      
      resnorm     = sqrt(                         &
                          vis_resnorm      ** 2 + &
                          turb_resnorm     ** 2   &
                        )

      resnorm_abs = sqrt(                         &
                          vis_resnorm_abs  ** 2 + &
                          turb_resnorm_abs ** 2   &
                        )

    end subroutine compute_resnorm

    subroutine setup_resnorm()
      implicit none
      integer :: i
      integer :: count=0
      character(len=LONG_BUFFER_STRING) :: header

      call read_permission_to_write()

      write_num = sum(write_permission)

      call alloc(write_data, 1, write_num,&
                errmsg='Unable to allocate memory to write_data in resnorm module')
      call alloc(write_name, 1, write_num,&
                errmsg='Unable to allocate memory to write_data in resnorm module')

      if (write_permission(1) == 1) then
        count = count + 1
        write_name(count) = 'resnorm '
        write_data(count) => resnorm
      end if
      if (write_permission(2) == 1)
        count = count + 1
        write_name(count) = 'resnorm_abs '
        write_data(count) => resnorm_abs
      end if
      if (write_permission(3) == 1)
        count = count + 1
        write_name(count) = 'vis_resnorm '
        write_data(count) => vis_resnorm
      end if
      if (write_permission(4) == 1)
        count = count + 1
        write_name(count) = 'vis_resnorm_abs '
        write_data(count) => vis_resnorm_abs
      end if
      if (write_permission(5) == 1)
        count = count + 1
        write_name(count) = 'turb_resnorm '
        write_data(count) => turb_resnorm
      end if
      if (write_permission(6) == 1)
        count = count + 1
        write_name(count) = 'turb_resnorm_abs '
        write_data(count) => turb_resnorm_abs
      end if
      if (write_permission(7) == 1)
        count = count + 1
        write_name(count) = 'cont_resnorm '
        write_data(count) => cont_resnorm
      end if
      if (write_permission(8) == 1)
        count = count + 1
        write_name(count) = 'x_mom_resnorm '
        write_data(count) => x_mom_resnorm
      end if
      if (write_permission(9) == 1)
        count = count + 1
        write_name(count) = 'y_mom_resnorm '
        write_data(count) => y_mom_resnorm
      end if
      if (write_permission(10) == 1)
        count = count + 1
        write_name(count) = 'z_mom_resnorm '
        write_data(count) => z_mom_resnorm
      end if
      if (write_permission(11) == 1)
        count = count + 1
        write_name(count) = 'energy_resnorm '
        write_data(count) => energy_resnorm
      end if
      if (write_permission(12) == 1)
        count = count + 1
        write_name(count) = 'TKE_resnorm '
        write_data(count) => TKE_resnorm
      end if
      if (write_permission(13) == 1)
        count = count + 1
        write_name(count) = 'omega_resnorm '
        write_data(count) => omega_resnorm
      end if

      header = ""
      do i = 1, write_num
        write(header,'(A)') header//trim(write_name(i))
      end do


      open(RESNORM_FILE_UNIT, file=resnorm_file ACTION='write')
      write(RESNORM_FILE_UNIT, '(A)') header

      end subroutine setup_resnorm

      subroutine destroy_resnorm()
        
        call dealloc(write_name)
        call dealloc(write_data)
        close(RESNORM_FILE_UNIT)

      end subroutine destroy_resnorm



end module resnorm

