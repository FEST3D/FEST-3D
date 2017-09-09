module read_multi_level_vtk
  
  !---------------------------------------------------------
  ! This module read state + other variable in output file
  !---------------------------------------------------------
  use global     , only :    IN_FILE_UNIT
  use global     , only : OUTIN_FILE_UNIT
  use global     , only : outin_file

  use global_vars, only : read_data_format
  use global_vars, only : read_file_format
  use global_vars, only : start_from
  use global_vars, only : process_id
  use global_vars, only : imax => imx
  use global_vars, only : jmax => jmx
  use global_vars, only : kmax => kmx
  use global_vars, only : grid_x
  use global_vars, only : grid_y
  use global_vars, only : grid_z
  use global_vars, only : density
  use global_vars, only : x_speed
  use global_vars, only : y_speed
  use global_vars, only : z_speed
  use global_vars, only : pressure
  use global_vars, only : tk
  use global_vars, only : tw
  use global_vars, only : tkl
  use global_vars, only : density_inf
  use global_vars, only : x_speed_inf
  use global_vars, only : y_speed_inf
  use global_vars, only : z_speed_inf
  use global_vars, only : pressure_inf 
  use global_vars, only : mu
  use global_vars, only : mu_t
  use global_vars, only : gm
  use global_vars, only : dist
  use global_vars, only : turbulence
  use global_vars, only : vis_resnorm
  use global_vars, only : cont_resnorm
  use global_vars, only : x_mom_resnorm
  use global_vars, only : y_mom_resnorm
  use global_vars, only : z_mom_resnorm
  use global_vars, only : energy_resnorm
  use global_vars, only : resnorm
  use global_vars, only :   mass_residue
  use global_vars, only :  x_mom_residue
  use global_vars, only :  y_mom_residue
  use global_vars, only :  z_mom_residue
  use global_vars, only : energy_residue
  use global_vars, only :        resnorm_0
  use global_vars, only :    vis_resnorm_0
  use global_vars, only :   turb_resnorm_0
  use global_vars, only :   cont_resnorm_0
  use global_vars, only :  x_mom_resnorm_0
  use global_vars, only :  y_mom_resnorm_0
  use global_vars, only :  z_mom_resnorm_0
  use global_vars, only : energy_resnorm_0
  use global_vars, only :    TKE_resnorm_0
  use global_vars, only :  omega_resnorm_0

  use global_vars, only : read_level
  use global_vars, only : mu_ref
  use global_vars, only : r_count
  use global_vars, only : r_list
  use global_vars, only : previous_flow_type

  use utils
  use string

  implicit none
  private
  integer :: i,j,k
  integer :: imx,jmx, kmx

  public :: read_file

  contains

    subroutine read_file()
      implicit none
      integer :: n

      imx = (imax-1)/2 + 1
      jmx = (jmax-1)/2 + 1
      kmx = (kmax-1)/2 + 1

      call read_header()
      call read_grid()
      do n = 1,r_count

        select case (trim(r_list(n)))
        
          case('Velocity')
            call read_velocity()

          case('Density')
            call read_scalar(density, 'Density', -2)
          
          case('Pressure')
            call read_scalar(pressure, 'Pressure', -2)
            
          case('TKE')
            call read_scalar(tk, 'TKE', -2)

          case('Omega')
            call read_scalar(tw, 'Omega', -2)

          case('Kl')
            call read_scalar(tkl, 'Kl', -2)

          case('do not read')
            !skip 
            continue

          case Default
            Print*, "read error: list var : "//trim(r_list(n))

        end select
      end do

    end subroutine read_file

    subroutine read_header()
      implicit none

      call dmsg(1, 'read_output_vtk', 'read_header')
      read(IN_FILE_UNIT, *) !'# vtk DataFile Version 3.1'
      read(IN_FILE_UNIT, *) !'cfd-iitm output'   ! comment line
      read(IN_FILE_UNIT, *) !trim(read_data_format)
      read(IN_FILE_UNIT, *) !'DATASET STRUCTURED_GRID'
      read(IN_FILE_UNIT, *)


    end subroutine read_header

    subroutine read_grid()
      implicit none

      ! read grid point coordinates
      call dmsg(1, 'read_output_vtk', 'read_grid')
      read(IN_FILE_UNIT, * ) !'DIMENSIONS ', imx, ' ', jmx, ' ', kmx
      read(IN_FILE_UNIT, * ) !'POINTS ', imx*jmx*kmx, ' DOUBLE'
      do k = 1, kmx
       do j = 1, jmx
        do i = 1, imx
          read(IN_FILE_UNIT, *) !grid_x(i, j, k), ' ', grid_y(i, j, k), ' ', grid_z(i, j, k)
        end do
       end do
      end do
      read(IN_FILE_UNIT, *) 

    end subroutine read_grid

    subroutine read_velocity()
      implicit none

      call dmsg(1, 'read_output_vtk', 'read_velocity')
      read(IN_FILE_UNIT, *) !'CELL_DATA ', (imx-1)*(jmx-1)*(kmx-1)
      read(IN_FILE_UNIT, *) !'VECTORS Velocity FLOAT'
      do k = 1, kmx - 1
       do j = 1, jmx - 1
        do i = 1, imx - 1
          print*, i,j,k
          read(IN_FILE_UNIT, *) x_speed(2*i, 2*j, 2*k), y_speed(2*i, 2*j, 2*k), z_speed(2*i, 2*j, 2*k)
          !x_speed(2*i-1,2*j-1,2*k-1) = x_speed(2*i,2*j,2*k)!0.5*(x_speed(2*i-1,2*j-2,2*k-2) + x_speed(2*i,2*j,2*k))
          !y_speed(2*i-1,2*j-1,2*k-1) = y_speed(2*i,2*j,2*k)!0.5*(y_speed(2*i-1,2*j-2,2*k-2) + y_speed(2*i,2*j,2*k))
          !z_speed(2*i-1,2*j-1,2*k-1) = z_speed(2*i,2*j,2*k)!0.5*(z_speed(2*i-1,2*j-2,2*k-2) + z_speed(2*i,2*j,2*k))
          !x_speed(2*i-0,2*j-1,2*k-1) = x_speed(2*i,2*j,2*k)!0.5*(x_speed(2*i-1,2*j-2,2*k-2) + x_speed(2*i,2*j,2*k))
          !y_speed(2*i-0,2*j-1,2*k-1) = y_speed(2*i,2*j,2*k)!0.5*(y_speed(2*i-1,2*j-2,2*k-2) + y_speed(2*i,2*j,2*k))
          !z_speed(2*i-0,2*j-1,2*k-1) = z_speed(2*i,2*j,2*k)!0.5*(z_speed(2*i-1,2*j-2,2*k-2) + z_speed(2*i,2*j,2*k))
          !x_speed(2*i-1,2*j-0,2*k-1) = x_speed(2*i,2*j,2*k)!0.5*(x_speed(2*i-1,2*j-2,2*k-2) + x_speed(2*i,2*j,2*k))
          !y_speed(2*i-1,2*j-0,2*k-1) = y_speed(2*i,2*j,2*k)!0.5*(y_speed(2*i-1,2*j-2,2*k-2) + y_speed(2*i,2*j,2*k))
          !z_speed(2*i-1,2*j-0,2*k-1) = z_speed(2*i,2*j,2*k)!0.5*(z_speed(2*i-1,2*j-2,2*k-2) + z_speed(2*i,2*j,2*k))
          !x_speed(2*i-1,2*j-1,2*k-0) = x_speed(2*i,2*j,2*k)!0.5*(x_speed(2*i-1,2*j-2,2*k-2) + x_speed(2*i,2*j,2*k))
          !y_speed(2*i-1,2*j-1,2*k-0) = y_speed(2*i,2*j,2*k)!0.5*(y_speed(2*i-1,2*j-2,2*k-2) + y_speed(2*i,2*j,2*k))
          !z_speed(2*i-1,2*j-1,2*k-0) = z_speed(2*i,2*j,2*k)!0.5*(z_speed(2*i-1,2*j-2,2*k-2) + z_speed(2*i,2*j,2*k))
          !x_speed(2*i-0,2*j-1,2*k-0) = x_speed(2*i,2*j,2*k)!0.5*(x_speed(2*i-1,2*j-2,2*k-2) + x_speed(2*i,2*j,2*k))
          !y_speed(2*i-0,2*j-1,2*k-0) = y_speed(2*i,2*j,2*k)!0.5*(y_speed(2*i-1,2*j-2,2*k-2) + y_speed(2*i,2*j,2*k))
          !z_speed(2*i-0,2*j-1,2*k-0) = z_speed(2*i,2*j,2*k)!0.5*(z_speed(2*i-1,2*j-2,2*k-2) + z_speed(2*i,2*j,2*k))
          !x_speed(2*i-1,2*j-0,2*k-0) = x_speed(2*i,2*j,2*k)!0.5*(x_speed(2*i-1,2*j-2,2*k-2) + x_speed(2*i,2*j,2*k))
          !y_speed(2*i-1,2*j-0,2*k-0) = y_speed(2*i,2*j,2*k)!0.5*(y_speed(2*i-1,2*j-2,2*k-2) + y_speed(2*i,2*j,2*k))
          !z_speed(2*i-1,2*j-0,2*k-0) = z_speed(2*i,2*j,2*k)!0.5*(z_speed(2*i-1,2*j-2,2*k-2) + z_speed(2*i,2*j,2*k))
          !x_speed(2*i-0,2*j-0,2*k-1) = x_speed(2*i,2*j,2*k)!0.5*(x_speed(2*i-1,2*j-2,2*k-2) + x_speed(2*i,2*j,2*k))
          !y_speed(2*i-0,2*j-0,2*k-1) = y_speed(2*i,2*j,2*k)!0.5*(y_speed(2*i-1,2*j-2,2*k-2) + y_speed(2*i,2*j,2*k))
          !z_speed(2*i-0,2*j-0,2*k-1) = z_speed(2*i,2*j,2*k)!0.5*(z_speed(2*i-1,2*j-2,2*k-2) + z_speed(2*i,2*j,2*k))
        end do
       end do
      end do
      read(IN_FILE_UNIT, *) 

    end subroutine read_velocity

    subroutine read_scalar(var, name,  index)
      implicit none
      integer, intent(in) :: index
      real, dimension(index:imx-index,index:jmx-index,index:kmx-index), intent(out) :: var
      character(len=*), intent(in) :: name

      call dmsg(1, 'read_output_vtk', trim(name))
      read(IN_FILE_UNIT, *) !'SCALARS '//trim(name)//' FLOAT'
      read(IN_FILE_UNIT, *) !'LOOKUP_TABLE default'
      do k = 1, kmx - 1
       do j = 1, jmx - 1
        do i = 1, imx - 1
          read(IN_FILE_UNIT, *) var(2*i, 2*j, 2*k)
          !var(2*i-1,2*j-1,2*k-1) = 0.5*(var(2*i-2,2*j-2,2*k-2) + var(2*i,2*j,2*k))
          var(2*i-1,2*j-1,2*k-1) = var(2*i,2*j,2*k)
          var(2*i-1,2*j-0,2*k-1) = var(2*i,2*j,2*k)
          var(2*i-0,2*j-1,2*k-1) = var(2*i,2*j,2*k)
          var(2*i-1,2*j-1,2*k-0) = var(2*i,2*j,2*k)
          var(2*i-1,2*j-0,2*k-0) = var(2*i,2*j,2*k)
          var(2*i-0,2*j-1,2*k-0) = var(2*i,2*j,2*k)
          var(2*i-0,2*j-0,2*k-1) = var(2*i,2*j,2*k)
        end do
       end do
      end do
      read(IN_FILE_UNIT, *)

    end subroutine read_scalar


end module read_multi_level_vtk
