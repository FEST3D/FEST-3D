module read_output_vtk
  
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
  use global_vars, only : imx
  use global_vars, only : jmx
  use global_vars, only : kmx
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

  use global_vars, only : mu_ref
  use global_vars, only : n_write
  use global_vars, only : rw_list
  use global_vars, only : previous_flow_type

  use utils
  use string

  implicit none
  private
  integer :: i,j,k
  real    :: speed_inf
  integer :: counter
  integer :: read_control_file_flag=0
  character(len=8) :: file_format
  character(len=16) :: data_format
  character(len=16) :: read_flow_type
  

  public :: read_file

  contains

    subroutine read_file()
      implicit none
      integer :: n

      call read_header()
      call read_grid()
      do n = 1,n_write

        select case (trim(rw_list(n)))
        
          case('U')
            call read_velocity()

          case('Density')
            call read_density()
          
          case('Pressure')
            call read_pressure()
            
          case('Mu')
            if (mu_ref/=0.0) then
              call read_mu()
            else
              print*, "Read error: Asked to read non-existing variable- "//trim(rw_list(n))
            end if
            
          case('Mu_t')
            if (turbulence/='none') then
              call read_mu_t()
            else
              print*, "Read error: Asked to read non-existing variable- "//trim(rw_list(n))
            end if
            
          case('TKE')
            if(turbulence=="sst" .and. previous_flow_type=="sst")then
            call read_TKE()
            else
              print*, "Read error: Asked to read non-existing variable- "//trim(rw_list(n))
            end if

          case('Omega')
            if(turbulence=="sst" .and. previous_flow_type=="sst") then
            call read_omega()
            else
              print*, "Read error: Asked to read non-existing variable- "//trim(rw_list(n))
            end if

          case('Wall_distance')
            if(turbulence/="none" .and. previous_flow_type/="viscous") then
            call read_dist()
            else
              print*, "Read error: Asked to read non-existing variable- "//trim(rw_list(n))
            end if

          case('Resnorm')
            call read_resnorm()

          case Default
            print*, "read_error: cannot read variable "//trim(rw_list(n))//" to file"

        end select
      end do

    end subroutine read_file

!
!    subroutine read_turbulence_variables()
!      implicit none
!
!      select case (read_flow_type)
!
!        case ('viscous')
!          !do nothing
!          !restart turbulent varibale with infinity condition !todo
!          continue
!        case ('sst')
!          call read_TKE()
!          call read_omega()
!        case DEFAULT
!          call dmsg(5,'read_output_vtk', 'read_turbulence_variables',&
!            'ERROR: Read flow type not recognised')
!          STOP
!      end select
!
!    end subroutine read_turbulence_variables
!
!    subroutine read_viscosity()
!      implicit none
!
!      if (mu_ref/=0.0) then
!        call read_mu()
!      end if
!
!      if (turbulence/='none') then
!        call read_mu_t()
!      end if
!
!    end subroutine read_viscosity

    subroutine read_header()
      implicit none

      call dmsg(1, 'read_output_vtk', 'read_header')

      if (read_data_format == "ASCII") then
        read(IN_FILE_UNIT, *) !'# vtk DataFile Version 3.1'
        read(IN_FILE_UNIT, *) !'cfd-iitm output'   ! comment line
        read(IN_FILE_UNIT, *) !trim(read_data_format)
        read(IN_FILE_UNIT, *) !'DATASET STRUCTURED_GRID'
        read(IN_FILE_UNIT, *)
      elseif (read_data_format == 'BINARY') then
        read(IN_FILE_UNIT) !'# vtk DataFile Version 3.1'
        read(IN_FILE_UNIT) !'cfd-iitm output'
        read(IN_FILE_UNIT) !trim(read_data_format)
        read(IN_FILE_UNIT) !'DATASET STRUCTURED_GRID'
        read(IN_FILE_UNIT)
      end if


    end subroutine read_header

    subroutine read_grid()
      implicit none

      ! read grid point coordinates
      call dmsg(1, 'read_output_vtk', 'read_grid')
      if (read_data_format == "ASCII") then
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
      elseif (read_data_format == 'BINARY') then
        read(IN_FILE_UNIT) !'DIMENSIONS ', imx, ' ', jmx, ' ', kmx
        read(IN_FILE_UNIT) !'POINTS ', imx*jmx*kmx, ' DOUBLE'
        do k = 1, kmx
         do j = 1, jmx
          do i = 1, imx
              read(IN_FILE_UNIT) !grid_x(i, j, k), ' ', grid_y(i, j, k), ' ', grid_z(i, j, k)
          end do
         end do
        end do
        read(IN_FILE_UNIT) 
      end if

    end subroutine read_grid

    subroutine read_velocity()
      implicit none
      call dmsg(1, 'read_output_vtk', 'read_velocity')

        ! Cell data
        ! Writing Velocity
      if (read_data_format == "ASCII") then
        read(IN_FILE_UNIT, *) !'CELL_DATA ', (imx-1)*(jmx-1)*(kmx-1)
        read(IN_FILE_UNIT, *) !'VECTORS Velocity FLOAT'
        do k = 1, kmx - 1
         do j = 1, jmx - 1
          do i = 1, imx - 1
            read(IN_FILE_UNIT, *) x_speed(i, j, k), y_speed(i, j, k), z_speed(i, j, k)
          end do
         end do
        end do
        read(IN_FILE_UNIT, *) 
      elseif (read_data_format == 'BINARY') then
        read(IN_FILE_UNIT) !'CELL_DATA ', (imx-1)*(jmx-1)*(kmx-1)
        read(IN_FILE_UNIT) !'VECTORS Velocity DOUBLE'
        do k = 1, kmx - 1
         do j = 1, jmx - 1
          do i = 1, imx - 1
            read(IN_FILE_UNIT) x_speed(i, j, k), y_speed(i, j, k), z_speed(i, j, k)
          end do
         end do
        end do
        read(IN_FILE_UNIT) 
      end if

    end subroutine read_velocity

    subroutine read_density()
      implicit none

      call dmsg(1, 'read_output_vtk', 'read_density')
      ! Writing Density
      if (read_data_format == "ASCII") then
        read(IN_FILE_UNIT, *) !'SCALARS Density FLOAT'
        read(IN_FILE_UNIT, *) !'LOOKUP_TABLE default'
        do k = 1, kmx - 1
         do j = 1, jmx - 1
          do i = 1, imx - 1
            read(IN_FILE_UNIT, *) density(i, j, k)
          end do
         end do
        end do
        read(IN_FILE_UNIT, *) 
      elseif (read_data_format == 'BINARY') then
        read(IN_FILE_UNIT) !'SCALARS Density DOUBLE'
        read(IN_FILE_UNIT) !'LOOKUP_TABLE default'
        do k = 1, kmx - 1
         do j = 1, jmx - 1
          do i = 1, imx - 1
            read(IN_FILE_UNIT) density(i, j, k)
          end do
         end do
        end do
        read(IN_FILE_UNIT)
      end if

    end subroutine read_density

    subroutine read_pressure()
      implicit none

      call dmsg(1, 'read_output_vtk', 'read_pressure')
      ! Writing Pressure
      if (read_data_format == "ASCII") then
        read(IN_FILE_UNIT, *) !'SCALARS Pressure FLOAT'
        read(IN_FILE_UNIT, *) !'LOOKUP_TABLE default'
        do k = 1, kmx - 1
         do j = 1, jmx - 1
          do i = 1, imx - 1
            read(IN_FILE_UNIT, *) pressure(i, j, k)
          end do
         end do
        end do
        read(IN_FILE_UNIT, *) 
      elseif (read_data_format == 'BINARY') then
        read(IN_FILE_UNIT) !'SCALARS Pressure DOUBLE'
        read(IN_FILE_UNIT) !'LOOKUP_TABLE default'
        do k = 1, kmx - 1
         do j = 1, jmx - 1
          do i = 1, imx - 1
            read(IN_FILE_UNIT) pressure(i, j, k)
          end do
         end do
        end do
        read(IN_FILE_UNIT) 
      end if

    end subroutine read_pressure

    subroutine read_TKE()
      implicit none

      call dmsg(1, 'read_output_vtk', 'read_TKE')
      ! Writing Pressure
      if (read_data_format == "ASCII") then
        read(IN_FILE_UNIT, *) !'SCALARS k FLOAT'
        read(IN_FILE_UNIT, *) !'LOOKUP_TABLE default'
        do k = 1, kmx - 1
         do j = 1, jmx - 1
          do i = 1, imx - 1
            read(IN_FILE_UNIT, *) tk(i, j, k)
          end do
         end do
        end do
        read(IN_FILE_UNIT, *) 
      elseif (read_data_format == 'BINARY') then
        read(IN_FILE_UNIT) !'SCALARS k DOUBLE'
        read(IN_FILE_UNIT) !'LOOKUP_TABLE default'
        do k = 1, kmx - 1
         do j = 1, jmx - 1
          do i = 1, imx - 1
            read(IN_FILE_UNIT) tk(i, j, k)
          end do
         end do
        end do
        read(IN_FILE_UNIT) 
      end if

    end subroutine read_TKE

    subroutine read_omega()
      implicit none

      call dmsg(1, 'read_output_vtk', 'read_omega')
      ! Writing Pressure
      if (read_data_format == "ASCII") then
        read(IN_FILE_UNIT, *) !'SCALARS Omega FLOAT'
        read(IN_FILE_UNIT, *) !'LOOKUP_TABLE default'
        do k = 1, kmx - 1
         do j = 1, jmx - 1
          do i = 1, imx - 1
            read(IN_FILE_UNIT, *) tw(i, j, k)
          end do
         end do
        end do
        read(IN_FILE_UNIT, *) 
      elseif (read_data_format == 'BINARY') then
        read(IN_FILE_UNIT) !'SCALARS Omega DOUBLE'
        read(IN_FILE_UNIT) !'LOOKUP_TABLE default'
        do k = 1, kmx - 1
         do j = 1, jmx - 1
          do i = 1, imx - 1
            read(IN_FILE_UNIT) tw(i, j, k)
          end do
         end do
        end do
        read(IN_FILE_UNIT) 
      end if

    end subroutine read_omega

    subroutine read_mu()
      implicit none

      call dmsg(1, 'read_output_vtk', 'read_mu')
      ! Writing Pressure
      if (read_data_format == "ASCII") then
        read(IN_FILE_UNIT, *) !'SCALARS mu FLOAT'
        read(IN_FILE_UNIT, *) !'LOOKUP_TABLE default'
        do k = 1, kmx - 1
         do j = 1, jmx - 1
          do i = 1, imx - 1
            read(IN_FILE_UNIT, *) mu(i, j, k)
          end do
         end do
        end do
        read(IN_FILE_UNIT, *) 
      elseif (read_data_format == 'BINARY') then
        read(IN_FILE_UNIT) !'SCALARS mu DOUBLE'
        read(IN_FILE_UNIT) !'LOOKUP_TABLE default'
        do k = 1, kmx - 1
         do j = 1, jmx - 1
          do i = 1, imx - 1
            read(IN_FILE_UNIT) mu(i, j, k)
          end do
         end do
        end do
        read(IN_FILE_UNIT) 
      end if

    end subroutine read_mu

    subroutine read_mu_t()
      implicit none

      call dmsg(1, 'read_output_vtk', 'read_mu_t')
      ! Writing Pressure
      if (read_data_format == "ASCII") then
        read(IN_FILE_UNIT, *) !'SCALARS mu_t FLOAT'
        read(IN_FILE_UNIT, *) !'LOOKUP_TABLE default'
        do k = 1, kmx - 1
         do j = 1, jmx - 1
          do i = 1, imx - 1
            read(IN_FILE_UNIT, *) mu_t(i, j, k)
          end do
         end do
        end do
        read(IN_FILE_UNIT, *) 
      elseif (read_data_format == 'BINARY') then
        read(IN_FILE_UNIT) !'SCALARS mu_t DOUBLE'
        read(IN_FILE_UNIT) !'LOOKUP_TABLE default'
        do k = 1, kmx - 1
         do j = 1, jmx - 1
          do i = 1, imx - 1
            read(IN_FILE_UNIT) mu_t(i, j, k)
          end do
         end do
        end do
        read(IN_FILE_UNIT) 
      end if

    end subroutine read_mu_t

! Auxilary subroutine for special case
    subroutine read_dist()
      implicit none

      call dmsg(1, 'read_output_vtk', 'read_dist')
      ! Writing wall distance for each cell
      if (read_data_format == "ASCII") then
        read(IN_FILE_UNIT, *) !'SCALARS dist FLOAT'
        read(IN_FILE_UNIT, *) !'LOOKUP_TABLE default'
        do k = 1, kmx - 1
         do j = 1, jmx - 1
          do i = 1, imx - 1
            read(IN_FILE_UNIT, *) dist(i, j, k)
          end do
         end do
        end do
        read(IN_FILE_UNIT, *) 
      elseif (read_data_format == 'BINARY') then
        read(IN_FILE_UNIT) !'SCALARS dist DOUBLE'
        read(IN_FILE_UNIT) !'LOOKUP_TABLE default'
        do k = 1, kmx - 1
         do j = 1, jmx - 1
          do i = 1, imx - 1
            read(IN_FILE_UNIT) dist(i, j, k)
          end do
         end do
        end do
        read(IN_FILE_UNIT) 
      end if

    end subroutine read_dist

    subroutine read_resnorm()
      implicit none

      call dmsg(1, 'read_output_vtk', 'read_resnorm')
      ! Writing resnorm for each cell
      if (read_data_format == "ASCII") then
        read(IN_FILE_UNIT, *) !'SCALARS Resnorm FLOAT'
        read(IN_FILE_UNIT, *) !'LOOKUP_TABLE default'
        do k = 1, kmx - 1
         do j = 1, jmx - 1
          do i = 1, imx - 1
            read(IN_FILE_UNIT, *) vis_resnorm
          end do
         end do
        end do
        read(IN_FILE_UNIT, *) 
      elseif (read_data_format == 'BINARY') then
        read(IN_FILE_UNIT) !'SCALARS Resnorm DOUBLE'
        read(IN_FILE_UNIT) !'LOOKUP_TABLE default'
        do k = 1, kmx - 1
         do j = 1, jmx - 1
          do i = 1, imx - 1
            read(IN_FILE_UNIT) vis_resnorm
          end do
         end do
        end do
        read(IN_FILE_UNIT)
      end if

    end subroutine read_resnorm

end module read_output_vtk
