module read_output_vtk
  
  !---------------------------------------------------------
  ! This module read state + other variable in output file
  !---------------------------------------------------------
  use global     , only :      IN_FILE_UNIT
  use global     , only : RESTART_FILE_UNIT
  use global_vars, only :      infile
  use global_vars, only : restartfile

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


  use utils
  use string

  implicit none
  private
  integer :: i,j,k
  real    :: speed_inf
  character(len=8) :: file_format
  character(len=16) :: data_format
  character(len=16) :: read_flow_type

  public :: read_file

  contains

    subroutine read_file()
      implicit none
      call setup_file
      call open_file(infile)
      call read_restart_file()
      call read_header()
      call read_grid()
      call read_velocity()
      call read_density()
      call read_pressure()
      call read_viscosity()
      call read_turbulence_variables()
 !     call read_resnorm()
      call close_file()
    end subroutine read_file

    subroutine read_turbulence_variables()
      implicit none

      select case (read_flow_type)

        case ('viscous')
          !do nothing
          !restart turbulent varibale with infinity condition !todo
          continue
        case ('sst')
          call read_TKE()
          call read_omega()
        case DEFAULT
          call dmsg(5,'read_output_vtk', 'read_turbulence_variables',&
            'ERROR: Read flow type not recognised')
          STOP
      end select

    end subroutine read_turbulence_variables

    subroutine read_viscosity()
      implicit none

      call read_mu()

      if (turbulence/='none') then
        call read_mu_t()
      end if

    end subroutine read_viscosity

    subroutine setup_file()
      implicit none
      call dmsg(1, 'read_output_vtk', 'setup_file')
      if (read_file_format == "vtk") then
        file_format = ".vtk"
      elseif (read_file_format == "tecplot") then
        file_format = ".dat"
      else
        print*, "File format not recoganised. Accepted formats are"
        print*, "'vtk' and 'tecplot' "
      end if

      if (read_data_format == "ASCII") then
        data_format = "formatted"
      elseif (read_data_format == "BINARY") then
        data_format = "unformatted"
      else
        print*, "Data format not recoganised. Accepted formats are"
        print*, "'ASCII' and 'BINARY' "
      end if

      !write(infile,'(a,i4.4,a,i2.2)') 'time_directories/',start_from,'process_',process_id

    end subroutine setup_file

    subroutine open_file(filename)
      implicit none
      character(len=*), intent(in) :: filename 
      call dmsg(1, 'read_output_vtk', 'open_file')

      write(restartfile, '(A,I4.4,A,I2.2)') 'time_directories/',start_from,&
                          '/restart/process_', process_id
      open(IN_FILE_UNIT, file=trim(filename)//trim(file_format))!, form=trim(data_format))
      open(RESTART_FILE_UNIT, file=restartfile, status='old')

    end subroutine open_file

    subroutine close_file()
      implicit none

      call dmsg(1, 'read_output_vtk', 'close_files')
      close(IN_FILE_UNIT)
      close(RESTART_FILE_UNIT)

    end subroutine close_file


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

    subroutine read_restart_file()
      implicit none
      read(RESTART_FILE_UNIT, *) read_flow_type

      read(RESTART_FILE_UNIT, *)        resnorm_0
      read(RESTART_FILE_UNIT, *)    vis_resnorm_0
      read(RESTART_FILE_UNIT, *)   turb_resnorm_0
      read(RESTART_FILE_UNIT, *)   cont_resnorm_0
      read(RESTART_FILE_UNIT, *)  x_mom_resnorm_0
      read(RESTART_FILE_UNIT, *)  y_mom_resnorm_0
      read(RESTART_FILE_UNIT, *)  z_mom_resnorm_0
      read(RESTART_FILE_UNIT, *) energy_resnorm_0
      read(RESTART_FILE_UNIT, *)    TKE_resnorm_0
      read(RESTART_FILE_UNIT, *)  omega_resnorm_0
    end subroutine read_restart_file

end module read_output_vtk
