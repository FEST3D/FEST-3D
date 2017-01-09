module read_output_vtk
  
  !---------------------------------------------------------
  ! This module read state + other variable in output file
  !---------------------------------------------------------
  use global     , only : IN_FILE_UNIT
  use global_vars, only : intfile

  use global_vars, only : read_data_format
  use global_vars, only : read_file_format
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
  use global_vars, only : density_inf
  use global_vars, only : x_speed_inf
  use global_vars, only : y_speed_inf
  use global_vars, only : z_speed_inf
  use global_vars, only : pressure_inf 
  use global_vars, only : gm
  use global_vars, only : dist
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

  use utils
  use string

  implicit none
  private
  integer :: i,j,k
  real    :: speed_inf
  character(len=8) :: file_format
  character(len=16) :: data_format

  public :: read_file

  contains

    subroutine read_file()
      implicit none
      call setup_file
      call open_file(infile)
      call read_header()
      call read_grid()
      call read_velocity()
      call read_density()
      call read_pressure()
 !     call read_resnorm()
      call close_file(infile)
    end subroutine read_file

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

    end subroutine setup_file

    subroutine open_file(filename)
      implicit none
      character(len=*), intent(in) :: filename 
      call dmsg(1, 'read_output_vtk', 'open_file')

      open(IN_FILE_UNIT, file=trim(filename)//trim(file_format), form=trim(data_format))

    end subroutine open_file

    subroutine close_file(filename)
      implicit none

      character(len=*), intent(in) :: filename 
      call dmsg(1, 'read_output_vtk', 'close_file')
      close(IN_FILE_UNIT)
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
