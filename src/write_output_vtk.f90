module write_output_vtk
  !---------------------------------------------------------
  ! This module write state + other variable in output file
  !---------------------------------------------------------
  use global     , only : OUT_FILE_UNIT
  use global_vars, only : outfile

  use global_vars, only : write_data_format
  use global_vars, only : write_file_format
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

  use global_vars, only : turbulence

  use utils
  use string

  implicit none
  private
  integer :: i,j,k
  real    :: speed_inf
  character(len=8) :: file_format
  character(len=16) :: data_format

  public :: write_file

  contains

    subroutine write_file()
      implicit none

      call setup_file
      call open_file(outfile)
      call write_header()
      call write_grid()
      call write_velocity()
      call write_density()
      call write_pressure()
      call write_turbulent_variables()
      call write_resnorm()
      call close_file(outfile)

    end subroutine write_file

    subroutine write_turbulent_variables()
      implicit none

      select case (turbulence)
        
        case ('none')
          !do nothing
          continue
        case ('sst')
          call write_TKE()
          call write_omega()
        case DEFAULT
          call dmsg(5, 'write_output_vtk',' write_turbulent_variables',&
              'ERROR: Turbulence model not recongnised')
          STOP
      end select

    end subroutine write_turbulent_variables

    subroutine setup_file()
      implicit none
      call dmsg(1, 'write_output_vtk', 'setup_file')
      if (write_file_format == "vtk") then
        file_format = ".vtk"
      elseif (write_file_format == "tecplot") then
        file_format = ".dat"
      else
        print*, "File format not recoganised. Accepted formats are"
        print*, "'vtk' and 'tecplot' "
      end if

      if (write_data_format == "ASCII") then
        data_format = "formatted"
      elseif (write_data_format == "BINARY") then
        data_format = "unformatted"
      else
        print*, "Data format not recoganised. Accepted formats are"
        print*, "'ASCII' and 'BINARY' "
      end if

    end subroutine setup_file

    subroutine open_file(filename)
      implicit none
      character(len=*), intent(in) :: filename 
      call dmsg(1, 'write_output_vtk', 'open_file')

      open(OUT_FILE_UNIT, file=trim(filename)//trim(file_format) + '.part', form=trim(data_format))

    end subroutine open_file

    subroutine close_file(filename)
      implicit none

      character(len=*), intent(in) :: filename 
      call dmsg(1, 'write_output_vtk', 'close_file')
      call rename(trim(filename)//trim(file_format) + '.part', trim(filename)//trim(file_format))
      close(OUT_FILE_UNIT)
    end subroutine close_file


    subroutine write_header()
      implicit none

      call dmsg(1, 'write_output_vtk', 'write_header')

      if (Write_data_format == "ASCII") then
        write(OUT_FILE_UNIT, fmt='(a)') '# vtk DataFile Version 3.1'
        write(OUT_FILE_UNIT, '(a)') 'cfd-iitm output'   ! comment line
        write(OUT_FILE_UNIT, '(a)') trim(Write_data_format)
        write(OUT_FILE_UNIT, '(a)') 'DATASET STRUCTURED_GRID'
        write(OUT_FILE_UNIT, *)
      elseif (write_data_format == 'BINARY') then
        write(OUT_FILE_UNIT) '# vtk DataFile Version 3.1'
        write(OUT_FILE_UNIT) 'cfd-iitm output'
        write(OUT_FILE_UNIT) trim(Write_data_format)
        write(OUT_FILE_UNIT) 'DATASET STRUCTURED_GRID'
        write(OUT_FILE_UNIT)
      end if


    end subroutine write_header

    subroutine write_grid()
      implicit none

      ! write grid point coordinates
      call dmsg(1, 'write_output_vtk', 'write_grid')
      if (Write_data_format == "ASCII") then
        write(OUT_FILE_UNIT, fmt='(a, i0, a, i0, a, i0)') &
            'DIMENSIONS ', imx, ' ', jmx, ' ', kmx
        write(OUT_FILE_UNIT, fmt='(a, i0, a)') &
            'POINTS ', imx*jmx*kmx, ' DOUBLE'
        do k = 1, kmx
         do j = 1, jmx
          do i = 1, imx
              write(OUT_FILE_UNIT, fmt='(f0.16, a, f0.16, a, f0.16)') &
                  grid_x(i, j, k), ' ', grid_y(i, j, k), ' ', grid_z(i, j, k)
          end do
         end do
        end do
        write(OUT_FILE_UNIT, *) 
      elseif (write_data_format == 'BINARY') then
        write(OUT_FILE_UNIT) &
            'DIMENSIONS ', imx, ' ', jmx, ' ', kmx
        write(OUT_FILE_UNIT) &
            'POINTS ', imx*jmx*kmx, ' DOUBLE'
        do k = 1, kmx
         do j = 1, jmx
          do i = 1, imx
              write(OUT_FILE_UNIT) &
                  grid_x(i, j, k), ' ', grid_y(i, j, k), ' ', grid_z(i, j, k)
          end do
         end do
        end do
        write(OUT_FILE_UNIT) 
      end if

    end subroutine write_grid

    subroutine write_velocity()
      implicit none
      call dmsg(1, 'write_output_vtk', 'write_velocity')

        ! Cell data
        ! Writing Velocity
      if (Write_data_format == "ASCII") then
        write(OUT_FILE_UNIT, fmt='(a, i0)') 'CELL_DATA ', (imx-1)*(jmx-1)*(kmx-1)
        write(OUT_FILE_UNIT, '(a)') 'VECTORS Velocity FLOAT'
        do k = 1, kmx - 1
         do j = 1, jmx - 1
          do i = 1, imx - 1
            write(OUT_FILE_UNIT, fmt='(f0.16, a, f0.16, a, f0.16)') &
                x_speed(i, j, k), ' ', y_speed(i, j, k), ' ', z_speed(i, j, k)
          end do
         end do
        end do
        write(OUT_FILE_UNIT, *) 
      elseif (write_data_format == 'BINARY') then
        write(OUT_FILE_UNIT) 'CELL_DATA ', (imx-1)*(jmx-1)*(kmx-1)
        write(OUT_FILE_UNIT) 'VECTORS Velocity DOUBLE'
        do k = 1, kmx - 1
         do j = 1, jmx - 1
          do i = 1, imx - 1
            write(OUT_FILE_UNIT) &
                x_speed(i, j, k), ' ', y_speed(i, j, k), ' ', z_speed(i, j, k)
          end do
         end do
        end do
        write(OUT_FILE_UNIT) 
      end if

    end subroutine write_velocity

    subroutine write_density()
      implicit none

      call dmsg(1, 'write_output_vtk', 'write_density')
      ! Writing Density
      if (Write_data_format == "ASCII") then
        write(OUT_FILE_UNIT, '(a)') 'SCALARS Density FLOAT'
        write(OUT_FILE_UNIT, '(a)') 'LOOKUP_TABLE default'
        do k = 1, kmx - 1
         do j = 1, jmx - 1
          do i = 1, imx - 1
            write(OUT_FILE_UNIT, fmt='(f0.16)') density(i, j, k)
          end do
         end do
        end do
        write(OUT_FILE_UNIT, *) 
      elseif (write_data_format == 'BINARY') then
        write(OUT_FILE_UNIT) 'SCALARS Density DOUBLE'
        write(OUT_FILE_UNIT) 'LOOKUP_TABLE default'
        do k = 1, kmx - 1
         do j = 1, jmx - 1
          do i = 1, imx - 1
            write(OUT_FILE_UNIT) density(i, j, k)
          end do
         end do
        end do
        write(OUT_FILE_UNIT)
      end if

    end subroutine write_density

    subroutine write_pressure()
      implicit none

      call dmsg(1, 'write_output_vtk', 'write_pressure')
      ! Writing Pressure
      if (Write_data_format == "ASCII") then
        write(OUT_FILE_UNIT, '(a)') 'SCALARS Pressure FLOAT'
        write(OUT_FILE_UNIT, '(a)') 'LOOKUP_TABLE default'
        do k = 1, kmx - 1
         do j = 1, jmx - 1
          do i = 1, imx - 1
            write(OUT_FILE_UNIT, fmt='(f0.16)') pressure(i, j, k)
          end do
         end do
        end do
        write(OUT_FILE_UNIT, *) 
      elseif (write_data_format == 'BINARY') then
        write(OUT_FILE_UNIT) 'SCALARS Pressure DOUBLE'
        write(OUT_FILE_UNIT) 'LOOKUP_TABLE default'
        do k = 1, kmx - 1
         do j = 1, jmx - 1
          do i = 1, imx - 1
            write(OUT_FILE_UNIT) pressure(i, j, k)
          end do
         end do
        end do
        write(OUT_FILE_UNIT) 
      end if

    end subroutine write_pressure

    subroutine write_TKE()
      implicit none

      call dmsg(1, 'write_output_vtk', 'write_TKE')
      ! Writing Pressure
      if (Write_data_format == "ASCII") then
        write(OUT_FILE_UNIT, '(a)') 'SCALARS k FLOAT'
        write(OUT_FILE_UNIT, '(a)') 'LOOKUP_TABLE default'
        do k = 1, kmx - 1
         do j = 1, jmx - 1
          do i = 1, imx - 1
            write(OUT_FILE_UNIT, fmt='(f0.16)') tk(i, j, k)
          end do
         end do
        end do
        write(OUT_FILE_UNIT, *) 
      elseif (write_data_format == 'BINARY') then
        write(OUT_FILE_UNIT) 'SCALARS k  DOUBLE'
        write(OUT_FILE_UNIT) 'LOOKUP_TABLE default'
        do k = 1, kmx - 1
         do j = 1, jmx - 1
          do i = 1, imx - 1
            write(OUT_FILE_UNIT) tk(i, j, k)
          end do
         end do
        end do
        write(OUT_FILE_UNIT) 
      end if

    end subroutine write_TKE

    subroutine write_omega()
      implicit none

      call dmsg(1, 'write_output_vtk', 'write_omega')
      ! Writing Pressure
      if (Write_data_format == "ASCII") then
        write(OUT_FILE_UNIT, '(a)') 'SCALARS Omega FLOAT'
        write(OUT_FILE_UNIT, '(a)') 'LOOKUP_TABLE default'
        do k = 1, kmx - 1
         do j = 1, jmx - 1
          do i = 1, imx - 1
            write(OUT_FILE_UNIT, fmt='(f0.16)') tw(i, j, k)
          end do
         end do
        end do
        write(OUT_FILE_UNIT, *) 
      elseif (write_data_format == 'BINARY') then
        write(OUT_FILE_UNIT) 'SCALARS Omega DOUBLE'
        write(OUT_FILE_UNIT) 'LOOKUP_TABLE default'
        do k = 1, kmx - 1
         do j = 1, jmx - 1
          do i = 1, imx - 1
            write(OUT_FILE_UNIT) tw(i, j, k)
          end do
         end do
        end do
        write(OUT_FILE_UNIT) 
      end if

    end subroutine write_omega

    subroutine write_dist()
      implicit none

      call dmsg(1, 'write_output_vtk', 'write_dist')
      ! Writing wall distance for each cell
      if (Write_data_format == "ASCII") then
        write(OUT_FILE_UNIT, '(a)') 'SCALARS dist FLOAT'
        write(OUT_FILE_UNIT, '(a)') 'LOOKUP_TABLE default'
        do k = 1, kmx - 1
         do j = 1, jmx - 1
          do i = 1, imx - 1
            write(OUT_FILE_UNIT, fmt='(f0.16)') dist(i, j, k)
          end do
         end do
        end do
        write(OUT_FILE_UNIT, *) 
      elseif (write_data_format == 'BINARY') then
        write(OUT_FILE_UNIT) 'SCALARS dist DOUBLE'
        write(OUT_FILE_UNIT) 'LOOKUP_TABLE default'
        do k = 1, kmx - 1
         do j = 1, jmx - 1
          do i = 1, imx - 1
            write(OUT_FILE_UNIT) dist(i, j, k)
          end do
         end do
        end do
        write(OUT_FILE_UNIT) 
      end if

    end subroutine write_dist

    subroutine write_resnorm()
      implicit none

      call dmsg(1, 'write_output_vtk', 'write_resnorm')
      ! Writing resnorm for each cell
      if (Write_data_format == "ASCII") then
        write(OUT_FILE_UNIT, '(a)') 'SCALARS Resnorm FLOAT'
        write(OUT_FILE_UNIT, '(a)') 'LOOKUP_TABLE default'
        speed_inf = sqrt(x_speed_inf**2 + y_speed_inf**2 &
                    + z_speed_inf**2)
        do k = 1, kmx - 1
         do j = 1, jmx - 1
          do i = 1, imx - 1

            energy_resnorm = (                                      &
                              (                                     &
                               energy_residue(i, j, k) /            &
                              (density_inf * speed_inf *            &
                              ((0.5 * speed_inf * speed_inf) +      &
                              (gm/(gm-1)*pressure_inf/density_inf)))&
                              ) ** 2                                &
                            )                                       

            x_mom_resnorm = (                                      &
                              (x_mom_residue(i, j, k) /            &
                              (density_inf * speed_inf ** 2)) ** 2 &
                             )                                     
              
            y_mom_resnorm = (                                      &
                              (y_mom_residue(i, j, k) /            &
                              (density_inf * speed_inf ** 2)) ** 2 &
                             )                                     
              
            z_mom_resnorm = (                                      &
                              (z_mom_residue(i, j, k) /            &
                              (density_inf * speed_inf ** 2)) ** 2 &
                             )                                     
            cont_resnorm =(                                       &
                             (mass_residue(i, j, k) /             &
                             (density_inf * speed_inf)) ** 2      &
                            )                                     
            vis_resnorm =    sqrt(                    &
                                    cont_resnorm    + &
                                    x_mom_resnorm   + &
                                    y_mom_resnorm   + &
                                    z_mom_resnorm   + &
                                    energy_resnorm    &
                                 )
                
            write(OUT_FILE_UNIT, fmt='(f0.16)') vis_resnorm
          end do
         end do
        end do
        write(OUT_FILE_UNIT, *) 
      elseif (write_data_format == 'BINARY') then
        write(OUT_FILE_UNIT) 'SCALARS Resnorm DOUBLE'
        write(OUT_FILE_UNIT) 'LOOKUP_TABLE default'
        speed_inf = sqrt(x_speed_inf**2 + y_speed_inf**2 &
                    + z_speed_inf**2)
        do k = 1, kmx - 1
         do j = 1, jmx - 1
          do i = 1, imx - 1

            energy_resnorm = (                                      &
                              (                                     &
                               energy_residue(i, j, k) /            &
                              (density_inf * speed_inf *            &
                              ((0.5 * speed_inf * speed_inf) +      &
                              (gm/(gm-1)*pressure_inf/density_inf)))&
                              ) ** 2                                &
                            )                                       

            x_mom_resnorm = (                                      &
                              (x_mom_residue(i, j, k) /            &
                              (density_inf * speed_inf ** 2)) ** 2 &
                             )                                     
              
            y_mom_resnorm = (                                      &
                              (y_mom_residue(i, j, k) /            &
                              (density_inf * speed_inf ** 2)) ** 2 &
                             )                                     
              
            z_mom_resnorm = (                                      &
                              (z_mom_residue(i, j, k) /            &
                              (density_inf * speed_inf ** 2)) ** 2 &
                             )                                     
            cont_resnorm =(                                       &
                             (mass_residue(i, j, k) /             &
                             (density_inf * speed_inf)) ** 2      &
                            )                                     
            vis_resnorm =    sqrt(                    &
                                    cont_resnorm    + &
                                    x_mom_resnorm   + &
                                    y_mom_resnorm   + &
                                    z_mom_resnorm   + &
                                    energy_resnorm    &
                                 )
                
            write(OUT_FILE_UNIT) vis_resnorm
          end do
         end do
        end do
        write(OUT_FILE_UNIT)
      end if

    end subroutine write_resnorm

end module write_output_vtk
