module write_output_vtk
  !---------------------------------------------------------
  ! This module write state + other variable in output file
  !---------------------------------------------------------
  use global     , only : OUT_FILE_UNIT
  use global     , only : OUTIN_FILE_UNIT
  use global     , only : outin_file

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
  use global_vars, only : mu 
  use global_vars, only : mu_t 
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
  use global_vars, only : mu_ref
  use global_vars, only : current_iter
  use global_vars, only : max_iters
  use global_vars, only : n_write
  use global_vars, only : rw_list

  use global_sst , only : sst_F1
  use global_vars, only : gradu_x
  use global_vars, only : gradu_y
  use global_vars, only : gradu_z
  use global_vars, only : gradv_x
  use global_vars, only : gradv_y
  use global_vars, only : gradv_z
  use global_vars, only : gradw_x
  use global_vars, only : gradw_y
  use global_vars, only : gradw_z
  use global_vars, only : gradT_x
  use global_vars, only : gradT_y
  use global_vars, only : gradT_z
  use global_vars, only : gradtk_x
  use global_vars, only : gradtk_y
  use global_vars, only : gradtk_z
  use global_vars, only : gradtw_x
  use global_vars, only : gradtw_y
  use global_vars, only : gradtw_z

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
      integer :: n
      character(len=*), parameter :: err="Write error: Asked to write non-existing variable- "

      call write_header()
      call write_grid()

      do n = 1,n_write

        select case (trim(rw_list(n)))
        
          case('U')
            call write_velocity()

          case('Density')
            call write_scalar(density ,"Density")
          
          case('Pressure')
            call write_scalar(pressure ,"Pressure")
            
          case('Mu')
            if (mu_ref/=0.0) then
              call write_scalar(mu ,"Mu")
            else
              print*, err//trim(rw_list(n))
            end if
            
          case('Mu_t')
            if (turbulence/='none') then
              call write_scalar(mu_t, "Mu_t")
            else
              print*, err//trim(rw_list(n))
            end if
            
          case('TKE')
            if(turbulence=="sst")then
            call write_scalar(tk, "TKE")
            else
              print*, err//trim(rw_list(n))
            end if

          case('Omega')
            if(turbulence=="sst") then
            call write_scalar(tw, "Omega")
            else
              print*, err//trim(rw_list(n))
            end if

          case('Wall_distance')
            if(turbulence/="none") then
            call write_scalar(dist, "dist")
            else
              print*, err//trim(rw_list(n))
            end if

          case('Resnorm')
            call write_resnorm()

          case('F1')
            call write_scalar(sst_F1 ,"F1")

          case('Dudx')
            call write_scalar(gradu_x ,"dudx ")

          case('Dudy')
            call write_scalar(gradu_y ,"dudy ")

          case('Dudz')
            call write_scalar(gradu_z ,"dudz ")

          case('Dvdx')
            call write_scalar(gradv_x ,"dvdx ")

          case('Dvdy')
            call write_scalar(gradv_y ,"dvdy ")

          case('Dvdz')
            call write_scalar(gradv_z ,"dvdz ")

          case('Dwdx')
            call write_scalar(gradw_x ,"dwdx ")

          case('Dwdy')
            call write_scalar(gradw_y ,"dwdy ")

          case('Dwdz')
            call write_scalar(gradw_z ,"dwdz ")

          case('DTdx')
            call write_scalar(gradT_x ,"dTdx ")

          case('DTdy')
            call write_scalar(gradT_y ,"dTdy ")

          case('DTdz')
            call write_scalar(gradT_z ,"dTdz ")

          case('Dtkdx')
            call write_scalar(gradtk_x,"dtkdx")

          case('Dtkdy')
            call write_scalar(gradtk_y,"dtkdy")

          case('Dtkdz')
            call write_scalar(gradtk_z,"dtkdz")

          case('Dtwdx')
            call write_scalar(gradtw_x,"dtwdx")

          case('Dtwdy')
            call write_scalar(gradtw_y,"dtwdy")

          case('Dtwdz')
            call write_scalar(gradtw_z,"dtwdz")

          case Default
            print*, err//trim(rw_list(n))//" to file"

        end select
      end do


    end subroutine write_file


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


    subroutine write_scalar(var, name)
      implicit none
      real, dimension(:,:,:), intent(in) :: var
      character(len=*),       intent(in):: name
      character                          :: newline=achar(10)
      character(len=32)                  :: line

      call dmsg(1, 'write_output_vtk', trim(name))

      if (Write_data_format == "ASCII") then
        write(OUT_FILE_UNIT, '(a)') 'SCALARS '//trim(name)//' FLOAT'
        write(OUT_FILE_UNIT, '(a)') 'LOOKUP_TABLE default'
        do k = 1, kmx - 1
         do j = 1, jmx - 1
          do i = 1, imx - 1
            write(OUT_FILE_UNIT, fmt='(f0.16)') var(i, j, k)
          end do
         end do
        end do
        write(OUT_FILE_UNIT, *)
      elseif (write_data_format == 'BINARY') then
        write(OUT_FILE_UNIT) 'SCALARS '//trim(name)//' FLOAT'//newline
        write(OUT_FILE_UNIT) 'LOOKUP_TABLE default'//newline
        do k = 1, kmx - 1
         do j = 1, jmx - 1
          do i = 1, imx - 1
            write(line, "(f0.16)") var(i,j,k)
            write(OUT_FILE_UNIT) trim(line)//newline
          end do
         end do
        end do
        write(OUT_FILE_UNIT)  newline
      end if

    end subroutine write_scalar


end module write_output_vtk
