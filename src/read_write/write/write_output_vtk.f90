  !< VTK module to write the solution in the vtk format
module write_output_vtk
  !< VTK module to write the solution in the vtk format
  !---------------------------------------------------------
  ! This module write state + other variable in output file
  !---------------------------------------------------------
#include "../../debug.h"
#include "../../error.h"
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
  use global_vars, only : tkl
  use global_vars, only : tv
  use global_vars, only : tgm
  use global_vars, only : te
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
  use global_vars, only : TKE_residue
  use global_vars, only : Tv_residue
  use global_vars, only : intermittency
  use global_vars, only : ExtraVar1
  use global_vars, only : ExtraVar2
  use global_vars, only : ExtraVar3
  use global_vars, only : ExtraVar4
  use global_vars, only : ExtraVar5

  use global_vars, only : process_id
  use global_vars, only : turbulence
  use global_vars, only : mu_ref
  use global_vars, only : current_iter
  use global_vars, only : max_iters
  use global_vars, only : w_count
  use global_vars, only : w_list

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
  character                          :: newline=achar(10)

  public :: write_file

  contains

    subroutine write_file()
      !< Write the header and variables in the file "process_xx.dat".
      implicit none
      integer :: n
      character(len=*), parameter :: err="Write error: Asked to write non-existing variable- "
      DebugCall("write_file")

      call write_header()
      call write_grid()

      do n = 1,w_count

        select case (trim(w_list(n)))

          case('Velocity')
            call write_velocity()

          case('Density')
            call write_scalar(density ,"Density", -2)

          case('Pressure')
            call write_scalar(pressure ,"Pressure", -2)

          case('Mu')
            call write_scalar(mu ,"Mu", -2)

          case('Mu_t')
            call write_scalar(mu_t, "Mu_t", -2)

          case('TKE')
            call write_scalar(tk, "TKE", -2)

          case('Omega')
            call write_scalar(tw, "Omega", -2)

          case('Kl')
            call write_scalar(tkl, "Kl", -2)

          case('tv')
            call write_scalar(tv, "tv", -2)

          case('tgm')
            call write_scalar(tgm, "tgm", -2)

          case('Dissipation')
            call write_scalar(te, "Dissipation", -2)

          case('Wall_distance')
            call write_scalar(dist, "dist", -2)

          case('Resnorm')
            call write_resnorm()

          case('TKE_residue')
            call write_scalar(TKE_residue ,"TKE_residue", 1)

          case('Tv_residue')
            call write_scalar(Tv_residue ,"Tv_residue", 1)

          case('F1')
            call write_scalar(sst_F1 ,"F1", -2)

          case('Dudx')
            call write_scalar(gradu_x ,"dudx ", 0)

          case('Dudy')
            call write_scalar(gradu_y ,"dudy ", 0)

          case('Dudz')
            call write_scalar(gradu_z ,"dudz ", 0)

          case('Dvdx')
            call write_scalar(gradv_x ,"dvdx ", 0)

          case('Dvdy')
            call write_scalar(gradv_y ,"dvdy ", 0)

          case('Dvdz')
            call write_scalar(gradv_z ,"dvdz ", 0)

          case('Dwdx')
            call write_scalar(gradw_x ,"dwdx ", 0)

          case('Dwdy')
            call write_scalar(gradw_y ,"dwdy ", 0)

          case('Dwdz')
            call write_scalar(gradw_z ,"dwdz ", 0)

          case('DTdx')
            call write_scalar(gradT_x ,"dTdx ", 0)

          case('DTdy')
            call write_scalar(gradT_y ,"dTdy ", 0)

          case('DTdz')
            call write_scalar(gradT_z ,"dTdz ", 0)

          case('Dtkdx')
            call write_scalar(gradtk_x,"dtkdx", 0)

          case('Dtkdy')
            call write_scalar(gradtk_y,"dtkdy", 0)

          case('Dtkdz')
            call write_scalar(gradtk_z,"dtkdz", 0)

          case('Dtwdx')
            call write_scalar(gradtw_x,"dtwdx", 0)

          case('Dtwdy')
            call write_scalar(gradtw_y,"dtwdy", 0)

          case('Dtwdz')
            call write_scalar(gradtw_z,"dtwdz", 0)

          case('y-mom-residue')
            call write_scalar(y_mom_residue, 'Y_mom_residue', 1)

          case('Intermittency')
            call write_scalar(intermittency, "Intermittency", -2)
          
          case('extravar1')
            if(allocated(ExtraVar1))then
              call write_scalar(ExtraVar1, "ExtraVar1", -2)
            else
              Issue_warning
            end if
          
          case('extravar2')
            if(allocated(ExtraVar2))then
              call write_scalar(ExtraVar2, "ExtraVar2", -2)
            else
              Issue_warning
            end if
          
          case('extravar3')
            if(allocated(ExtraVar3))then
              call write_scalar(ExtraVar3, "ExtraVar3", -2)
            else
              Issue_warning
            end if
          
          case('extravar4')
            if(allocated(ExtraVar4))then
              call write_scalar(ExtraVar4, "ExtraVar4", -2)
            else
              Issue_warning
            end if
          
          case('extravar5')
            if(allocated(ExtraVar5))then
              call write_scalar(ExtraVar5, "ExtraVar5", -2)
            else
              Issue_warning
            end if
          
          case('do not write')
            ! do nothing
            continue

          case Default
            print*, err//trim(w_list(n))//" to file"

        end select
      end do


    end subroutine write_file


    subroutine write_header()
      !< Write the header in the output file in the tecplot format
      implicit none

      DebugCall("write_header")

      if (Write_data_format == "ASCII") then
        write(OUT_FILE_UNIT, fmt='(a)') '# vtk DataFile Version 3.1'
        write(OUT_FILE_UNIT, '(a)') 'cfd-iitm output'   ! comment line
        write(OUT_FILE_UNIT, '(a)') trim(Write_data_format)
        write(OUT_FILE_UNIT, '(a)') 'DATASET STRUCTURED_GRID'
        !write(OUT_FILE_UNIT, *)
      elseif (write_data_format == 'BINARY') then
        write(OUT_FILE_UNIT) '# vtk DataFile Version 3.1'//newline
        write(OUT_FILE_UNIT) 'cfd-iitm output'//newline
        write(OUT_FILE_UNIT) trim(Write_data_format)//newline
        write(OUT_FILE_UNIT) 'DATASET STRUCTURED_GRID'//newline
        write(OUT_FILE_UNIT) newline
      end if


    end subroutine write_header

    subroutine write_grid()
      !< Write the grid information in the output file
      implicit none

      ! write grid point coordinates
      DebugCall("write_grid")
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
                  grid_x(i, j, k), ' ', grid_y(i, j, k), ' ', grid_z(i, j, k), newline
          end do
         end do
        end do
        write(OUT_FILE_UNIT)
      end if

    end subroutine write_grid

    subroutine write_velocity()
      !< write the velocity vector in the output file
      implicit none
      DebugCall("write_velocity")

        ! Cell data
        ! Writing Velocity
      if (Write_data_format == "ASCII") then
        write(OUT_FILE_UNIT, fmt='(a, i0)') 'CELL_DATA ', (imx-1)*(jmx-1)*(kmx-1)
        write(OUT_FILE_UNIT, '(a)') 'VECTORS Velocity FLOAT'
        do k = 1, kmx - 1
         do j = 1, jmx - 1
          do i = 1, imx - 1
            write(OUT_FILE_UNIT, fmt='(ES27.16E4, a, ES27.16E4, a, ES27.16E4)') &
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
                x_speed(i, j, k), ' ', y_speed(i, j, k), ' ', z_speed(i, j, k), newline
          end do
         end do
        end do
        write(OUT_FILE_UNIT)
      end if

    end subroutine write_velocity

    subroutine write_resnorm()
      !< write the residual information in the output file
      implicit none

      DebugCall("write_resnorm")
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


    subroutine write_scalar(var, name, index)
      !< write the scalar variable in the output file
      implicit none
      integer, intent(in) :: index
      real, dimension(index:imx-index,index:jmx-index,index:kmx-index), intent(in) :: var
      character(len=*),       intent(in):: name
      character(len=128)                  :: line

      DebugCall("write_scalar: "//trim(name))

      if (Write_data_format == "ASCII") then
        write(OUT_FILE_UNIT, '(a)') 'SCALARS '//trim(name)//' FLOAT'
        write(OUT_FILE_UNIT, '(a)') 'LOOKUP_TABLE default'
        do k = 1, kmx - 1
         do j = 1, jmx - 1
          do i = 1, imx - 1
            write(OUT_FILE_UNIT, fmt='(ES25.16E4)') var(i, j, k)
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
            write(line, "(ES28.16E4)") var(i,j,k)
            write(OUT_FILE_UNIT) trim(line)//newline
          end do
         end do
        end do
        write(OUT_FILE_UNIT)  newline
      end if

    end subroutine write_scalar


end module write_output_vtk
