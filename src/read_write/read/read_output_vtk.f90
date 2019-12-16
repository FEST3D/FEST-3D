  !< Read the restart file in the vtk format
module read_output_vtk
  !< Read the restart file in the vtk format
  !---------------------------------------------------------
  ! This module read state + other variable in output file
  !---------------------------------------------------------
#include "../../debug.h"
#include "../../error.h"
  use vartypes
  use global     , only :    IN_FILE_UNIT
  use global     , only : OUTIN_FILE_UNIT
  use global     , only : outin_file

  use global_vars, only : process_id
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
  use global_vars, only : intermittency
  use global_vars, only : dist

!  use global_vars, only : r_count
!  use global_vars, only : r_list

  use utils
!  use string

  implicit none
  private
  integer :: imx, jmx, kmx
  integer :: i,j,k

  public :: read_file

  contains

    subroutine read_file(control, dims)
      !< Read all the variable for the vtk restart file
      implicit none
      type(controltype), intent(in) :: control
      type(extent), intent(in) :: dims
      integer :: n
      imx = dims%imx
      jmx = dims%jmx
      kmx = dims%kmx

      call read_header()
      call read_grid()
      do n = 1,control%r_count

        select case (trim(control%r_list(n)))
        
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

          case('tv')
            call read_scalar(tv, 'tv', -2)

          case('tgm')
            call read_scalar(tgm, 'tgm', -2)

          case('Intermittecny')
            call read_scalar(intermittency, 'Intermittecny', -2)

          case('do not read')
            call skip_scalar()

          case Default
            Print*, "read error: list var : "//trim(control%r_list(n))

        end select
      end do

    end subroutine read_file

    subroutine read_header()
      !< Skip read the header in the vtk file
      implicit none

      DebugCall('read_output_vtk: read_header')
      read(IN_FILE_UNIT, *) !'# vtk DataFile Version 3.1'
      read(IN_FILE_UNIT, *) !'cfd-iitm output'   ! comment line
      read(IN_FILE_UNIT, *) !trim(read_data_format)
      read(IN_FILE_UNIT, *) !'DATASET STRUCTURED_GRID'
      !read(IN_FILE_UNIT, *)


    end subroutine read_header

    subroutine read_grid()
      !< Skip the grid read in the restart file
      implicit none

      ! read grid point coordinates
      DebugCall('read_output_vtk: read_grid')
      read(IN_FILE_UNIT, * ) !'DIMENSIONS ', imx, ' ', jmx, ' ', kmx
      read(IN_FILE_UNIT, * ) !'POINTS ', imx*jmx*kmx, ' DOUBLE'
      do k = 1, kmx
       do j = 1, jmx
        do i = 1, imx
          read(IN_FILE_UNIT, *) !point(i, j, k)%x, ' ', point(i, j, k)%y, ' ', point(i, j, k)%z
        end do
       end do
      end do
      read(IN_FILE_UNIT, *) 

    end subroutine read_grid

    subroutine read_velocity()
      !< Read velocity vector from the vtk file
      implicit none

      DebugCall('read_output_vtk: read_velocity')
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

    end subroutine read_velocity

    subroutine read_scalar(var, name,  index)
      !< Read scalar from the vtk file
      implicit none
      integer, intent(in) :: index
      real, dimension(index:imx-index,index:jmx-index,index:kmx-index), intent(out) :: var
      character(len=*), intent(in) :: name

      DebugCall('read_output_vtk'//trim(name))
      read(IN_FILE_UNIT, *) !'SCALARS '//trim(name)//' FLOAT'
      read(IN_FILE_UNIT, *) !'LOOKUP_TABLE default'
      do k = 1, kmx - 1
       do j = 1, jmx - 1
        do i = 1, imx - 1
          read(IN_FILE_UNIT, *) var(i, j, k)
        end do
       end do
      end do
      read(IN_FILE_UNIT, *)

    end subroutine read_scalar

    subroutine skip_scalar()
      !< Skip read scalar from the vtk file
      implicit none

      DebugCall('read_output_vtk: skip_scalar')
      read(IN_FILE_UNIT, *) !'SCALARS '//trim(name)//' FLOAT'
      read(IN_FILE_UNIT, *) !'LOOKUP_TABLE default'
      do k = 1, kmx - 1
       do j = 1, jmx - 1
        do i = 1, imx - 1
          read(IN_FILE_UNIT, *)
        end do
       end do
      end do
      read(IN_FILE_UNIT, *)

    end subroutine skip_scalar

end module read_output_vtk
