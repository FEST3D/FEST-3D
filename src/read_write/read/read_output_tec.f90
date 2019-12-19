  !< Read the restart file in the tecplot format
module read_output_tec
  !< Read the restart file in the tecplot format
  !---------------------------------------------------------
  ! This module read state + other variable in output file
  !---------------------------------------------------------
#include "../../debug.h"
#include "../../error.h"
  use vartypes
!  use global     , only : IN_FILE_UNIT
!  use global     , only : OUTIN_FILE_UNIT
!  use global     , only : outin_file
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
  use global_vars, only : mu 
  use global_vars, only : mu_t 
  use global_vars, only : dist
  use global_vars, only : intermittency

!  use global_vars, only : r_count
!  use global_vars, only : r_list

  use global_vars, only : process_id

  use utils
!  use string

  implicit none
  private
  integer :: IN_FILE_UNIT
  integer :: imx, jmx, kmx
  integer :: i,j,k
  public :: read_file

  contains

    subroutine read_file(file_handler, control, dims)
      !< Read all the variable for the tecplot restart file
      implicit none
      integer, intent(in) :: file_handler
      type(controltype), intent(in) :: control
      type(extent), intent(in) :: dims
      integer :: n

      IN_FILE_UNIT = file_handler
      imx = dims%imx
      jmx = dims%jmx
      kmx = dims%kmx

      call read_header(control)
      call read_grid()

      do n = 1,control%r_count

        select case (trim(control%r_list(n)))
        
          case('Velocity')
            call read_scalar(x_speed, "u", -2)
            call read_scalar(y_speed, "v", -2)
            call read_scalar(z_speed, "w", -2)

          case('Density')
            call read_scalar(density, "Density", -2)
          
          case('Pressure')
            call read_scalar(pressure, "Pressure", -2)
            
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

          case('Intermittency')
            call read_scalar(intermittency, 'Intermittency', -2)

          case('do not read')
            call skip_scalar()

          case Default
            Print*, "read error: list var : "//trim(control%r_list(n))

        end select
      end do

    end subroutine read_file


    subroutine read_header(control)
      !< Skip read the header in the tecplot file
      implicit none
      type(controltype), intent(in) :: control
      integer :: n

      DebugCall('read_output_tec: read_header')
      read(IN_FILE_UNIT, *) !"variables = x y z "

      do n = 1,control%r_count
        read(IN_FILE_UNIT, *) !trim(w_list(n))
      end do

      read(IN_FILE_UNIT, *) ! "zone T=block" ...
      read(IN_FILE_UNIT, *) !"Varlocation=([1-3]=Nodal)"
      read(IN_FILE_UNIT, *) !"Varlocation=([4-",total,"]=CELLCENTERED)"
      read(IN_FILE_UNIT, *) !"STRANDID"
      read(IN_FILE_UNIT, *) !"SolutionTime"

    end subroutine read_header


    subroutine read_grid()
      !< Skip the grid read in the restart file
      implicit none
      real :: dummy

      ! read grid point coordinates
      DebugCall('read_output_tec: read_grid')
      read(IN_FILE_UNIT, *) (((dummy,i=1,imx), j=1,jmx), k=1,kmx)
      read(IN_FILE_UNIT, *) (((dummy,i=1,imx), j=1,jmx), k=1,kmx)
      read(IN_FILE_UNIT, *) (((dummy,i=1,imx), j=1,jmx), k=1,kmx)

    end subroutine read_grid

    subroutine read_scalar(var, name, index)
      !< Read scalar from the tecplot file
      implicit none
      integer, intent(in) :: index
      real, dimension(index:imx-index,index:jmx-index,index:kmx-index), intent(out) :: var
      character(len=*),       intent(in):: name

      DebugCall('read_output_tec'//trim(name))
      read(IN_FILE_UNIT, *) (((var(i, j, k),i=1,imx-1), j=1,jmx-1), k=1,kmx-1)

    end subroutine read_scalar

    subroutine skip_scalar()
      !< Skip read scalar from the tecplot file
      implicit none
      real :: dummy

      DebugCall('read_output_tec: skip_scalar')
      read(IN_FILE_UNIT, *) (((dummy ,i=1,imx-1), j=1,jmx-1), k=1,kmx-1)

    end subroutine skip_scalar

end module read_output_tec
