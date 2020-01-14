  !< Read the restart file in the tecplot format
module read_output_tec
  !< Read the restart file in the tecplot format
  !---------------------------------------------------------
  ! This module read state + other variable in output file
  !---------------------------------------------------------
#include "../../debug.h"
#include "../../error.h"
  use vartypes
  use utils

  implicit none
  private
  integer :: IN_FILE_UNIT
  integer :: imx, jmx, kmx
  integer :: i,j,k
  real(wp), dimension(:, :, :), pointer :: density      
   !< Rho pointer, point to slice of qp (:,:,:,1)
  real(wp), dimension(:, :, :), pointer :: x_speed      
   !< U pointer, point to slice of qp (:,:,:,2) 
  real(wp), dimension(:, :, :), pointer :: y_speed      
   !< V pointer, point to slice of qp (:,:,:,3) 
  real(wp), dimension(:, :, :), pointer :: z_speed      
   !< W pointer, point to slice of qp (:,:,:,4)
  real(wp), dimension(:, :, :), pointer :: pressure     
   !< P pointer, point to slice of qp (:,:,:,5)
  real(wp), dimension(:, :, :), pointer :: tk        
  !< TKE, point to slice of qp (:,:,:,6)
  real(wp), dimension(:, :, :), pointer :: tw        
  !< Omega, point to slice of qp (:,:,:,7)
  real(wp), dimension(:, :, :), pointer :: te        
  !< Dissipation, point to slice of qp (:,:,:,7)
  real(wp), dimension(:, :, :), pointer :: tv        
  !< SA visocity, point to slice of qp (:,:,:,6)
  real(wp), dimension(:, :, :), pointer :: tkl       
  !< KL K-KL method, point to slice of qp (:,:,:,7)
  real(wp), dimension(:, :, :), pointer :: tgm       
  !< Intermittency of LCTM2015, point to slice of qp (:,:,:,8)
  public :: read_file

  contains

    subroutine read_file(file_handler, state, control, scheme, dims)
      !< Read all the variable for the tecplot restart file
      implicit none
      integer, intent(in) :: file_handler
      type(controltype), intent(in) :: control
      type(schemetype), intent(in) :: scheme
      type(extent), intent(in) :: dims
      real(wp), dimension(-2:dims%imx+2, -2:dims%jmx+2, -2:dims%kmx+2, 1:dims%n_var), intent(inout), target :: state
      integer :: n

      IN_FILE_UNIT = file_handler
      imx = dims%imx
      jmx = dims%jmx
      kmx = dims%kmx

      density(-2:imx+2, -2:jmx+2, -2:kmx+2) => state(:, :, :, 1)
      x_speed(-2:imx+2, -2:jmx+2, -2:kmx+2) => state(:, :, :, 2)
      y_speed(-2:imx+2, -2:jmx+2, -2:kmx+2) => state(:, :, :, 3)
      z_speed(-2:imx+2, -2:jmx+2, -2:kmx+2) => state(:, :, :, 4)
      pressure(-2:imx+2, -2:jmx+2, -2:kmx+2) => state(:, :, :, 5)

      select case (trim(scheme%turbulence))
          case ("none")
              !include nothing
              continue
          
          case ("sst", "sst2003", "bsl", "des-sst", "kw")
              tk(-2:imx+2, -2:jmx+2, -2:kmx+2) => state(:, :, :, 6)
              tw(-2:imx+2, -2:jmx+2, -2:kmx+2) => state(:, :, :, 7)

          case ("kkl")
              tk(-2:imx+2, -2:jmx+2, -2:kmx+2) => state(:, :, :, 6)
              tkl(-2:imx+2, -2:jmx+2, -2:kmx+2) => state(:, :, :, 7)

          case ("sa", "saBC")
              tv(-2:imx+2, -2:jmx+2, -2:kmx+2) => state(:, :, :, 6)

          case ("ke")
              tk(-2:imx+2, -2:jmx+2, -2:kmx+2) => state(:, :, :, 6)
              te(-2:imx+2, -2:jmx+2, -2:kmx+2) => state(:, :, :, 7)

          case DEFAULT
            Fatal_error
      end select

      ! Transition modeling
      select case(trim(scheme%transition))
        case('lctm2015')
          tgm(-2:imx+2, -2:jmx+2, -2:kmx+2) => state(:, :, :, 8)

        case('bc', 'none')
          !do nothing
          continue

        case DEFAULT
          Fatal_error
      end Select
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
      real(wp) :: dummy

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
      real(wp), dimension(index:imx-index,index:jmx-index,index:kmx-index), intent(out) :: var
      character(len=*),       intent(in):: name

      DebugCall('read_output_tec'//trim(name))
      read(IN_FILE_UNIT, *) (((var(i, j, k),i=1,imx-1), j=1,jmx-1), k=1,kmx-1)

    end subroutine read_scalar

    subroutine skip_scalar()
      !< Skip read scalar from the tecplot file
      implicit none
      real(wp) :: dummy

      DebugCall('read_output_tec: skip_scalar')
      read(IN_FILE_UNIT, *) (((dummy ,i=1,imx-1), j=1,jmx-1), k=1,kmx-1)

    end subroutine skip_scalar

end module read_output_tec
